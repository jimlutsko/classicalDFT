#ifndef __LUTSKO_DFT_FACTORY__
#define __LUTSKO_DFT_FACTORY__

#include "options.h"
#include "TimeStamp.h"

#include "Grace.h"
#include "DFT.h"
#include "Log.h"
#include "myColor.h"
#include "Minimizer.h"

template <class DensityType>
class DFT_Factory
{
 public:
 DFT_Factory(int argc, char** argv) : argc_(argc), argv_(argv)
    {
      options_.addOption("nCores", &nCores_);
      options_.addOption("kT", &kT_);

      options_.addOption("HSD", &hsd1_);

      options_.addOption("eps1",   &eps1_);
      options_.addOption("sigma1", &sigma1_);
      options_.addOption("rcut1",  &rcut1_);

      options_.addOption("Lx", L_);
      options_.addOption("Ly", L_+1);
      options_.addOption("Lz", L_+2);

      options_.addOption("Dx", &dx1_);

      options_.addOption("FixedBackground", &fixed_background_);
      
      options_.addOption("MaxIterations", &maxIterations_);
      options_.addOption("Tolerence", &tol_);
      //      options_.addOption("DensityInputFile", &SourceInput_);
      options_.addOption("DensityInputFile", &infile_);
      options_.addOption("DensityOutputFile", &outfile_);
    }

  ~DFT_Factory()
    {
      if(theLog_) delete theLog_;
      if(potential1_) delete potential1_;
      if(theDensity_) delete theDensity_;
      if(species1_) delete species1_;
      if(interaction1_) delete interaction1_;
      if(fmt_) delete fmt_;
      if(dft_) delete dft_;
    }
    

  // Either of these works ... the first requires C++20
  //void addOption(const char  *name, auto *variable) { options_.addOption(name,variable);}
  template<typename Q> void addOption(const char *name, Q *variable) { options_.addOption(name,variable);}
  
  void initialize()
    {
      options_.read(argc_, argv_);

      if(SourceInput_.empty() == false)
	{
	  cout << "Reading " << SourceInput_.c_str() << endl;
	  options_.read(SourceInput_.c_str());
	  cout << "Finished" << endl;
	}
      theLog_ = new Log("log.dat");
      *theLog_ << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;

      *theLog_ << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;

      options_.write(*theLog_);
      *theLog_ <<  myColor::GREEN << "=================================" <<  myColor::RESET << endl;

      if(L_[0] < 0) L_[0] = dx1_;
      if(L_[1] < 0) L_[1] = dx1_;
      if(L_[2] < 0) L_[2] = dx1_;        

#ifdef USE_OMP    
      omp_set_dynamic(0);
      omp_set_num_threads(nCores_);

      int fftw_init_threads();
      fftw_plan_with_nthreads(omp_get_max_threads());

      *theLog_ << "OMP initialized with " << omp_get_num_threads() << endl;
#endif
      
      //////////////////////////////////////
      ////// Construct DFT

      potential1_ = new LJ(sigma1_, eps1_, rcut1_);
      theDensity_ = new DensityType(dx1_, L_, true);
      fmt_        = new esFMT(1,0);
  
      if(hsd1_ < 0) hsd1_ = potential1_->getHSD(kT_);
      
      species1_     = new FMT_Species(*theDensity_,hsd1_,1);      
      interaction1_ = new Interaction_Interpolation_QF(species1_,species1_,potential1_,kT_);

      dft_ = new DFT(species1_);

      species1_->setFixedBackground(fixed_background_);
      
      dft_->addHardCoreContribution(fmt_);  
      dft_->addInteraction(interaction1_);

      if(infile_.empty() == false)
	theDensity_->readDensity(infile_.c_str());

      /////////////////////////////////////////////////////
      // Report
      *theLog_ <<  myColor::GREEN << "=================================" <<  myColor::RESET << endl;
      *theLog_ <<  myColor::GREEN << "Summary:" <<  myColor::RESET << endl << endl;  
      *theLog_ << "\tomp max threads  : " << omp_get_max_threads() << endl;
      *theLog_ << "\tFMT is           : " << dft_->get_fmt_name() << endl;
      *theLog_ << "\tThe potential is : " << dft_->get_potential_name() << endl;
      *theLog_ << endl;
      *theLog_ << "\tHSD1                        = " << hsd1_ << endl;
      *theLog_ << "\tVDW parameter (potential)   = " << potential1_->getVDW_Parameter(kT_) << endl;      
      *theLog_ << "\tVDW parameter (interaction) = " << 0.5*interaction1_->getVDWParameter() << endl;
      *theLog_ << endl;

      is_initialized_ = true;
    }

  DFT& get_DFT() { check(); return *dft_;}

  void get_thermodynamics()
  {
    check();
    /////////////////////////////////////////////////////
    // Thermodynamics
    *theLog_ <<  myColor::GREEN << "=================================" <<  myColor::RESET << endl;
    *theLog_ <<  myColor::GREEN << "Thermodynamics:" <<  myColor::RESET << endl << endl;
    
    xv_ = xl_ = xs1_ = xs2_ = -1;
    dft_->findSpinodal(1.0, 1e-4, xs1_, xs2_, 1e-8);
    dft_->findCoex(1.0, 1e-4, xv_, xl_,1e-8);				  

    *theLog_ << "\tOmega/(V kT) = " << dft_->omega_times_beta_over_volume(xl_) << endl;  
    *theLog_ << "\tmu/kT        = " << dft_->mu_times_beta(xl_) << endl;

    *theLog_ << "\tSpinodal   : " << xs1_ << " " << xs2_ << endl;
    *theLog_ << "\tCoexistence: " << xv_ << " " << xl_ << endl;
    *theLog_ << endl;      
  }


  Log&              get_log()         { return *theLog_;}
  Options&          get_options()     { return options_;}
  DensityType&      get_density()     { return *theDensity_;}
  FMT&              get_FMT()         { return *fmt_;}
  Interaction_Base& get_interaction() { return *interaction1_;}
  Species&          get_species()     { return *species1_;}

  double get_liq_coex_density() const { return xl_;}
  double get_vap_coex_density() const { return xv_;}
  string get_infile()           const { return infile_;}
  string get_outfile()          const { return outfile_;}

  void check() const {    if(!is_initialized_) throw std::runtime_error("DFT factory not initialized");}

 private:
  int argc_;
  char **argv_;
  Options options_;

  Log         *theLog_;
  Potential1  *potential1_;
  DensityType *theDensity_;
  Species     *species1_;

  Interaction_Base *interaction1_;
  
  FMT *fmt_;
  DFT *dft_;

  double xv_  = -1;
  double xl_  = -1;
  double xs1_ = -1;
  double xs2_ = -1;

  bool is_initialized_ = false;
  
 public:
  int nCores_  = 6;     
  double L_[3] = {10,10,10};
  double dx1_  = 0.1;

  double kT_   = 1;
  
  double hsd1_ = -1;
  
  double eps1_   = 1;
  double sigma1_ = 1;
  double rcut1_  = 3;

  string infile_;
  string outfile_;
  string SourceInput_;

  bool fixed_background_ = false;
  
  double maxIterations_ = 1000;
  double tol_           = 1e-8;   
};
#endif // __LUTSKO_DFT_FACTORY__
