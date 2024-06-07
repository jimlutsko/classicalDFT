#ifndef __LUTSKO_DFT_FACTORY__
#define __LUTSKO_DFT_FACTORY__

#include "options.h"
#include "TimeStamp.h"

#ifdef USE_GRACE
#include "Grace.h"
#endif

#ifdef USE_OMP
#include <omp.h>
#endif

#include "DFT.h"
#include "Log.h"
#include "myColor.h"
#include "Minimizer.h"


static string LJ_name = "LJ";
static string WHDF_name = "WHDF";

template <class DensityType>
class DFT_Factory
{
public:

  DFT_Factory(int argc, char** argv, bool verbose = true)
  {
    add_options(argc, argv, verbose);
  }
  DFT_Factory(){}

  void add_options(int argc, char** argv, bool verbose = true)
  {
    argc_    = argc;
    argv_    = argv;
    verbose_ = verbose;
      
    options_.addOption("nCores", &nCores_);
    options_.addOption("kT", &kT_);

    options_.addOption("HSD", &hsd1_);

    options_.addOption("eps1",   &eps1_);
    options_.addOption("sigma1", &sigma1_);
    options_.addOption("rcut1",  &rcut1_);

    options_.addOption("Lx", L_);
    options_.addOption("Ly", L_+1);
    options_.addOption("Lz", L_+2);

    options_.addOption("Cell_Size",&cellsize_);
    
    options_.addOption("Dx", &dx_);
    options_.addOption("Dy", &dy_);
    options_.addOption("Dz", &dz_);

    options_.addOption("HomogeneousBoundary", &homogeneous_boundary_);
    options_.addOption("FixedBoundary", &fixed_background_);
    options_.addOption("FixedBackground", &fixed_background_);
    options_.addOption("Open_System", &fixed_background_);
    options_.addOption("BoundaryWidth",&boundary_width_);
    
    options_.addOption("MaxIterations", &maxIterations_);
    options_.addOption("Tolerence", &tol_);

    options_.addOption("DensityInputFile", &infile_);
    options_.addOption("DensityInputStream", &instream_);
    options_.addOption("DensityOutputFile", &outfile_);
    options_.addOption("DensityOutputStream", &outstream_);
      
    options_.addOption("ShowGraphics", &show_graphics_);
    options_.addOption("Include_HS", &include_hs_);
    options_.addOption("Include_Interaction", &include_interaction_);

    options_.addOption("Potential", &potential_name_);
    options_.addOption("EOS_Correction", &eos_correction_);
    options_.addOption("D_EOS", &D_EOS_);
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

  void set_density_input_file(string str) {infile_ = str;}
  void set_density_input_stream(string str) {instream_ = str;}

  void set_density_output_file(string str) {outfile_ = str;}
  
  void set_option_no_log()      { include_log_     = false;}
  void set_option_no_dynamics() { include_density_ = include_interaction_ = include_hs_ = false;}
  
  void set_L(double lx, double ly, double lz) { L_[0] = lx; L_[1] = ly; L_[2] = lz;}
  
  // Either of these works ... the first requires C++20
  //void addOption(const char  *name, auto *variable) { options_.addOption(name,variable);}
  template<typename Q> void addOption(const char *name, Q *variable) { options_.addOption(name,variable);}
  
  void initialize()
  {
    options_.read(argc_, argv_, verbose_);
    
    if(SourceInput_.empty() == false)
      options_.read(SourceInput_.c_str(), verbose_);

    if(cellsize_ > 0)
      {
	stringstream ss;
	ss << "log_" << cellsize_;
	log_file_name_ = ss.str();
      }
    
    if(include_log_ && theLog_ == NULL) theLog_ = new Log(log_file_name_.c_str(),-1,-1,NULL, -1, true, verbose_);
    if(verbose_ && theLog_ != NULL) *theLog_ << myColor::GREEN << "#=================================" << myColor::RESET << endl << "#" << endl;

    if(verbose_ && theLog_ != NULL) *theLog_ << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;

    if(verbose_ && theLog_ != NULL) options_.write(*theLog_);
    if(verbose_ && theLog_ != NULL) *theLog_ <<  myColor::GREEN << "#=================================" <<  myColor::RESET << endl;
    
    if (dy_<0) dy_ = dx_;
    if (dz_<0) dz_ = dx_;
    
    if(L_[0] < 0) L_[0] = (cellsize_ > 0 ? cellsize_ : dx_);
    if(L_[1] < 0) L_[1] = (cellsize_ > 0 ? cellsize_ : dy_);
    if(L_[2] < 0) L_[2] = (cellsize_ > 0 ? cellsize_ : dz_);

#ifdef USE_OMP
    omp_set_dynamic(0);
    omp_set_num_threads(nCores_);
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif
      
    //////////////////////////////////////
    ////// Construct DFT

    fmt_ = NULL;
    if(include_hs_)
      fmt_ = new esFMT(1,0);

    interaction1_ = NULL;
    if(include_interaction_)
      {
	if(potential_name_ == LJ_name)
	  potential1_ = new LJ(sigma1_, eps1_, rcut1_);
	else if(potential_name_ == WHDF_name)
	  potential1_ = new WHDF(sigma1_, eps1_, rcut1_);
	else throw std::runtime_error("Requested potential " + potential_name_ + " unknown to DFT_Factory");
	  
	if(hsd1_ < 0) hsd1_ = potential1_->getHSD(kT_);
      } else if(hsd1_ < 1) hsd1_ = 1;

    if(include_density_)
      {
#ifdef DENSITY_CONSTRUCTOR_DX_DY_DZ
	theDensity_ = new DensityType(dx_, dy_, dz_, L_, show_graphics_);
#else
	theDensity_ = new DensityType(dx_, L_, show_graphics_);
#endif

	theDensity_->set_boundary_width(boundary_width_);

	if(eos_correction_ == LJ_JZG_EOS)   eos_ = new LJ_JZG(kT_, rcut1_); // need to add no-shift option
	//	if(eos_correction_ == LJ_MECKE_EOS) eos_ = new LJ_Mecke(kT_, rcut1_); // need to add no-shift option
	
	double avdw = 0;
	if(potential1_) avdw = potential1_->getVDW_Parameter(kT_);
	  
	if(eos_ == NULL)	  
	  species1_   = new FMT_Species(*theDensity_,hsd1_,0.0,verbose_, 0);
	else 
	  species1_   = new FMT_Species_EOS(D_EOS_, *eos_, avdw, *theDensity_,hsd1_,0.0,0);

	    
	dft_        = new DFT(species1_);

	species1_->set_fixed_background(fixed_background_);
	species1_->set_homogeneous_boundary(homogeneous_boundary_);
      
	if(include_hs_) dft_->addHardCoreContribution(fmt_);

	if(include_interaction_)
	  {
	    interaction1_ = new Interaction_Interpolation_QF(species1_,species1_,potential1_,kT_,verbose_);
	    dft_->addInteraction(interaction1_);
	  }

	if(infile_.empty() == false)
	  theDensity_->readDensity(infile_.c_str());
	if(instream_.empty() == false)
	  {
	    ifstream in(instream_.c_str());
	    if(in.good()) in >> *theDensity_;      
	    else throw std::runtime_error("Input file stream no good  ... aborting");	  	    
	  }
      }
    /////////////////////////////////////////////////////
    // Report
    if(verbose_ && theLog_ != NULL) *theLog_ <<  myColor::GREEN << "#=================================" <<  myColor::RESET << endl;
    if(verbose_ && theLog_ != NULL) *theLog_ <<  myColor::GREEN << "#Summary:" <<  myColor::RESET << endl << endl;  
#ifdef USE_OMP
    if(verbose_ && theLog_ != NULL) *theLog_ << "\tomp max threads  : " << omp_get_max_threads() << endl;
#endif      
    if(verbose_ && theLog_ != NULL) *theLog_ << "\tFMT is           : " << (dft_ == NULL ? "empty" : dft_->get_fmt_name()) << endl;
    if(verbose_ && theLog_ != NULL) *theLog_ << "\tThe potential is : " << (dft_ == NULL ? "empty" : dft_->get_potential_name()) << endl;
    if(verbose_ && theLog_ != NULL) *theLog_ << endl;
    if(verbose_ && theLog_ != NULL) *theLog_ << "\tHSD1                        = " << hsd1_ << endl;
    if(potential1_ && interaction1_)
      {
	if(verbose_ && theLog_ != NULL) *theLog_ << "\tVDW parameter (potential)   = " << potential1_->getVDW_Parameter(kT_) << endl;      
	if(verbose_ && theLog_ != NULL) *theLog_ << "\tVDW parameter (interaction) = " << 0.5*interaction1_->getVDWParameter() << endl;
      }
    if(verbose_ && theLog_ != NULL) *theLog_ << endl;

    is_initialized_ = true;
  }

  // this resets the temperature and recreates the DFT (and Species and Field ...). 
  // This does NOT re-read the input file and it does not change the density or the log
  void set_temperature(double kT)
  {
    kT_ = kT;
    hsd1_ = 1;
      
    // Eliminate old objects:
    // deleting the dft has no consequences for any of the other objects
    if(dft_)          delete dft_;                   dft_ = NULL;
    if(fmt_)          delete fmt_;                   fmt_ = NULL;
    if(species1_)     delete species1_;         species1_ = NULL;
    if(interaction1_) delete interaction1_; interaction1_ = NULL;
    if(potential1_)   delete potential1_;     potential1_ = NULL;
      
    if(include_hs_)
      fmt_ = new esFMT(1,0);

    if(include_interaction_)
      {
	if(potential_name_ == LJ_name)
	  potential1_ = new LJ(sigma1_, eps1_, rcut1_);
	else if(potential_name_ == WHDF_name)
	  potential1_ = new WHDF(sigma1_, eps1_, rcut1_);
	else throw std::runtime_error("Requested potential " + potential_name_ + " unknown to DFT_Factory");
	  
	hsd1_ = potential1_->getHSD(kT_);
      }

    species1_ = new FMT_Species(*theDensity_,hsd1_,0.0,verbose_, 0);
    dft_      = new DFT(species1_);

    species1_->set_fixed_background(fixed_background_);
    species1_->set_homogeneous_boundary(homogeneous_boundary_);
      
    if(include_hs_) dft_->addHardCoreContribution(fmt_);

    if(include_interaction_)
      {
	interaction1_ = new Interaction_Interpolation_QF(species1_,species1_,potential1_,kT_,verbose_);
	dft_->addInteraction(interaction1_);
      }

    /////////////////////////////////////////////////////
    // Report
    if(verbose_ && theLog_ != NULL) *theLog_ <<  myColor::GREEN << "#=================================" <<  myColor::RESET << endl;
    if(verbose_ && theLog_ != NULL) *theLog_ <<  myColor::GREEN << "Temperature set to " << kT_ << " giving:" <<  myColor::RESET << endl << endl;  
    if(verbose_ && theLog_ != NULL) *theLog_ << "\tHSD1                        = " << hsd1_ << endl;
    if(potential1_ && interaction1_)
      {
	if(verbose_ && theLog_ != NULL) *theLog_ << "\tVDW parameter (potential)   = " << potential1_->getVDW_Parameter(kT_) << endl;      
	if(verbose_ && theLog_ != NULL) *theLog_ << "\tVDW parameter (interaction) = " << 0.5*interaction1_->getVDWParameter() << endl;
      }
    if(verbose_ && theLog_ != NULL) *theLog_ << endl;

    is_initialized_ = true;
  }  

  DFT& get_DFT() { check(); return *dft_;}
  EOS &get_eos() { if(eos_ == NULL) throw std::runtime_error("No EOS found"); return *eos_;}

  double get_cell_size() const { return cellsize_;}
  
  void get_thermodynamics(bool verbose_ = true)
  {
    check();
    /////////////////////////////////////////////////////
    // Thermodynamics
    if(verbose_ && theLog_ != NULL) *theLog_ <<  myColor::GREEN << "#=================================" <<  myColor::RESET << endl;
    if(verbose_ && theLog_ != NULL) *theLog_ <<  myColor::GREEN << "#Thermodynamics:" <<  myColor::RESET << endl << endl;
    
    xv_ = xl_ = xs1_ = xs2_ = -1;
    dft_->findSpinodal(1.0, 1e-4, xs1_, xs2_, 1e-8);
    dft_->findCoex(1.5, 1e-4, xv_, xl_,1e-8);				  

    if(verbose_ && theLog_ != NULL) *theLog_ << "\tkT = " << kT_ << endl;
    if(verbose_ && theLog_ != NULL) *theLog_ << "\tOmega/(V kT) = " << dft_->omega_times_beta_over_volume(xl_) << endl;  
    if(verbose_ && theLog_ != NULL) *theLog_ << "\tmu/kT        = " << dft_->mu_times_beta(xl_) << endl;

    if(verbose_ && theLog_ != NULL) *theLog_ << "\tSpinodal   : " << xs1_ << " " << xs2_ << endl;
    if(verbose_ && theLog_ != NULL) *theLog_ << "\tCoexistence: " << xv_ << " " << xl_ << endl;
    if(verbose_ && theLog_ != NULL) *theLog_ << endl;      
  }


  Log&              get_log()         { return *theLog_;}
  Options&          get_options()     { return options_;}
  DensityType&      get_density()     { return *theDensity_;}
  FMT&              get_FMT()         { return *fmt_;}
  Interaction_Base& get_interaction() { return *interaction1_;}
  Species&          get_species()     { return *species1_;}

  double get_temperature() const { return kT_;}
  double get_liq_coex_density() const { return xl_;}
  double get_vap_coex_density() const { return xv_;}
  double get_liq_spinodal_density() const { return xs2_;}
  double get_vap_spinodal_density() const { return xs1_;}
  string get_infile()           const { return infile_;}
  string get_outfile()          const { return outfile_;}
  string get_outstream()          const { return outstream_;}

  bool   get_homogeneous_boundary() const { return homogeneous_boundary_;}
  bool   get_fixed_background() const { return fixed_background_;}
  
  void check() const {    if(!is_initialized_) throw std::runtime_error("DFT factory not initialized");}

  void set_show_graphics(bool show) { show_graphics_ = show;}
  void set_log_file_name(string name) { log_file_name_ = name;}

  bool get_include_interactions() const { return include_interaction_;}
protected:
  int argc_;
  char **argv_;
  Options options_;

  Log         *theLog_     = NULL;
  Potential1  *potential1_ = NULL;
  DensityType *theDensity_ = NULL;
  Species     *species1_   = NULL;

  Interaction_Base *interaction1_ = NULL;
  
  FMT *fmt_ = NULL;
  DFT *dft_ = NULL;

  double xv_  = -1;
  double xl_  = -1;
  double xs1_ = -1;
  double xs2_ = -1;

  bool is_initialized_ = false;
  
public:
  int nCores_  = 6;     
  double L_[3] = {10,10,10};
  double dx_  = 0.1;
  double dy_  = -1; // if unspecified, we will copy dx_
  double dz_  = -1; // if unspecified, we will copy dx_

  double cellsize_ = -1;
  
  double kT_   = 1;
  
  double hsd1_ = -1;
  
  double eps1_   = 1;
  double sigma1_ = 1;
  double rcut1_  = 3;

  string infile_;
  string instream_;
  string outfile_;
  string outstream_;
  string SourceInput_;

  string eos_correction_;
  double D_EOS_;
  EOS *eos_ = NULL;
  
  bool fixed_background_ = false;
  bool homogeneous_boundary_ = false;
  
  double maxIterations_ = 1000;
  double tol_           = 1e-8;   
  
  bool show_graphics_        = true;

  bool include_hs_           = true;
  bool include_interaction_  = true;
  bool include_density_      = true;
  bool include_log_          = true;
  
  bool verbose_  = true;

  int boundary_width_        = 0;
  
  string potential_name_     = "LJ";
  string log_file_name_      = "log.dat";
};
#endif // __LUTSKO_DFT_FACTORY__
