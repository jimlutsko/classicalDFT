#ifndef __LUTSKO_CONJUG2__
#define __LUTSKO_CONJUG2__


#include <chrono>
using namespace std;

#include "DFT.h"

/**
  *  @brief Minimizer base class
  *
  */  

class Minimizer
{
 public:
  Minimizer(DFT *dft) : dft_(dft), forceLimit_(0.1),  bFrozenBoundary_(false)
  {
    x_.resize(dft->getNumberOfSpecies());

    for(auto &x: x_) x.resize(dft_->get_lattice().Ntot()); 
  }
  Minimizer() {}
  
  void setFrozenBoundaryFlag(bool f) {bFrozenBoundary_ = f;}
  
  virtual void initialize();

  void run(long maxSteps = -1);
  void resume(long maxSteps = -1);

  virtual double step() = 0;

  virtual void draw_after() {};  // Display something after the minimization

  virtual void reportMessage(string message){}

  
  int getCalls() const { return calls_;}
  double getF() const { return F_;}

  double getForceTerminationCriterion() const {return forceLimit_;}
  void   setForceTerminationCriterion(double v) {forceLimit_ = v;}
  
  virtual double getDF_DX();
  virtual double get_convergence_monitor() const { return dft_->get_convergence_monitor();}

  void setMinDensity(double m) { minDensity_ = m;}
  
 protected:
  DFT *dft_ = NULL;

  vector<DFT_Vec> x_;

  int calls_ = 0;
  int step_counter_ = 0;
  double F_ = 0;

  double forceLimit_;
  double f_abs_max_; // max absolute value of dF_
  bool bFrozenBoundary_;

  double vv_ = 0;

  double minDensity_ = -1;

  int image_number_ = 0;
  
  std::chrono::duration<double> elapsed_seconds_;

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & dft_;
    ar & x_;
    ar & calls_;
    ar & step_counter_;
    ar & F_;
    ar & forceLimit_;
    ar & f_abs_max_;
    ar & bFrozenBoundary_;
    ar & vv_;
    ar & minDensity_;
  }

  
};



/**
  *  @brief Minimizer using FIRE2 algorithm
  *
  */

class fireMinimizer2 : public Minimizer 
{
 public:
  fireMinimizer2(DFT *dft);
  fireMinimizer2(): Minimizer(){}

  virtual void   initialize();
  virtual double step();

  double getRMS_Force() const { return rms_force_;}

  void onlyRelaxSpecies(int j) { onlyRelax_ = j;}  
  
  void setTimeStep(double dt)    { dt_ = dt;}
  void setTimeStepMax(double dt) { dt_max_ = dt;}  
  void setAlphaStart(double a)   { alpha_start_ = a;}
  void setAlphaFac(double a)     { f_alf_ = a;}
  void setBacktrackFac(double a) { f_back_ = a;}  

  virtual double get_convergence_monitor() const { return fabs(vv_/dft_->get_lattice().getVolume());}
  
 protected:
  void SemiImplicitEuler(int begin, int end);
  
 protected:
 vector<DFT_Vec> v_;

  int onlyRelax_ = -1;

  double f_dec_    = 0.5;
  double f_inc_    = 1.1;
  double f_alf_    = 0.99;
  double f_back_   = 0.8;

  double dt_max_  = 1.0;
  double dt_      = 1.0;
  double dt_min_  = 0.0;        ///< minimum allowed timestep
  double dt_best_ = 0.0;
  int backtracks_ = 0;
  
  unsigned it_;     ///< Loop counter
  unsigned cut_;    ///< used to keep track of number of steps P >0
  int N_delay_ = 5;

  double alpha_;
  double alpha_start_ = 0.1;

  
  int N_P_positive_ = 0;
  int N_P_negative_ = 0;
  int N_P_negative_max_ = 20;
  
  double rms_force_ = 0.0; // for reporting
  double vnorm_     = 0.0; // for reporting
  
  bool initial_delay_ = true; ///< flag to allow for a warmup period
  double threshold_ = 0.0;
  double fmax_ = -1.0;

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Minimizer>(*this);
    ar & v_;
    ar & onlyRelax_;

    ar & f_dec_;
    ar & f_inc_;
    ar & f_alf_;
    ar & f_back_;

    ar & dt_min_;
    ar & dt_max_;
    ar & dt_;

    ar & it_;
    ar & cut_;

    ar & alpha_start_;
    ar & alpha_;

    ar & N_delay_;
    ar & N_P_positive_;
    ar & N_P_negative_;
    ar & N_P_negative_;

    ar & rms_force_;
    ar & vnorm_;
    ar & initial_delay_;
    ar & threshold_;
    ar & fmax_;


    boost::serialization::void_cast_register<fireMinimizer2, Minimizer>(static_cast<fireMinimizer2 *>(NULL),static_cast<Minimizer *>(NULL));
  }  

};

/**
  *  @brief Base class for a family of  finite elements ddft integrators. 
  *
  */  
//  TODO: Update
class DDFT : public Minimizer
{
 public:
 DDFT(DFT *dft, bool showGraphics = true)
   : Minimizer(dft), show_(showGraphics) , tolerence_fixed_point_(1e-4), successes_(0)
  {
    double dx = dft_->get_lattice().getDX();
    dt_ = 10*0.1*dx*dx;
    dt_ = 0.0001*dx*dx;
    dtMax_ = 1; //1*dx*dx;
  }
  ~DDFT() {}

  virtual void initialize();

  void set_tolerence_fixed_point(double e) { tolerence_fixed_point_ = e;}
  void set_max_time_step(double t) { dtMax_ = t;}

  void   setTimeStep(double dt) { dt_ = dt;}    
  double getTimeStep() const { return dt_;}

  virtual double step() = 0;

  void Display(double F, double dFmin, double dFmax, double N);



  
  //  double F_string(Density &d, double *fmax = NULL);
  //  void reverseForce(DFT_Vec *tangent);
  //  virtual double step_string(double &dt, Density &d, unsigned &time_den, bool verbose = true) = 0;

 protected:

  bool show_ = true;
  
  double dt_;
  double tolerence_fixed_point_ = 1e-4;

  // control of adaptive time step
  int successes_ = 0;
  double dtMax_  = 1;
};

/**
  *  @brief DDFT minimizer Class using integrating factor
  *
  *  @detailed This integrates the pure diffusion part of the dynamics exactly (using FFTs) and treats the rest implicitly via a Crank-Nicholson type method.
  *
  */  
class DDFT_IF : public DDFT
{
 public:
 DDFT_IF(DFT *dft, bool showGraphics = true)
   : DDFT(dft) {}

  ~DDFT_IF() {}

  virtual void initialize();
  virtual double step();

  void set_is_closed(bool val) { is_closed_ = val;}

  void determine_unstable_eigenvector(vector<DFT_FFT> &eigen_vector, double& eigen_value, double shift = 0.0) const;
  void Hessian_dot_v(vector<DFT_FFT> &eigen_vector, vector<DFT_Vec>& d2F) const;  
  
 protected:
  virtual double fftDiffusion(DFT_Vec &d1) = 0;
  virtual void   calcNonlinearTerm(const DFT_Vec &density, Species *species, bool use_R0) = 0;
  virtual void   finish_nonlinear_calc(DFT_Vec& d0, DFT_Vec& d1) = 0;
  
  void calcNonlinearTerm_intern(const DFT_Vec &density, DFT_Vec &dF, DFT_Vec &RHS1);
  //  virtual void update_forces_fixed_background(const Density &density,const DFT_Vec &d2, DFT_Vec &dF, const double D[3]);
  void A_dot_x(DFT_Vec& x, DFT_Vec& Ax, const Density &density, double D[], bool do_subtract_ideal = false) const; 
  
protected:
  vector<double> Lamx_;
  vector<double> Lamy_;
  vector<double> Lamz_;

  int Nx_ = 0;
  int Ny_ = 0;
  int Nz_ = 0;

  double dx_ = 0.0;
  double dy_ = 0.0;
  double dz_ = 0.0;

  bool is_closed_ = false;
};

/**
  *  @brief DDFT minimizer Class using integrating factor
  *
  *  @detailed This integrates the pure diffusion part of the dynamics exactly (using FFTs) and treats the rest implicitly via a Crank-Nicholson type method.
  *
  */  
/* TODO : Update*/
class DDFT_IF_Periodic : public DDFT_IF
{
 public:
 DDFT_IF_Periodic(DFT *dft, bool showGraphics = true)
   : DDFT_IF(dft,showGraphics) {}

  ~DDFT_IF_Periodic(){}

  virtual void initialize();


 protected:
  virtual double fftDiffusion(DFT_Vec &new_density);
  virtual void   calcNonlinearTerm(const DFT_Vec &density, Species *species, bool use_R0);
  virtual void   finish_nonlinear_calc(DFT_Vec& d0, DFT_Vec& d1);

  void restore_values_on_border(const Lattice &lattice, const DFT_Vec &d0, DFT_Vec &density);  

 protected:
  DFT_FFT RHS0;
  DFT_FFT RHS1;
};
	  

/**
  *  @brief DDFT minimizer Class using integrating factor. This is for fixed density at the boundaries - a model for an open system.
  *
  *  @detailed This integrates the pure diffusion part of the dynamics exactly (using FFTs) and treats the rest implicitly via a Crank-Nicholson type method.
  */  
/*TODO: Update */
class DDFT_IF_Fixed_Border : public DDFT_IF
{
 public:
  DDFT_IF_Fixed_Border(DFT *dft, double background,  double border_forces, bool showGraphics = true)
    : DDFT_IF(dft,showGraphics), background_(background), border_forces_(border_forces), sin_in_(NULL), sin_out_(NULL)
    {}
  ~DDFT_IF_Fixed_Border() {if(sin_in_) delete sin_in_; if(sin_out_) delete sin_out_; if(RHS0_sin_transform_) delete RHS0_sin_transform_; if(RHS1_sin_transform_) delete RHS1_sin_transform_;}

  virtual void initialize();
  
protected:
  virtual double fftDiffusion(DFT_Vec &d1);
  virtual void   calcNonlinearTerm(const DFT_Vec &density, Species *species, bool use_R0);  
  virtual void   finish_nonlinear_calc(DFT_Vec& d0, DFT_Vec& d1);

  void update_forces_fixed_background(const Density &density,const DFT_Vec &d2, DFT_Vec &dF, const double DD[3]);
  
  void pack_for_sin_transform(const double *x, double val);  
  void unpack_after_transform(double *x, double val);

  
 protected:
  DFT_Vec RHS0;
  DFT_Vec RHS1;

  double *RHS0_sin_transform_ = NULL;
  double *RHS1_sin_transform_ = NULL;  
  
  double background_;
  
  fftw_plan sin_plan_;
  unsigned  sin_Ntot_;
  unsigned  sin_Norm_;
  double   *sin_in_;
  double   *sin_out_;

  double border_forces_ = 0.0;
};
		    

/* Currently not working
class picardMinimizer : public Minimizer
{
 public:
 picardMinimizer(DFT &dft,  double alpha) :  Minimizer(dft), alpha_(alpha) {}

  virtual void   initialize()
  {
    step_counter_ = 0;
    calls_ = 0;
    
    try {
      F_ = dft_.calculateFreeEnergyAndDerivatives(false);   
    } catch( Eta_Too_Large_Exception &e) {
      throw e;
    }
  }

  virtual double step()
  {
    int numSpecies = dft_.getNumberOfSpecies();
    long Ntot = x_[0].size();
    double dV = dft_.getDensity(0).dV();

    vector<DFT_Vec> x_rem(numSpecies);

    for(int Jspecies = 0; Jspecies < numSpecies; Jspecies++)
      x_rem[Jspecies].set(dft_.getDensity(Jspecies).getDensity());

    bool bSuccess = false;
    do{
      try {
	for(int Jspecies = 0; Jspecies < numSpecies; Jspecies++)
	  for(long i = 0; i<Ntot; i++)
	    {
	      double x = dft_.getDensity(Jspecies).getDensity(i);
	      double df = dft_.getDF(Jspecies).get(i)/dV;	      	      
	      df -= log(x);	      
	      dft_.set_density(Jspecies, i, (1-alpha_)*x + alpha_*exp(-df));	    
	    }	
	calls_++;
	F_ = dft_.calculateFreeEnergyAndDerivatives(false);
	bSuccess = true;
      } catch (Eta_Too_Large_Exception &e) {
	alpha_ /= 2;
	stringstream s;
	s << "Backtrack .. alpha_ is now " << alpha_;
	reportMessage(s.str());
	
	for(int Jspecies = 0; Jspecies < numSpecies; Jspecies++)
	  dft_.set_density(Jspecies, x_rem[Jspecies]);	
      }
    } while(!bSuccess);
    return F_;
  }

  double getRMS_Force() const { return rms_force_;}

 protected:
  double alpha_ = 0.1;
  double rms_force_ = 0.0; // for reporting
};

*/

#endif // sentinal
