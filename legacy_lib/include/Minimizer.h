#ifndef __LUTSKO_CONJUG2__
#define __LUTSKO_CONJUG2__


#include <chrono>
using namespace std;

#include "DFT.h"

#include <armadillo> // for Arnoldi stuff

// brief Minimizer base class
class Minimizer
{
 public:
  Minimizer(DFT *dft);
  Minimizer() {}
  
  void         run(long maxSteps = -1); // set counters to zero and run = reset() + resume()
  void         resume(long maxSteps = -1); // run without restting counters 
  virtual void reset();   // reset internal state

  int    getCalls() const { return calls_;}  
  double getF()     const { return F_;}
  double getForceTerminationCriterion()    const {return forceLimit_;}
  virtual double get_convergence_monitor() const { return dft_->get_convergence_monitor();}
  
  void setMinDensity(double m) { minDensity_ = m;}
  void setForceTerminationCriterion(double v) {forceLimit_ = v;}  
  void setVerbose(bool verbose) { verbose_ = verbose;}
  void set_fixed_direction(const DFT_Vec& fixed) { fixed_direction_ = fixed;}
  
  // report activity
  virtual void   draw_after() {};  // Display something after the minimization
  virtual void   reportMessage(string message){}
  virtual void   cleanup() {} // called when minimization is finished to allow decendents to clean up user output.
  
protected:
  double get_energy_and_forces() {calls_++;  return dft_->calculateFreeEnergyAndDerivatives(false);}
  virtual double getDF_DX();
  virtual double step() = 0;
  
 protected:
  vector<DFT_Vec> x_; // independent variables
  DFT *dft_ = NULL; // calculates energy and forces

  bool verbose_ = false;
  
  int    calls_        = 0;
  int    step_counter_ = 0;
  double F_            = 0.0;

  double forceLimit_;
  double f_abs_max_; // max absolute value of dF_
  double vv_         = 0;
  double minDensity_ = -1;
  int image_number_  = 0;

  DFT_Vec fixed_direction_;
  
  mutable std::chrono::duration<double> elapsed_seconds_;

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

  virtual void reset();
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
class DDFT : public Minimizer, public Dynamical_Matrix
{
 public:
  DDFT(DFT *dft, bool showGraphics = true, bool central_differences = false);
  ~DDFT() {}

  double step();
  void set_tolerence_fixed_point(double e) { tolerence_fixed_point_ = e;}
  void set_max_time_step(double t) { dtMax_ = t;}

  void   setTimeStep(double dt) { dt_ = dt;}    
  double getTimeStep() const { return dt_;}

  double get_time()                const { return time_;}
  double get_convergence_monitor() const { return RHS_max_;}

  void Display(double F, double dFmin, double dFmax, double N);

  //  double F_string(Density &d, double *fmax = NULL);
  //  void reverseForce(DFT_Vec *tangent);
  //  virtual double step_string(double &dt, Density &d, unsigned &time_den, bool verbose = true) = 0;

  // Dynamical_Matrix interface
  virtual unsigned get_dimension(int direction) const {return dft_->get_dimension(direction);}
  virtual long     get_Nboundary()              const {return dft_->get_Nboundary();}
  virtual bool     get_next_boundary_point(int &ix, int &iy, int &iz) const {return dft_->get_next_boundary_point(ix,iy,iz);}
  virtual bool     get_next_boundary_point(long &pos) const {return dft_->get_next_boundary_point(pos);}
  virtual long     boundary_pos_2_pos(int p)    const {return dft_->boundary_pos_2_pos(p);}    
  virtual bool     is_boundary_point(long p) const {return dft_->is_boundary_point(p);}
  virtual bool     is_fixed_boundary() const {return dft_->is_fixed_boundary();}
  virtual void     matrix_dot_v_intern(const vector<DFT_FFT> &v, vector<DFT_Vec> &result, void *param) const;
  virtual bool     is_dynamic() const { return true;}
  
protected:
  void   g_dot_x(const DFT_Vec& x, DFT_Vec& gx) const;
  double get_neighbors(const DFT_Vec &x, int species, long pos, int stride,
			     double &xpx, double &xmx, double &xpy, double &xmy, double &xpz, double &xmz) const;

  void change_timestep(double dt);
  
  double apply_diffusion_propogator(DFT_Vec &d1);
  double calculate_excess_RHS(const Species *species, DFT_FFT& RHS) const;
  void   subtract_ideal_gas(const DFT_Vec &density, DFT_Vec& RHS) const;
 protected:

  bool show_graphics_       = true;
  bool central_differences_ = false;
  
  double dt_;
  double time_ = 0;
  double tolerence_fixed_point_ = 1e-4;

  // control of adaptive time step
  int successes_ = 0;
  double dtMax_  = 1;
  
 protected:
  DFT_FFT RHS0_;
  DFT_FFT RHS1_;

  double RHS_max_ = 1;

  vector<double> Lamx_;
  vector<double> Lamy_;
  vector<double> Lamz_;

  vector<double> fx_;
  vector<double> fy_;
  vector<double> fz_;
  
  int Nx_ = 0;
  int Ny_ = 0;
  int Nz_ = 0;

  double dx_ = 0.0;
  double dy_ = 0.0;
  double dz_ = 0.0;  
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
  DDFT_IF(DFT *dft, bool showGraphics = true);
  ~DDFT_IF() {}

  
  // ***********************************************************
  // The rest is the old interface which should eventually be removed
  
public:
  void set_is_closed(bool val) { is_closed_ = val;}

  double determine_unstable_eigenvector(vector<DFT_FFT> &eigen_vector, bool fixed_boundary, double shift, string Filename, bool dynamic = true, long maxSteps = 1000, double tol = 1e-8) const;
  
  bool check_eigenvectors(arma::cx_mat eigvec, arma::cx_vec eigval, double shift, bool fixed_boundary, bool dynamic = true, double tol=1e-8) const;
  bool check_factorisation(arma::cx_mat V, arma::cx_mat H, arma::cx_vec f, double shift, bool fixed_boundary, bool dynamic = true, double tol=1e-8) const;
  void extend_arnoldi_factorisation(arma::cx_mat &V, arma::cx_mat &H, arma::cx_vec &f, const int k, const int p, double shift, bool fixed_boundary, bool dynamic = true, double tol=1e-8) const;
  double determine_unstable_eigenvector_IRArnoldi(vector<DFT_FFT> &eigen_vector, bool fixed_boundary, double shift, string Filename, bool dynamic = true, int k=10, int p=15, long maxSteps = 1000, double tol = 1e-8) const;
  
  void Hessian_dot_v(const vector<DFT_FFT> &v, vector<DFT_Vec>& result, bool fixed_boundary, bool dynamic) const;  
  void Hessian_dot_v(arma::cx_vec v, arma::cx_vec& result, double shift, bool fixed_boundary, bool dynamic) const;

 protected:
  //  virtual void update_forces_fixed_background(const Density &density,const DFT_Vec &d2, DFT_Vec &dF, const double D[3]);
  void A_dot_x(const DFT_Vec& x, DFT_Vec& Ax, const Density &density, const double D[], bool do_subtract_ideal = false) const; 

protected:
  bool is_closed_ = false;
};

/**
  *  @brief DDFT minimizer Class using integrating factor
  *
  *  @detailed This integrates the pure diffusion part of the dynamics exactly (using FFTs) and treats the rest implicitly via a Crank-Nicholson type method.
  *
  */  
/* TODO : Update*/
/*
class DDFT_IF_Periodic : public DDFT_IF
{
 public:
  DDFT_IF_Periodic(DFT *dft, bool showGraphics = true);
  ~DDFT_IF_Periodic(){}

 protected:
  virtual double apply_diffusion_propogator(DFT_Vec &new_density);
  virtual void   calcNonlinearTerm(const DFT_Vec &density, Species *species, bool use_R0);

  void restore_values_on_border(const Lattice &lattice, const DFT_Vec &d0, DFT_Vec &density);  

  // Dynamical_Matrix interface
  virtual bool is_fixed_boundary() const{ return false;}  
  
 protected:
  DFT_FFT RHS0;
  DFT_FFT RHS1;
};
*/	  

/**
  *  @brief DDFT minimizer Class using integrating factor. This is for fixed density at the boundaries - a model for an open system.
  *
  *  @detailed This integrates the pure diffusion part of the dynamics exactly (using FFTs) and treats the rest implicitly via a Crank-Nicholson type method.
  */  
/*TODO: Update */
/*
class DDFT_IF_Fixed_Border : public DDFT_IF
{
 public:
  DDFT_IF_Fixed_Border(DFT *dft,  bool showGraphics = true);
  ~DDFT_IF_Fixed_Border();

protected:
  virtual double apply_diffusion_propogator(DFT_Vec &d1);
  virtual void   calcNonlinearTerm(const DFT_Vec &density, Species *species, bool use_R0);  

  //  void update_forces_fixed_background(const Density &density,const DFT_Vec &d2, DFT_Vec &dF, const double DD[3]);
  
  void pack_for_sin_transform(const double *x);  
  void unpack_after_transform(double *x);

  // Dynamical_Matrix interface
  virtual bool is_fixed_boundary() const{ return true;}  
  
 protected:
  //  DFT_Vec RHS0;
  //  DFT_Vec RHS1;

  DFT_Vec RHS_Boundary_;

  double *RHS0_sin_transform_ = NULL;
  double *RHS1_sin_transform_ = NULL;  

  fftw_plan sin_plan_;
  unsigned  sin_Ntot_;
  unsigned  sin_Norm_;
  double   *sin_in_;
  double   *sin_out_;
};
		    
*/
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
