#ifndef __LUTSKO_CONJUG2__
#define __LUTSKO_CONJUG2__

using namespace std;

#include "DFT.h"

/**
  *  @brief Minimizer base class
  *
  */  

class Minimizer
{
 public:
 Minimizer(DFT &dft, ostream& log) : dft_(dft), forceLimit_(0.1),  bFrozenBoundary_(false), log_(log)
  {
    x_.resize(dft.getNumberOfSpecies());

    for(auto &x: x_) x.resize(dft_.lattice().Ntot()); 
  }

  void setFrozenBoundaryFlag(bool f) {bFrozenBoundary_ = f;}
  
  virtual void initialize();

  void run(long maxSteps = -1);

  virtual double step() = 0;

  virtual int draw_before() // Display something before the next step
  {
    stringstream ts;
    ts << "calls = " << calls_ << " F = " << F_ << " err = " << f_abs_max_;  //" N = " << density_.getNumberAtoms();

    string title = ts.str();
    string file("image_current.png");
    dft_.doDisplay(title, file);
    return 1;
  }
  virtual int draw_during() = 0;  // Display something during the minimization
  virtual void draw_after()  = 0;  // Display something after the minimization

  int getCalls() const { return calls_;}
  double getF() const { return F_;}

  double getForceTerminationCriterion() const {return forceLimit_;}
  void   setForceTerminationCriterion(double v) {forceLimit_ = v;}
  
  virtual double getDF_DX();
  virtual double get_convergence_monitor() const { return dft_.get_convergence_monitor();}
  
 protected:
  DFT &dft_;

  vector<DFT_Vec> x_;

  int calls_ = 0;
  int step_counter_ = 0;
  double F_ = 0;

  double forceLimit_;
  double f_abs_max_; // max absolute value of dF_
  bool bFrozenBoundary_;
  ostream &log_;

  double vv_ = 0;

};

/**
  *  @brief Base class for a family of  finite elements ddft integrators. 
  *
  */  

class DDFT : public Minimizer
{
 public:
 DDFT(DFT &dft, ostream &log,   bool showGraphics = true)
   : Minimizer(dft, log), show_(showGraphics) , tolerence_fixed_point_(1e-4), successes_(0), fixedBorder_(false), modified_(false)
  {
    double dx = dft_.lattice().getDX();
    dt_ = 10*0.1*dx*dx;
    dt_ = 0.0001*dx*dx;
    dtMax_ = 1*dx*dx;
  }
  ~DDFT() {}

  virtual void initialize();
  
  virtual int draw_during(){ return 1;}  // Display something during the minimization
  virtual void draw_after()  // Display something after the minimization
  {
    cout << "After DDFT step " << step_counter_ 
	 << " F = " << F_ 
	 << " max(dF) = " << f_abs_max_
	 << " and N = " << dft_.getDensity(0).getNumberAtoms() 
	 << endl;
  }

  void Display(double F, double dFmin, double dFmax, double N);
    
  void set_tolerence_fixed_point(double e) { tolerence_fixed_point_ = e;}
  void set_max_time_step(double t) { dtMax_ = t;}
 void setTimeStep(double dt) { dt_ = dt;}

  void setFixedBoundary() {fixedBorder_ = true;}
  void setModified() {modified_ = true;}
    
  double getTimeStep() const { return dt_;}
  
  double F_string(Density &d, double *fmax = NULL);
  void reverseForce(DFT_Vec *tangent);

  virtual double step() = 0;
  virtual double step_string(double &dt, Density &d, unsigned &time_den, bool verbose = true) = 0;

 protected:

  bool show_;
  
  double dt_;
  double tolerence_fixed_point_;

  // control of adaptive time step
  int successes_;
  double dtMax_;

  bool fixedBorder_;
  bool modified_; // for certain derived classes
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
 DDFT_IF(DFT &dft, ostream &log,  bool showGraphics = true)
   : DDFT(dft,  log, showGraphics) {}

  ~DDFT_IF() {}

  virtual void initialize();

  virtual double step();

  virtual double step_string(double &dt, Density &d, unsigned &time_den, bool verbose = true);

  double fftDiffusion(const Density &density, DFT_Vec &d1, const DFT_FFT &RHS0, const DFT_FFT &RHS1);
  void calcNonlinearTerm(const DFT_Vec &d2, const DFT_Vec &dF, DFT_Vec &RHS1);
  void restore_values_on_border(DFT_Vec& d1, const DFT_Vec &d0);

  virtual double get_convergence_monitor() const { return dft_.get_convergence_monitor();}

 protected:

  double largest_change_on_border_; // for reporting effect of restore_values_on_border()
};

/**
  *  @brief DDFT minimizer Class using integrating factor. This is for fixed density at the boundaries - a model for an open system.
  *
  *  @detailed This integrates the pure diffusion part of the dynamics exactly (using FFTs) and treats the rest implicitly via a Crank-Nicholson type method.
  */  

class DDFT_IF_Open : public DDFT
{
 public:
 DDFT_IF_Open(DFT &dft, double background, ostream &log,  bool showGraphics = true)
   : DDFT(dft, log, showGraphics), background_(background), sin_in_(NULL), sin_out_(NULL)
    {}
  ~DDFT_IF_Open() {if(sin_in_) delete sin_in_; if(sin_out_) delete sin_out_;}

  virtual void initialize();

  virtual double step();

  virtual double step_string(double &dt, Density &d, unsigned &time_den, bool verbose = true);

  double fftDiffusion(const Density &density, DFT_Vec &d1, const DFT_FFT &RHS0, const DFT_FFT &RHS1) {throw std::runtime_error("Need to adapt fftDiffusion for non-string application");}
  double fftDiffusion(DFT_Vec &d1, const double *RHS0_sin_transform, const double *RHS1_sin_transform);
  void   calcNonlinearTerm(const DFT_Vec &d2, const DFT_Vec &dF, DFT_Vec &RHS1);

  void pack_for_sin_transform(const double *x, double val);
  
  void unpack_after_transform(double *x, double val);

  
 protected:
  double background_;
  
  fftw_plan sin_plan_;
  unsigned sin_Ntot_;
  unsigned sin_Norm_;
  double *sin_in_;
  double *sin_out_;

  vector<double> Lamx;
  vector<double> Lamy;
  vector<double> Lamz;
};

/**
  *  @brief Minimizer using FIRE algorithm
  *
  */

class fireMinimizer_Mu : public Minimizer
{
 public:
 fireMinimizer_Mu(DFT &dft, ostream &log) :  Minimizer(dft, log)
  {
    v_.resize(dft_.getNumberOfSpecies());
    for(auto &v: v_) v.resize(dft_.lattice().Ntot());

    dt_ = 1e-3; // my guess
    dt_max_ = 10*dt_;
    
    N_delay_ = 5;
  
    alpha_start_ = 0.1;
    f_dec_ = 0.5;
    f_inc_ = 1.1;
    f_alf_ = 0.99;
  }

  void onlyRelaxSpecies(int j) { onlyRelax_ = j;}
  
  virtual void initialize();

  virtual double step();

  virtual int  draw_during();
  virtual void draw_after();

  void verlet();

  void setTimeStep(double dt)    { dt_ = dt;}
  void setTimeStepMax(double dt) { dt_max_ = dt;}  
  void setAlphaStart(double a)   { alpha_start_ = a;}
  void setAlphaFac(double a)     { f_alf_ = a;}

 protected:
  vector<DFT_Vec> v_;

  int onlyRelax_ = -1;

  double alpha_start_ = 0.1;
  double f_dec_    = 0.5;
  double f_inc_    = 1.1;
  double f_alf_    = 0.99;
  double dt_max_;
  double dt_;

  unsigned it_;     ///< Loop counter
  unsigned cut_;    ///< used to keep track of number of steps P >0
  int N_delay_ = 5;
  double alpha_;

};


/**
  *  @brief Minimizer using FIRE2 algorithm
  *
  */

class fireMinimizer2 : public Minimizer //fireMinimizer_Mu
{
 public:
  fireMinimizer2(DFT &dft, ostream &log);

  virtual void   initialize();
  virtual double step();

  double getRMS_Force() const { return rms_force_;}

  void onlyRelaxSpecies(int j) { onlyRelax_ = j;}  
  
  virtual int  draw_during(){ log_ << "\t1D minimization F-mu*N = " << F_  << " N = " << dft_.getNumberAtoms(0) << endl;  return 1;}
  //  virtual void draw_after() { cout << "dt_ = " << dt_ << endl;}

  void setTimeStep(double dt)    { dt_ = dt;}
  void setTimeStepMax(double dt) { dt_max_ = dt;}  
  void setAlphaStart(double a)   { alpha_start_ = a;}
  void setAlphaFac(double a)     { f_alf_ = a;}
  void setBacktrackFac(double a) { f_back_ = a;}  
  
 protected:
  void SemiImplicitEuler(int begin, int end);
  
 protected:
 vector<DFT_Vec> v_;

  int onlyRelax_ = -1;

  double alpha_start_ = 0.1;
  double f_dec_    = 0.5;
  double f_inc_    = 1.1;
  double f_alf_    = 0.99;
  double f_back_   = 0.8;
  double dt_max_;
  double dt_;

  unsigned it_;     ///< Loop counter
  unsigned cut_;    ///< used to keep track of number of steps P >0
  int N_delay_ = 5;
  double alpha_;
  
  int N_P_positive_ = 0;
  int N_P_negative_ = 0;

  double rms_force_ = 0.0; // for reporting
  double vnorm_     = 0.0; // for reporting
  
  double dt_min_ = 0.0;        ///< minimum allowed timestep
  bool initial_delay_ = true; ///< flag to allow for a warmup period
  int N_P_negative_max_ = 20;
  int steps_since_backtrack_ = 0;

  double threshold_ = 0.0;
  double fmax_ = 0.0;
};

class adamMinimizer : public Minimizer
{
 public:
 adamMinimizer(DFT &dft, ostream &log) :  Minimizer(dft, log)
    {
      m_.resize(dft_.getNumberOfSpecies());
      s_.resize(dft_.getNumberOfSpecies());      
      for(auto &m: m_) m.resize(dft_.lattice().Ntot());
      for(auto &s: s_) s.resize(dft_.lattice().Ntot());
    }

  virtual void   initialize()
  {
    Minimizer::initialize();
    
    F_ = getDF_DX();

    beta1_t_ = beta1_;
    beta2_t_ = beta2_;
    
    for(auto &m: m_)
      m.zeros(m.size());
    for(auto &s: s_)
      s.zeros(s.size());    
  }

  void setTimeStep(double dt) { dt_ = dt_start_ = dt;}

  virtual double step()
  {
    int numSpecies = dft_.getNumberOfSpecies();
    
    vector<DFT_Vec> x_rem(numSpecies);
    vector<DFT_Vec> m_rem(numSpecies);
    vector<DFT_Vec> s_rem(numSpecies);

    for(int Jspecies =0; Jspecies < numSpecies; Jspecies++)
      {
	x_rem[Jspecies].set(x_[Jspecies]);
	m_rem[Jspecies].set(m_[Jspecies]);
	s_rem[Jspecies].set(s_[Jspecies]);
      }

    long Ntot = x_[0].size();

    for(int Jspecies = 0; Jspecies < numSpecies; Jspecies++)
      for(long i = 0; i<Ntot; i++)
	{
	  double d = dft_.getDF(Jspecies).get(i);
	  double mm = m_[Jspecies].get(i)*beta1_+(1-beta1_)*d;
	  double ss = s_[Jspecies].get(i)*beta2_+(1-beta2_)*d*d;
	  
	  m_[Jspecies].set(i, mm);
	  s_[Jspecies].set(i, ss);
	}

    try {
      for(int Jspecies = 0; Jspecies < numSpecies; Jspecies++)
	for(long i = 0; i<Ntot; i++)
	  {
	    double mm = m_[Jspecies].get(i);
	    double ss = s_[Jspecies].get(i);	    
	    x_[Jspecies].set(i,x_[Jspecies].get(i) - dt_*(mm/(1-beta1_t_))/(small_eps_ + sqrt(ss/(1-beta2_t_))));
	  }
      F_ = getDF_DX();
      dt_ = min(dt_*2, dt_start_);
    } catch (Eta_Too_Large_Exception &e) {
      dt_ /= 2;
      log_ << "Backtrack .. dt_ is now " << dt_ << endl;
    }
    beta1_t_ *= beta1_;
    beta2_t_ *= beta2_;
    return F_;
  }
  virtual void draw_after() { cout << "dt_ = " << dt_ << endl;}
  double getRMS_Force() const { return rms_force_;}

int draw_during()
{
  log_ << "\t1D minimization F-mu*N = " << F_  << " N = " << dft_.getNumberAtoms(0) << endl;
  return 1;
}
  
 protected:
  double beta1_ = 0.9;
  double beta2_ = 0.99;

  double beta1_t_ = 0.9;
  double beta2_t_ = 0.99;
  
  double small_eps_ = 1e-8;
  vector<DFT_Vec> m_;
  vector<DFT_Vec> s_;

  double dt_ = 1;
  double dt_start_ = 1;
  
  double rms_force_ = 0.0; // for reporting
};


class picardMinimizer : public Minimizer
{
 public:
 picardMinimizer(DFT &dft, ostream &log, double alpha) :  Minimizer(dft, log), alpha_(alpha) {}

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
	log_ << "Backtrack .. alpha_ is now " << alpha_ << endl;
	
	for(int Jspecies = 0; Jspecies < numSpecies; Jspecies++)
	  dft_.set_density(Jspecies, x_rem[Jspecies]);	
      }
    } while(!bSuccess);
    return F_;
  }


  
  virtual void draw_after() { cout << "alpha_ = " << alpha_ << endl;}
  double getRMS_Force() const { return rms_force_;}

  int draw_during()
  {
    log_ << "\t1D minimization F-mu*N = " << F_  << " N = " << dft_.getNumberAtoms(0) << endl;
    return 1;
  }
  
 protected:
  double alpha_ = 0.1;
  double rms_force_ = 0.0; // for reporting
};



#endif // sentinal
