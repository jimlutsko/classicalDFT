#ifndef __LUTSKO_PICARD__
#define __LUTSKO_PICARD__

#include <mgl2/mgl.h>
#include <mgl2/fltk.h>

#include "DFT.h"

/**
  *  @brief Conjugate Gradients Class
  *
  */  

class Picard
{
 public:
 Picard(DFT &dft, Density &density, double mu) : dft_(dft), density_(density), mu_(mu), Mixing_(0.0)
  {
    x_.resize(density_.Ntot());
    for(long i=0;i<density_.Ntot();i++)
      x_(i) = sqrt(max(0.0, density_.getDensity(i) - SMALL_VALUE));
    dF_.resize(x_.size());
  }
  
  void run(string& logfile);
  void display();
  double step();


  void check(bool write);

  virtual void initialize();
  virtual void finish(const char *){}
  virtual double getDF_DX();

  virtual int draw_before() {} // Display what ever you want before the next step
  virtual int draw_during();  // Display what ever you want during the minimization
  virtual int draw_after() {    cout << "After picard step " << step_counter_ << " F-mu*N = " << F_ << " and N = " << density_.getNtotal() << endl;} // Display what ever you want  after the next step

  const Density & getDensity() const { return density_;}
  int getCalls() const { return calls_;}
  double getF() const { return F_;}

  void set_mixing_parameter(double Mixing) {Mixing_ = Mixing;}

  bool continuePicard() const { return 1;}

 protected:
  DFT &dft_;
  Density &density_;
  double mu_;

  // Working space for the minimization

  DFT_Vec dF_;
  DFT_Vec x_;

  int calls_ = 0;
  int step_counter_ = 0;
  double F_ = 0;

  double Mixing_;
};


class Picard_1D_Graphics : public Picard
{
 public:
 Picard_1D_Graphics(DFT &dft, Density &density, double mu, void (*display)(Picard &cg)) : Picard(dft, density, mu), display_(display)
  {}

  virtual int draw_before(){ (*display_)(*this);}
  virtual void finish(const char *) { }

 protected:
  void (*display_)(Picard &cg);
};

class Picard_2D_Graphics : public Picard
{
 public:
 Picard_2D_Graphics(DFT &dft, Density &density, double mu, bool show = true) : Picard(dft, density, mu), show_(show)
  {}

  virtual void initialize();

  virtual int draw_before();
  virtual void finish(const char *c) { gr_.WriteFrame(c);}

 protected:
  mglFLTK gr_;
  bool show_;
};


// Here we minimize  F[x*x*(N/sum(x*x*dV))] over x.
// There are a couple of problems.
// First,  the calculation of dF/dX is not completely correct - there are some constant factors missing
// Second, since the total mass is fixed, there are really only Ntot-1 independent variables. We should 
// really reduce the size of the variable array by one to account for this ... 

class Picard_Fixed : public Picard_2D_Graphics
{
 public:
 Picard_Fixed(DFT &dft, Density &density, double Ntarget, bool showGraphics = true) : Picard_2D_Graphics(dft, density, 0.0, showGraphics), N_fixed_target_(Ntarget)
  {}

  virtual int draw_during();
  virtual int draw_after();

  virtual double getDF_DX();

 protected:
  double N_fixed_target_;
  double mu_eff_;
};


#endif // sentinal
