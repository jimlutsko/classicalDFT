#ifndef __LUTSKO_POTENTIAL__
#define __LUTSKO_POTENTIAL__

#include "Integrator.h"

class Potential1
{
 public:
 Potential1(double sigma, double eps, double rcut) : sigma_(sigma), eps_(eps),rcut_(rcut)
  {
    shift_ = vr(rcut_);
    rmin_  = pow(2.0,1.0/6.0)*sigma;
    Vmin_  = V(rmin_);
  }

  double V(double r) const { return vr(r)-shift_;}

  double Watt(double r) const { return (r < rcut_ ? (r < rmin_ ? Vmin_ : V(r)) : 0);}
  double V0(double r)   const { return  (r < rmin_ ? V(r)-Vmin_ : 0.0);}
  
  double getHSD(double kT) const
  {
    kT_ = kT; 
    Integrator<const Potential1> II(this, &Potential1::dBH_Kernal, 1e-4, 1e-6);
    return II.integrateFinite(0.0,rmin_);
  }


  // This is $a = \frac{1}{2}\int Watt(r)/kT d{\bf r}$
  // I assume that $kT = k_{B}T/\eps$
  double getVDW_Parameter(double kT) const
  {
    double y0 = sigma_/rmin_;
    double yc = sigma_/rcut_;
    
    return (16*M_PI/9)*sigma_*sigma_*sigma_*(2*pow(y0,9)-3*pow(y0,3)+3*pow(yc,3)-2*pow(yc,9))/kT;
  }
  
 private:
  double vr(double r) const
  {
    double y = sigma_/r;
    double y3 = y*y*y;
    double y6 = y3*y3;
    return 4*eps_*(y6*y6-y6);
  }

  virtual double dBH_Kernal(double r) const {return (1.0-exp(-V0(r)/kT_));}
  
 private:
 double sigma_;
 double eps_;
 double rcut_;
 double rmin_;
 double shift_;
 mutable double kT_;
double Vmin_;
};

#endif // __POTENTIAL
