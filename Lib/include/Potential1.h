#ifndef __LUTSKO_POTENTIAL1__
#define __LUTSKO_POTENTIAL1__

#include "Integrator.h"

/**
  *  @brief Potential function - it knows how to split into hard sphere and attractive parts and to compute the hard sphere diameter.
  */
class Potential1
{
 public:
 Potential1(double sigma, double eps, double rcut) : sigma_(sigma), eps_(eps),rcut_(rcut), shift_(0), rmin_(0), Vmin_(0)
  { }

  virtual double getRmin() const = 0;
  virtual double getHardCore() const = 0;
  double getRcut() const { return rcut_;}
  
  double V(double r) const { return vr(r)-shift_;}

  double Watt(double r) const { return (r < rcut_ ? (r < rmin_ ? Vmin_ : V(r)) : 0);}
  double V0(double r)   const { return  (r < rmin_ ? V(r)-Vmin_ : 0.0);}
  
  double getHSD(double kT) const
  {
    kT_ = kT;
    double hc = getHardCore();
    
    Integrator<const Potential1> II(this, &Potential1::dBH_Kernal, 1e-4, 1e-6);
    return hc + II.integrateFinite(hc,rmin_);
  }

  // This is $a = \frac{1}{2}\int Watt(r)/kT d{\bf r}$
  // I assume that $kT = k_{B}T/\eps$
  virtual double getVDW_Parameter(double kT) const = 0;

  virtual string getIdentifier() const = 0;

  
 protected:
  virtual double vr(double r) const = 0;

  virtual double dBH_Kernal(double r) const {return (1.0-exp(-V0(r)/kT_));}
  
 protected:
 double sigma_;
 double eps_;
 double rcut_;
 double rmin_;
 double shift_;
 mutable double kT_;
 double Vmin_;
};

/**
  *  @brief Truncated Lennard-Jones potential
  */
class LJ : public Potential1
{
 public:
 LJ(double sigma, double eps, double rcut) : Potential1(sigma, eps, rcut)
  {
    shift_ = vr(rcut_);
    rmin_  = getRmin(); //pow(2.0,1.0/6.0)*sigma;
    Vmin_  = V(rmin_);
  }

  virtual double getRmin() const { return pow(2.0,1.0/6.0)*sigma_;}
  virtual double getHardCore() const { return 0.0;}
  virtual double getVDW_Parameter(double kT) const
  {
    double y0 = sigma_/rmin_;
    double yc = sigma_/rcut_;
    
    return (16*M_PI/9)*sigma_*sigma_*sigma_*(2*pow(y0,9)-3*pow(y0,3)+3*pow(yc,3)-2*pow(yc,9))/kT;
  }


  virtual string getIdentifier() const
  {
    stringstream ss;
    ss << "LJ_" << sigma_ << "_" << eps_ << "_" << rcut_;
    return ss.str();    
  }
  
 protected:
  virtual double vr(double r) const
  {
    double y = sigma_/r;
    double y3 = y*y*y;
    double y6 = y3*y3;
    return 4*eps_*(y6*y6-y6);
  }
};

/**
  *  @brief Truncated ten Wolde-Frenkel potential
  */

class tWF : public Potential1
{
 public:
 tWF(double sigma, double eps, double rcut, double alpha = 50.0) : Potential1(sigma, eps, rcut), alpha_(alpha)
  {
    shift_ = vr(rcut_);
    rmin_  = getRmin();
    Vmin_  = V(rmin_);

  }

  virtual double getRmin() const { return sigma_*sqrt(1+pow(2.0/alpha_,1.0/3.0));}
  virtual double getHardCore() const { return sigma_;}
  
  // This is $a = \frac{1}{2}\int Watt(r)/kT d{\bf r}$
  // I assume that $kT = k_{B}T/\eps$
  virtual double getVDW_Parameter(double kT) const
  {
    double xm = rmin_/sigma_;
    double xc = rcut_/sigma_;
    
    double ym = 1.0/(xm*xm-1);
    double yc = 1.0/(xc*xc-1);

    double gc = ((7+32*alpha_)/512)*log((xc-1)/(xc+1))-(1.0/3840)*(xc*yc*yc*yc*yc*yc*(105+xc*xc*(790+xc*xc*(-896+xc*xc*(490-105*xc*xc))))) + alpha_*xc*(xc*xc+1)*yc*yc/8;
    double gm = ((7+32*alpha_)/512)*log((xm-1)/(xm+1))-(1.0/3840)*(xm*ym*ym*ym*ym*ym*(105+xm*xm*(790+xm*xm*(-896+xm*xm*(490-105*xm*xm))))) + alpha_*xm*(xm*xm+1)*ym*ym/8;
    
    return (2*M_PI/kT)*sigma_*sigma_*sigma_*((xm*xm*xm/3)*vr(xm)-(xc*xc*xc/3)*vr(xc)+(4*eps_/(alpha_*alpha_))*(gc-gm));
  }

  virtual string getIdentifier() const
  {
    stringstream ss;
    ss << "TWF_" << sigma_ << "_" << eps_ << "_" << alpha_ << "_" << rcut_;
    return ss.str();    
  }
  
 protected:
  virtual double vr(double r) const
  {
      if(r < sigma_) return 1e50;

      double s = r/sigma_;
      double y = 1.0/(s*s-1);
      double y3 = y*y*y;
      
      return (4*eps_/(alpha_*alpha_))*(y3*y3-alpha_*y3);
  }

  double alpha_;
};

/**
  *  @brief Frenkel's potential from https://arxiv.org/pdf/1910.05746.pdf
  */
class WHDF : public Potential1
{
 public:
 WHDF(double sigma, double eps, double rcut) : Potential1(sigma, eps, rcut)
  {
    eps_ *= 2*pow(rcut/sigma,2)*pow(2*((rcut/sigma)*(rcut/sigma)-1)/3.0,-3.0);
    
    shift_ = vr(rcut_);
    rmin_  = getRmin(); 
    Vmin_  = V(rmin_);
    
  }

  virtual double getRmin() const { return rcut_*pow((1.0+2.0*(rcut_/sigma_)*(rcut_/sigma_))/3.0,-0.5);}
  virtual double getHardCore() const { return 0.0;}
  virtual double getVDW_Parameter(double kT) const
  {
    double yc = rcut_/sigma_;
    double ym = getRmin()/rcut_;
    
    return (4*M_PI*eps_/3)*sigma_*sigma_*sigma_*pow(yc/(yc*yc-1),3)*(54*(1+yc*yc)-2*(ym*ym*ym)*(1+2*yc*yc)*(1+2*yc*yc)*(5+yc*yc));
  }

  virtual string getIdentifier() const
  {
    stringstream ss;
    ss << "WHDF_" << sigma_ << "_" << eps_ << "_" << rcut_;
    return ss.str();    
  }  

 protected:
  virtual double vr(double r) const
  {
    double y = sigma_/r;
    double z = rcut_/r;
    return eps_*(y*y-1)*(z*z-1);
  }
};




#endif // __POTENTIAL
