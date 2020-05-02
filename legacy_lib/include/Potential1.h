#ifndef __LUTSKO_POTENTIAL1__
#define __LUTSKO_POTENTIAL1__

#include "Integrator.h"

/**
  *  @brief Potential function - it knows how to split into hard sphere and attractive parts and to compute the hard sphere diameter. 
  *
  *  @details By default, the split is done using a WCA prescription but this can be changed with the functions
  *           setBH and WCA_limit. 
  *
  *           The default WCA option writes v(r) = v0(r) + w(r) with v0(r) = 0 for r > Rmin, the minimum of the potential.
  *           and w(r) = v(Rmin) for r < Rmin, the minimum of the potential. 
  *
  *           If WCA_limit> 0, then w(r) is zero for r < WCA_limit: i.e. the continuation inside the core is truncated. 
  *
  *           If the BH flag is true, the potential is split at the point it reaches zero 
  *           and the attractive part is only non-zero for positions greater than this point.
  */
class Potential1
{
 public:
 Potential1(double sigma, double eps, double rcut) : sigma_(sigma), eps_(eps),rcut_(rcut), shift_(0), rmin_(0), Vmin_(0)
  {}

  virtual double getHardCore() const = 0; ///< Hard core diameter of the original potential (e.g. zero for LJ). 
  
  double getRcut() const { return rcut_;} ///< Returns the cutoff
  double getR0()   const { return r0_;} ///< Returns the position at which potential is zero
  virtual double getRmin() const = 0; ///< min of potential
  
  void set_WCA_limit(double r) { r_att_min_ = r;} ///< Distance from origin that attractive continuate extends
  void setBH() { bhFlag_ = true; r_att_min_ = getR0();} ///< Use the BH split of the potential
  
  double V(double r)   const { return vr(r)-shift_;}    ///< The cut and shifted potential at point r
  double V2(double r2) const { return vr2(r2)-shift_;}  ///< The cut and shifted potential at point r calculated from r^2

  double Watt(double r) const {return Watt2(r*r);} ///< The attractive tail 
  
    /*
  {
    double ret = 0.0;    
    if(r < rcut_)  
      if(bhFlag_) ret =  (r < getR0() ? 0.0 : V(r)); 
      else if(r < r_att_min_) ret = 0.0;  
      else if (r < rmin_) ret = Vmin_;    
      else ret = V(r);                    
    
    return ret;
  }
    */
  double Watt2(double r2) const ///< Attractive part calcuated from r2
  {
    double ret = 0.0;
    if(r2 < rcut_*rcut_)  // zero outside of cutoff      
      if(bhFlag_) ret =  (r2 < getR0()*getR0() ? 0.0 : V2(r2)); // BH case 
      else if(r2 < r_att_min_*r_att_min_) ret = 0.0; // Our "generalized" WCA
      else if (r2 < rmin_*rmin_) ret = Vmin_; // WCA continuation inside rmin
      else ret = V2(r2); // Just the potential outsize of rmin
    return ret;
  }  

  double V0(double r) const ///< The repulsive part of the potential
  {
    if(bhFlag_) return (r < getR0() ? V(r) : 0.0);
    return  (r < rmin_ ? V(r)-Vmin_ : 0.0);
  }
  
  double getHSD(double kT) const ///< Calculates HSD by numeric integration. 
  {
    kT_ = kT;
    double hc = getHardCore();

    double rlimit = rmin_;
    if(bhFlag_) rlimit = getR0();
    
    Integrator<const Potential1> II(this, &Potential1::dBH_Kernal, 1e-4, 1e-6);
    return hc + II.integrateFinite(hc,rlimit);
  }

  double getVDW_Parameter(double kT) const ///< Calculates the vDW parameter by numeric integration
  {
    kT_ = kT;
    Integrator<const Potential1> II(this, &Potential1::vdw_Kernal, 1e-6, 1e-8);
    if(bhFlag_) return (2*M_PI/kT)*II.integrateFinite(getR0(),rcut_);
    return (2*M_PI/kT)*(II.integrateFinite(r_att_min_,rmin_)+II.integrateFinite(rmin_,rcut_));
  }

  virtual string getIdentifier() const = 0; ///< Unique name of this potential
  
 protected:
  virtual double vr(double r)   const = 0; ///< The underlying potential
  virtual double vr2(double r2) const = 0; ///< The underlying potential

  virtual double dBH_Kernal(double r) const {return (1.0-exp(-V0(r)/kT_));}  ///< Kernal for calcuating hsd
  virtual double vdw_Kernal(double r) const {return r*r*Watt(r);}  ///< Kernel for calculating vDW parameter
  
 protected:
  double sigma_     =  1.0; ///< Length scale of potential
  double eps_       =  1.0;  ///< Energy scale of potential
  double rcut_      = -1.0;  ///< Cutoff
  double rmin_      =  1.0;   ///< Minimum of the potential
  double shift_     =  0.0;   ///< Amount to shift potential
  double Vmin_      = -1.0;   ///< Value of potential at its minimum
  double r_att_min_ =  0.0;   ///< parameter: wca continuation limit
  double r0_        =  1.0;   ///< Position at which potential goes to zero
  bool bhFlag_      = false;  ///< Flag to control which split to use
 
  mutable double kT_;         ///< Used to pass the temperature during evaluaton of numerical integrals. 
 
};

/**
  *  @brief Truncated Lennard-Jones potential
  */
class LJ : public Potential1
{
 public:
 LJ(double sigma, double eps, double rcut) : Potential1(sigma, eps, rcut)
  {
    shift_ = (rcut_ < 0 ? 0.0 : vr(rcut_));
    rmin_  = getRmin(); 
    Vmin_  = V(rmin_);
    r0_    = pow(0.5*sqrt(1+shift_)+0.5,-1.0/6.0);          
  }

  virtual double getRmin() const { return pow(2.0,1.0/6.0)*sigma_;}
  virtual double getHardCore() const { return 0.0;}
  /*
  virtual double getVDW_Parameter(double kT) const
  {
    double y0 = sigma_/rmin_;
    double yc = sigma_/rcut_;
    
    return (16*M_PI/9)*sigma_*sigma_*sigma_*(2*pow(y0,9)-3*pow(y0,3)+3*pow(yc,3)-2*pow(yc,9))/kT;
  }
  */

  virtual string getIdentifier() const
  {
    stringstream ss;
    ss << "LJ_" << sigma_ << "_" << eps_ << "_" << rcut_ << r_att_min_ << "_" << bhFlag_;
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
  virtual double vr2(double r2) const
  {
    double y2 = sigma_*sigma_/r2;
    double y6 = y2*y2*y2;
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
    shift_ = (rcut_ <= 0.0 ? 0.0 : vr(rcut_));
    rmin_  = getRmin();
    Vmin_  = V(rmin_);
    r0_    = sqrt(1+pow(25*sqrt(1+shift_)+25,-1.0/3.0));
  }

  virtual double getRmin() const { return sigma_*sqrt(1+pow(2.0/alpha_,1.0/3.0));}
  virtual double getHardCore() const { return sigma_;}
  
  // This is $a = \frac{1}{2}\int Watt(r)/kT d{\bf r}$
  // I assume that $kT = k_{B}T/\eps$
  /*  virtual double getVDW_Parameter(double kT) const
  {
    double xm = rmin_/sigma_;
    double xc = rcut_/sigma_;
    
    double ym = 1.0/(xm*xm-1);
    double yc = 1.0/(xc*xc-1);

    double gc = ((7+32*alpha_)/512)*log((xc-1)/(xc+1))-(1.0/3840)*(xc*yc*yc*yc*yc*yc*(105+xc*xc*(790+xc*xc*(-896+xc*xc*(490-105*xc*xc))))) + alpha_*xc*(xc*xc+1)*yc*yc/8;
    double gm = ((7+32*alpha_)/512)*log((xm-1)/(xm+1))-(1.0/3840)*(xm*ym*ym*ym*ym*ym*(105+xm*xm*(790+xm*xm*(-896+xm*xm*(490-105*xm*xm))))) + alpha_*xm*(xm*xm+1)*ym*ym/8;
    
    return (2*M_PI/kT)*sigma_*sigma_*sigma_*((xm*xm*xm/3)*vr(xm)-(xc*xc*xc/3)*vr(xc)+(4*eps_/(alpha_*alpha_))*(gc-gm));
  }
  */
  virtual string getIdentifier() const
  {
    stringstream ss;
    ss << "TWF_" << sigma_ << "_" << eps_ << "_" << alpha_ << "_" << rcut_ << r_att_min_ << "_" << bhFlag_;
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

  virtual double vr2(double r2) const
  {
      if(r2 < sigma_*sigma_) return 1e50;

      double s2 = r2/(sigma_*sigma_);
      double y = 1.0/(s2-1);
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
    
    shift_ = 0.0; // cutoff is built in
    rmin_  = getRmin(); 
    Vmin_  = V(rmin_);
    r0_    = 1.0;
  }

  virtual double getRmin() const { return rcut_*pow((1.0+2.0*(rcut_/sigma_)*(rcut_/sigma_))/3.0,-0.5);}
  virtual double getHardCore() const { return 0.0;}
  /*  virtual double getVDW_Parameter(double kT) const
  {
    double yc = rcut_/sigma_;
    double ym = getRmin()/rcut_;
    
    return 0.5*(4*M_PI*eps_/3)*(1.0/kT)*sigma_*sigma_*sigma_*pow(yc/(yc*yc-1),3)*(54*(1+yc*yc)-2*(ym*ym*ym)*(1+2*yc*yc)*(1+2*yc*yc)*(5+yc*yc));
  }
  */
  virtual string getIdentifier() const
  {
    stringstream ss;
    ss << "WHDF_" << sigma_ << "_" << eps_ << "_" << rcut_ << r_att_min_ << "_" << bhFlag_;
    return ss.str();    
  }  

 protected:
  virtual double vr(double r) const
  {
    double y = sigma_/r;
    double z = rcut_/r;
    return eps_*(y*y-1)*(z*z-1)*(z*z-1);
  }
  virtual double vr2(double r2) const
  {
    double y2 = sigma_*sigma_/r2;
    double z2 = rcut_*rcut_/r2;
    return eps_*(y2-1)*(z2-1)*(z2-1);
  }  
};




#endif // __POTENTIAL
