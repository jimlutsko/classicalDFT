#ifndef __LUTSKO_EOS__
#define __LUTSKO_EOS__

#include <iostream>
#include <vector>
#include <cmath>


#ifdef USE_GRACE
#include "Grace.h"
#endif

/**
  *  @brief Equation of state interface
  */  

class EOS
{
 public:
  /**
   *   @brief  Default  constructor for a FCC cell in a periodic box 
   *  
   *   @param  dx is lattice spacing: assumed to be the same in all directions
   *   @param  L[] are the dimensions of the physical region (FCC cube side)
   *   @return nothing 
   */  
  EOS() {}

   ~EOS(){}

   virtual char const * const getName() const =0;

   // F/N
   double freeEnergy(double density, double kT) const { return kT*(log(density)-1)+phix(density,kT);}

   // P : the pressure is n^2 dphi/dn = n^2d(f/n)/dn = n(df/dn)-f = n(dn*phi/dn)- n*phi= n^2(dphi/dn)
   double Pressure(double density, double kT) const { return kT*density+density*density*dphix_dx(density,kT);}

   double dP_dx(double x, double kT)    const { return kT+2*x*dphix_dx(x,kT) + x*x*dphix_dx2(x,kT);}
   double dP_dxdT(double x, double kT)  const { return 1.0+2*x*dphix_dxdT(x,kT) + x*x*dphix_dx2dT(x,kT);}
   double dP_dx2(double x, double kT)   const { return 2*dphix_dx(x,kT)+4*x*dphix_dx2(x,kT)+x*x*dphix_dx3(x,kT);}
   double dP_dx2dT(double x, double kT) const { return 2*dphix_dxdT(x,kT)+4*x*dphix_dx2dT(x,kT)+x*x*dphix_dx3dT(x,kT);}
   
   // excess free energy per atom:  F/N
   virtual double phix(double density, double kT)  const = 0;

   // derivative of excess free energy per atom :  d(F/N)/dx    
   virtual double dphix_dx(double density, double kT) const = 0;

   // 2nd derivative of excess free energy per atom :  d2(F/N)/dx2       
   virtual double dphix_dx2(double density, double kT) const = 0;

   // 3rd derivative of excess free energy per atom :  d3(F/N)/d32       
   virtual double dphix_dx3(double density, double kT) const = 0;

   virtual double dphix_dxdT(double density, double kT) const = 0;
   virtual double dphix_dx2dT(double density, double kT) const = 0;
   virtual double dphix_dx3dT(double density, double kT) const = 0;

   bool findCriticalPoint(double &density, double &kT) const;
   bool findCoexistence(double &x1, double &x2, double &kT) const;
   bool findSpinodal(double &x1, double &x2, double &kT) const;
   
};


class JZG : public EOS
{
 public:
 JZG(double rc = -1.0, bool doShift = true)
   : EOS(), gamma(3.0), da_(0.0)
    {
      x1 = 0.8623085097507421;
      x2 = 2.976218765822098;
      x3 = -8.402230115796038;
      x4 = 0.1054136629203555;
      x5 = -0.8564583828174598;
      x6 = 1.582759470107601;
      x7 = 0.7639421948305453;
      x8 = 1.753173414312048;
      x9 = 2.798291772190376e3;
      x10 = -4.8394220260857657e-2;
      x11 = 0.9963265197721935;
      x12 = -3.698000291272493e1;
      x13 = 2.084012299434647e1;
      x14 = 8.305402124717285e1;
      x15 = -9.574799715203068e2;
      x16 = -1.477746229234994e2;
      x17 = 6.398607852471505e1;
      x18 = 1.603993673294834e1;
      x19 = 6.805916615864377e1;
      x20 = -2.791293578795945e3;
      x21 = -6.245128304568454;
      x22 = -8.116836104958410e3;
      x23 = 1.488735559561229e1;
      x24 = -1.059346754655084e4;
      x25 = -1.131607632802822e2;
      x26 = -8.867771540418822e3;
      x27 = -3.986982844450543e1;
      x28 = -4.689270299917261e3;
      x29 = 2.593535277438717e2;
      x30 = -2.694523589434903e3;
      x31 = -7.218487631550215e2;
      x32 = 1.721802063863269e2;

      if(rc > 0.0)
	{
	  da_ = -(32.0/9.0)*M_PI*(pow(rc,-9.0)-1.5*pow(rc,-3.0));
	  if(doShift) da_ += (8.0/3.0)*M_PI*(pow(rc,-9.0)-pow(rc,-3.0));
	}
    }

  virtual ~JZG() {}

  virtual char const * const getName() const { return "JZG";}

  // excess free energy per atom and density derivatives
  virtual double phix(double density, double kT)  const 
  { 
      double f  = 0.0;
      for(int i=1;i<=8;i++) f += a(kT,i)*pow(density,i)/i;
      for(int i=1;i<=6;i++) f += b(kT,i)*G(density,i);
      f += da_*density;
      return f;
  }

  virtual double dphix_dx(double density, double kT) const 
  { 
      double F = exp(-gamma*density*density);

      double f  = 0.0;
      for(int i=1;i<=8;i++) f += a(kT,i)*pow(density,i-1);
      for(int i=1;i<=6;i++) f += b(kT,i)*F*pow(density,2*i-1);
      f += da_;
      return f;
  }

  virtual double dphix_dx2(double density, double kT) const
  { 
      double F = exp(-gamma*density*density);

      double f  = 0.0;
      for(int i=2;i<=8;i++) f += a(kT,i)*(i-1)*pow(density,i-2);
      for(int i=1;i<=6;i++) f += b(kT,i)*F*((2*i-1)*pow(density,2*i-2)-2*gamma*pow(density,2*i));
      return f;
  }
  virtual double dphix_dx3(double density, double kT) const 
  { 
      double F = exp(-gamma*density*density);

      double f  = 0.0;
      for(int i=3;i<=8;i++) f += a(kT,i)*(i-1)*(i-2)*pow(density,i-3);
      for(int i=1;i<=6;i++) f += b(kT,i)*(-2*gamma*density)*F*((2*i-1)*pow(density,2*i-2)-2*gamma*pow(density,2*i));
      for(int i=2;i<=6;i++) f += b(kT,i)*F*(2*i-1)*(2*i-2)*pow(density,2*i-3);
      for(int i=1;i<=6;i++) f += b(kT,i)*F*(-2*2*i*gamma*pow(density,2*i-1));
      return f;
  }
// excess free energy per atom and density derivatives
  virtual double dphix_dT(double density, double kT)  const 
  { 
      double f  = 0.0;
      for(int i=1;i<=8;i++) f += da_dT(kT,i)*pow(density,i)/i;
      for(int i=1;i<=6;i++) f += db_dT(kT,i)*G(density,i);
      f += da_*density;
      return f;
  }

  virtual double dphix_dxdT(double density, double kT) const 
  { 
      double F = exp(-gamma*density*density);

      double f  = 0.0;
      for(int i=1;i<=8;i++) f += da_dT(kT,i)*pow(density,i-1);
      for(int i=1;i<=6;i++) f += db_dT(kT,i)*F*pow(density,2*i-1);
      f += da_;
      return f;
  }

  virtual double dphix_dx2dT(double density, double kT) const
  { 
      double F = exp(-gamma*density*density);

      double f  = 0.0;
      for(int i=2;i<=8;i++) f += da_dT(kT,i)*(i-1)*pow(density,i-2);
      for(int i=1;i<=6;i++) f += db_dT(kT,i)*F*((2*i-1)*pow(density,2*i-2)-2*gamma*pow(density,2*i));
      return f;
  }
  virtual double dphix_dx3dT(double density, double kT) const 
  { 
      double F = exp(-gamma*density*density);

      double f  = 0.0;
      for(int i=3;i<=8;i++) f += da_dT(kT,i)*(i-1)*(i-2)*pow(density,i-3);
      for(int i=1;i<=6;i++) f += db_dT(kT,i)*(-2*gamma*density)*F*((2*i-1)*pow(density,2*i-2)-2*gamma*pow(density,2*i));
      for(int i=2;i<=6;i++) f += db_dT(kT,i)*F*(2*i-1)*(2*i-2)*pow(density,2*i-3);
      for(int i=1;i<=6;i++) f += db_dT(kT,i)*F*(-2*2*i*gamma*pow(density,2*i-1));
      return f;
  }  




 private:

  double a(double T, int i) const
  {
    double T1 = 1.0/T;
    double T2 = T1*T1;
    double T3 = T1*T2;
    double Ts = sqrt(T);
    
    if(i == 1) return x1*T+x2*Ts+x3+x4*T1+x5*T2;
    if(i == 2) return x6*T+x7+x8*T1+x9*T2;
    if(i == 3) return x10*T+x11+x12*T1;
    if(i == 4) return x13;
    if(i == 5) return x14*T1+x15*T2;
    if(i == 6) return x16*T1;
    if(i == 7) return x17*T1+x18*T2;
    if(i == 8) return x19*T2;
    throw std::runtime_error("Unknown values in LJ");
  }

  double da_dT(double T, int i) const
  {
    double T1 = 1.0/T;
    double T2 = T1*T1;
    double T3 = T1*T2;
    double T4 = T1*T3;
    double Ts = sqrt(T);

    double dT1 = -T2;
    double dT2 = -2*T3;
    double dT3 = -3*T4;
    double dTs = 0.5/Ts;
    
    if(i == 1) return x1+x2*dTs+x4*dT1+x5*dT2;
    if(i == 2) return x6+x8*dT1+x9*dT2;
    if(i == 3) return x10+x12*dT1;
    if(i == 4) return 0;
    if(i == 5) return x14*dT1+x15*dT2;
    if(i == 6) return x16*dT1;
    if(i == 7) return x17*dT1+x18*dT2;
    if(i == 8) return x19*dT2;
    throw std::runtime_error("Unknown values in LJ");
  }

  double b(double T, int i) const 
  {
    double T1 = 1.0/T;
    double T2 = T1*T1;
    double T3 = T1*T2;
    double T4 = T2*T2;

    if(i == 1) return x20*T2+x21*T3;
    if(i == 2) return x22*T2+x23*T4;
    if(i == 3) return x24*T2+x25*T3;
    if(i == 4) return x26*T2+x27*T4;
    if(i == 5) return x28*T2+x29*T3;
    if(i == 6) return x30*T2+x31*T3+x32*T4;
    throw std::runtime_error("Unknown values in LJ");
  }

  double db_dT(double T, int i) const 
  {
    double T1 = 1.0/T;
    double T2 = T1*T1;
    double T3 = T1*T2;
    double T4 = T2*T2;
    double T5 = T3*T2;
    
    double dT1 = -T2;
    double dT2 = -2*T3;
    double dT3 = -3*T4;
    double dT4 = -4*T5;

    if(i == 1) return x20*dT2+x21*dT3;
    if(i == 2) return x22*dT2+x23*dT4;
    if(i == 3) return x24*dT2+x25*dT3;
    if(i == 4) return x26*dT2+x27*dT4;
    if(i == 5) return x28*dT2+x29*dT3;
    if(i == 6) return x30*dT2+x31*dT3+x32*dT4;
    throw std::runtime_error("Unknown values in LJ");
  }  

  double G(double d, int i) const 
  {
    double F = exp(-gamma*d*d);

    double G1 = (1-F)/(2*gamma);
    if(i == 1) return G1;

    double G2 = -(F*pow(d,2)-2*G1)/(2*gamma);
    if(i == 2) return G2;

    double G3 = -(F*pow(d,4)-4*G2)/(2*gamma);
    if(i == 3) return G3;

    double G4 = -(F*pow(d,6)-6*G3)/(2*gamma);
    if(i == 4) return G4;

    double G5 = -(F*pow(d,8)-8*G4)/(2*gamma);
    if(i == 5) return G5;

    double G6 = -(F*pow(d,10)-10*G5)/(2*gamma);
    if(i == 6) return G6;

    throw std::runtime_error("Unknown values in LJ");
  }


  double gamma; // fixed parameter (could be adjusted???)
  double da_; // correction to fe due to shifted and truncated cutoff.
  double x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,x29,x30,x31,x32;
};




#endif // SLIT_PORE2__
