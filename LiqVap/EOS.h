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
   //   double dP_dxdT(double x, double kT)  const { return 1.0+2*x*dphix_dxdT(x,kT) + x*x*dphix_dx2dT(x,kT);}
   double dP_dx2(double x, double kT)   const { return 2*dphix_dx(x,kT)+4*x*dphix_dx2(x,kT)+x*x*dphix_dx3(x,kT);}
   //   double dP_dx2dT(double x, double kT) const { return 2*dphix_dxdT(x,kT)+4*x*dphix_dx2dT(x,kT)+x*x*dphix_dx3dT(x,kT);}
   
   // excess free energy per atom:  F/N
   virtual double phix(double density, double kT)  const = 0;

   // derivative of excess free energy per atom :  d(F/N)/dx    
   virtual double dphix_dx(double density, double kT) const = 0;

   // 2nd derivative of excess free energy per atom :  d2(F/N)/dx2       
   virtual double dphix_dx2(double density, double kT) const = 0;

   // 3rd derivative of excess free energy per atom :  d3(F/N)/d32       
   virtual double dphix_dx3(double density, double kT) const = 0;

   //   virtual double dphix_dxdT(double density, double kT) const = 0;
   //   virtual double dphix_dx2dT(double density, double kT) const = 0;
   //   virtual double dphix_dx3dT(double density, double kT) const = 0;

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

class WRDF_Colloid : public EOS
{
 public:
 WRDF_Colloid() : EOS()
    {
      a_[0][0] = -1.10009;
      a_[0][1] = 6.27642;
      a_[0][2] = -13.5539;
      a_[0][3] = 15.5901;
      a_[0][4] = -7.05497;
      a_[0][5] = 0.63096;

      a_[1][0] = 4.79833;
      a_[1][1] = -31.7208;
      a_[1][2] = 76.0381;
      a_[1][3] = -84.1657;
      a_[1][4] = 46.4876;
      a_[1][5] = -10.4596;      

      a_[2][0] = -0.47723;
      a_[2][1] = 19.8875;
      a_[2][2] = -62.2031;
      a_[2][3] = 78.4142;
      a_[2][4] = -55.7024;
      a_[2][5] = 18.4785;

      a_[3][0] = -7.60729;
      a_[3][1] = 39.2723;
      a_[3][2] = -131.249;
      a_[3][3] = 205.539;
      a_[3][4] = -110.462;
      a_[3][5] = 9.10103;

      a_[4][0] = -22.4105;
      a_[4][1] = 65.7709;
      a_[4][2] = 66.3596;
      a_[4][3] = -315.674;
      a_[4][4] = 244.739;
      a_[4][5] = -43.4275;

      a_[0][0] = 56.4536;
      a_[5][1] = -228.002;
      a_[5][2] = 243.158;
      a_[5][3] = 30.6737;
      a_[5][4] = -133.904;
      a_[5][5] = 33.7645;

      a_[6][0] = -29.5309;
      a_[6][1] = 127.279;
      a_[6][2] = -174.625;
      a_[6][3] = 68.3272;
      a_[6][4] = 18.5640;
      a_[6][5] = -9.79081;
    }

  virtual ~WRDF_Colloid() {}

  virtual char const * const getName() const { return "WRDF_Colloid";}

  // excess free energy per atom and density derivatives
  virtual double phix(double density, double kT)  const 
  { 
      double a  = 0.0;

      for(int n=2; n<=8;n++)
	for(int m=-3; m<=2; m++)
	  a += kT*a_[n-2][m+3]*pow(density,n-1)*pow(kT,-m);
      return a;
  }

  virtual double dphix_dx(double density, double kT) const 
  {
      double a  = 0.0;

      for(int n=2; n<=8;n++)
	for(int m=-3; m<=2; m++)
	  a += kT*a_[n-2][m+3]*(n-1)*(n-2 < 0 ? 0 : (n-2 == 0 ? 1 : pow(density,n-2)))*pow(kT,-m);
      return a;
  }

  virtual double dphix_dx2(double density, double kT) const
  {
      double a  = 0.0;

      for(int n=2; n<=8;n++)
	for(int m=-3; m<=2; m++)
	  a += kT*a_[n-2][m+3]*(n-1)*(n-2)*(n-3 < 0 ? 0 : (n-3 == 0 ? 1 : pow(density,n-3)))*pow(kT,-m);
      return a;    
  }
  virtual double dphix_dx3(double density, double kT) const 
  {
      double a  = 0.0;

      for(int n=2; n<=8;n++)
	for(int m=-3; m<=2; m++)
	  a += kT*a_[n-2][m+3]*(n-1)*(n-2)*(n-3)*(n-3 < 0 ? 0 : (n-3 == 0 ? 1 : pow(density,n-4)))*pow(kT,-m);
      return a;        
  }

  private:
  double a_[7][6];
};


class WRDF_LJ : public EOS
{
 public:
 WRDF_LJ() : EOS()
    {
      // ranges are actually n = 2, ... 8; m = -3 ... 2
      a_[0][0] =  7.690224838392756279e-01;
      a_[0][1] = -3.153910169943185959e+00;
      a_[0][2] =  2.886406073436768693e+00;
      a_[0][3] =  5.170749028341210085e+00;
      a_[0][4] = -1.044622980663629441e+01;
      a_[0][5] =  8.656382956332409062e-01;

      a_[1][0] =  3.269747272846510011e+00;
      a_[1][1] = -7.773917981058953330e+01;
      a_[1][2] =  3.685941276727217542e+02;
      a_[1][3] = -7.230317930487650528e+02;
      a_[1][4] =  6.465538365485742816e+02;
      a_[1][5] = -2.193018472879539331e+02;

      a_[2][0] = -8.293526173216596575e+00;
      a_[2][1] =  3.778727621371398868e+02;
      a_[2][2] = -1.929641574600540707e+03;
      a_[2][3] =  3.888378430601394939e+03;
      a_[2][4] = -3.522436436983288786e+03;
      a_[2][5] =  1.219511341403025654e+03;

      a_[3][0] = -3.188581072468004507e+01;
      a_[3][1] = -5.978745233526092306e+02;
      a_[3][2] =  3.874727611676720244e+03;
      a_[3][3] = -8.352099111519108192e+03;
      a_[3][4] =  7.764536367165250340e+03;
      a_[3][5] = -2.720454410705084683e+03;

      a_[4][0] =  1.022571045500056499e+02;
      a_[4][1] =  3.492984522056419223e+02;
      a_[4][2] = -3.783550535487431262e+03;
      a_[4][3] =  8.937046874892746928e+03;
      a_[4][4] = -8.548363392920171464e+03;
      a_[4][5] =  3.013089322044620531e+03;
      
      a_[5][0] = -9.717532400121930891e+01;
      a_[5][1] = -4.109700355637258440e+00;
      a_[5][2] =  1.783386606546049052e+03;
      a_[5][3] = -4.738687139283345459e+03;
      a_[5][4] =  4.672478415816356573e+03;
      a_[5][5] = -1.647794570675797786e+03;

      a_[6][0] =  3.112982250741344359e+01;
      a_[6][1] = -4.474045918417971279e+01;
      a_[6][2] = -3.155759031266126726e+02;
      a_[6][3] =  9.856584906924986171e+02;
      a_[6][4] = -1.004955570487645673e+03;
      a_[6][5] =  3.532964212447155319e+02;      
    }

  virtual ~WRDF_LJ() {}

  virtual char const * const getName() const { return "WRDF_LJ";}

  // excess free energy per atom and density derivatives
  virtual double phix(double density, double kT)  const 
  { 
      double a  = 0.0;

      for(int n=2; n<=8;n++)
	for(int m=-3; m<=2; m++)
	  a += kT*a_[n-2][m+3]*pow(density,n-1)*pow(kT,-m);
      return a;
  }

  virtual double dphix_dx(double density, double kT) const 
  {
      double a  = 0.0;

      for(int n=2; n<=8;n++)
	for(int m=-3; m<=2; m++)
	  a += kT*a_[n-2][m+3]*(n-1)*(n-2 < 0 ? 0 : (n-2 == 0 ? 1 : pow(density,n-2)))*pow(kT,-m);
      return a;
  }

  virtual double dphix_dx2(double density, double kT) const
  {
      double a  = 0.0;

      for(int n=2; n<=8;n++)
	for(int m=-3; m<=2; m++)
	  a += kT*a_[n-2][m+3]*(n-1)*(n-2)*(n-3 < 0 ? 0 : (n-3 == 0 ? 1 : pow(density,n-3)))*pow(kT,-m);
      return a;    
  }
  virtual double dphix_dx3(double density, double kT) const 
  {
      double a  = 0.0;

      for(int n=2; n<=8;n++)
	for(int m=-3; m<=2; m++)
	  a += kT*a_[n-2][m+3]*(n-1)*(n-2)*(n-3)*(n-3 < 0 ? 0 : (n-3 == 0 ? 1 : pow(density,n-4)))*pow(kT,-m);
      return a;        
  }

  private:
  double a_[7][6];
};




#endif // SLIT_PORE2__
