#ifndef __LUTSKO_EOS__
#define __LUTSKO_EOS__

#include <ostream>
#include "Enskog.h"

// Empirical equations of state.

const string LJ_JZG_EOS("LJ_JZG_EOS");
const string LJ_MECKE_EOS("LJ_MECKE_EOS");
const string HS_PY_EOS("HS_PY_EOS");
const string EOS_NULL("EOS_NULL");


class EOS
{
public:
  EOS(double kT): kT_(kT){}
  virtual ~EOS(){}

  // Basic functionality
  virtual double get_energy_per_atom(  double density) = 0;     // F/(NkT) = FE per unit atom (phi)
  virtual double get_energy_per_volume(double density) {return density*get_energy_per_atom(density);}
  virtual double get_pressure(double density) = 0;     // P/nkT

  // convenience helpers ...
  virtual double get_energy_per_atom(double density, double kT)   {kT_ = kT; return get_energy_per_atom(density);}
  virtual double get_energy_per_volume(double density, double kT) {          return density*get_energy_per_atom(density,kT);}
  virtual double get_pressure(double density, double kT)          {kT_ = kT; return get_pressure(density);}

  // Excess F_ex/(kT*V) and derivatives
  virtual double fex(double density)  {return density*phix(density);}
  virtual double f1ex(double density) {return phix(density) + density*phi1x(density);}
  virtual double f2ex(double density) {return 2*phi1x(density) + density*phi2x(density);}
  
  // excess free energy per atom and density derivatives
  virtual double phix(double density)  const = 0; //{ throw std::runtime_error("phix not implemented in EOS object");}
  virtual double phi1x(double density) const = 0; //{ throw std::runtime_error("phix not implemented in EOS object");}
  virtual double phi2x(double density) const = 0; //{ throw std::runtime_error("phix not implemented in EOS object");}
  virtual double phi3x(double density) const = 0; //{ throw std::runtime_error("phix not implemented in EOS object");}

  virtual bool isNull() {return false;}
  
  void set_temperature(double kT) {kT_ = kT;}
  
  virtual char const * const getName() const = 0;
  virtual double getEffHSD(){return 0.0;}
  /// Allow printing of the description
  friend std::ostream& operator<< (std::ostream& o, const EOS& ps)
  {
    o << "#EOS model is " << ps.getName() << std::endl;
    return o;
  }
  
protected:
  double kT_ = 1;
};

class EOS_NULL_ : public EOS
{
public:
  EOS_NULL_(double kT):  EOS(kT) {}
  virtual ~EOS_NULL_(){}

  // Basic functionality
  virtual double get_energy_per_atom(  double density) {return 0.0;}  // F/(NkT) = FE per unit atom (phi)
  virtual double get_energy_per_volume(double density) {return 0.0;}
  virtual double get_pressure(double density) {return 0.0;}     // P/nkT

  // convenience helpers ...
  virtual double get_energy_per_atom(double density, double kT)   {kT_ = kT; return get_energy_per_atom(density);}
  virtual double get_energy_per_volume(double density, double kT) {          return density*get_energy_per_atom(density,kT);}
  virtual double get_pressure(double density, double kT)          {kT_ = kT; return get_pressure(density);}

  // Excess F_ex/(kT*V) and derivatives
  virtual double fex(double density)  {return 0.0;} //density*phix(density);}
  virtual double f1ex(double density) {return 0.0;} //phix(density) + density*phi1x(density);}
  virtual double f2ex(double density) {return 0.0;} //2*phi1x(density) + density*phi2x(density);}
  
  // excess free energy per atom and density derivatives
  virtual double phix(double density)  const {return 0.0;} //{ throw std::runtime_error("phix not implemented in EOS object");}
  virtual double phi1x(double density) const {return 0.0;} //{ throw std::runtime_error("phix not implemented in EOS object");}
  virtual double phi2x(double density) const {return 0.0;} //{ throw std::runtime_error("phix not implemented in EOS object");}
  virtual double phi3x(double density) const {return 0.0;} //{ throw std::runtime_error("phix not implemented in EOS object");}

  void set_temperature(double kT) {kT_ = kT;}

  virtual bool isNull() {return true;}
  
  virtual char const * const getName() const {return "NULL";}
  virtual double getEffHSD(){return 0.0;}
  
protected:
  double kT_ = 1;
};

class PY_EOS : public EOS
{
 public:
  PY_EOS(double kT, double rc = -1.0, int noShift = 0)
   : EOS(kT) {}

  virtual ~PY_EOS() {}

  virtual char const * const getName() const { return "PY_EOS";}
  
  virtual double get_energy_per_atom(double density)
  {
    Enskog e(density);
    return e.freeEnergyPYC();
  }
  virtual double get_pressure(double density)
  {
    Enskog e(density);
    return e.pressurePYC();
  }
  // excess free energy per atom and density derivatives
  virtual double phix(double density)  const 
  {
    Enskog e(density);
    return e.exFreeEnergyPYC();
  }

  virtual double phi1x(double density) const 
  {
    Enskog e(density);
    return e.dexFreeEnergyPYCdRho();
  }

  virtual double phi2x(double density) const
  {
    Enskog e(density);
    return e.d2exFreeEnergyPYCdRho2();    
  }
  virtual double phi3x(double density) const 
  {
    Enskog e(density);
    return e.d3exFreeEnergyPYCdRho3();        
  }

 private:
};


class LJ_JZG : public EOS
{
 public:
  LJ_JZG(double kT, double rc = -1.0, int noShift = 0)
   : EOS(kT), gamma(3.0), da_(0.0)
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
	  if(noShift != 0) da_ += (8.0/3.0)*M_PI*(pow(rc,-9.0)-pow(rc,-3.0));
	}
    }

  virtual ~LJ_JZG() {}

  virtual char const * const getName() const { return "JZG";}
  
  virtual double get_energy_per_atom(double density)
  {
    double f  = kT_*(log(density)-1);
    for(int i=1;i<=8;i++) f += a(i)*pow(density,i)/i;
    for(int i=1;i<=6;i++) f += b(i)*G(density,i);
    f += da_*density;
    return f/kT_;
  }
  virtual double get_pressure(double density)
  {
    double P  = density*kT_;
    double F = exp(-gamma*density*density);
    for(int i=1;i<=8;i++) P += a(i)*pow(density,i+1);
    for(int i=1;i<=6;i++) P += F*b(i)*pow(density,2*i+1);
    P += da_*density*density;
    return P/(density*kT_);
  }
  // excess free energy per atom and density derivatives
  virtual double phix(double density)  const 
  { 
      double f  = 0.0;
      for(int i=1;i<=8;i++) f += a(i)*pow(density,i)/i;
      for(int i=1;i<=6;i++) f += b(i)*G(density,i);
      f += da_*density;
      return f/kT_;
  }

  virtual double phi1x(double density) const 
  { 
      double F = exp(-gamma*density*density);

      double f  = 0.0;
      for(int i=1;i<=8;i++) f += a(i)*pow(density,i-1);
      for(int i=1;i<=6;i++) f += b(i)*F*pow(density,2*i-1);
      f += da_;
      return f/kT_;
  }

  virtual double phi2x(double density) const
  { 
      double F = exp(-gamma*density*density);

      double f  = 0.0;
      for(int i=2;i<=8;i++) f += a(i)*(i-1)*pow(density,i-2);
      for(int i=1;i<=6;i++) f += b(i)*F*((2*i-1)*pow(density,2*i-2)-2*gamma*pow(density,2*i));
      return f/kT_;
  }
  virtual double phi3x(double density) const 
  { 
      double F = exp(-gamma*density*density);

      double f  = 0.0;
      for(int i=3;i<=8;i++) f += a(i)*(i-1)*(i-2)*pow(density,i-3);
      for(int i=1;i<=6;i++) f += b(i)*(-2*gamma*density)*F*((2*i-1)*pow(density,2*i-2)-2*gamma*pow(density,2*i));
      for(int i=2;i<=6;i++) f += b(i)*F*(2*i-1)*(2*i-2)*pow(density,2*i-3);
      for(int i=1;i<=6;i++) f += b(i)*F*(-2*2*i*gamma*pow(density,2*i-1));
      return f/kT_;
  }

 private:

  double a(int i) const
  {
    if(i == 1) return x1*kT_+x2*sqrt(kT_)+x3+(x4/kT_)+x5/(kT_*kT_);
    if(i == 2) return x6*kT_+x7+(x8/kT_)+x9/(kT_*kT_);
    if(i == 3) return x10*kT_+x11+(x12/kT_);
    if(i == 4) return x13;
    if(i == 5) return (x14/kT_)+x15/(kT_*kT_);
    if(i == 6) return x16/kT_;
    if(i == 7) return (x17/kT_)+x18/(kT_*kT_);
    if(i == 8) return x19/(kT_*kT_);
    throw std::runtime_error("Unknown values in LJ");
  }

  double b(int i) const 
  {
    double T1 = 1.0/kT_;
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



class LJ_Mecke : public EOS
{
 public:
  LJ_Mecke(double kT, double rc = -1.0, int noShift = 0)
   : EOS(kT), rhoc_(0.3107), kTc_(1.328), da_(0.0)
    {

      c_[0] =  0.33619760720e-05;
      c_[1] = -0.14707220591e+01;
      c_[2] = -0.11972121043e+00;
      c_[3] = -0.11350363539e-04;
      c_[4] = -0.26778688896e-04;

      c_[5] =  0.12755936511e-05;
      c_[6] =  0.40088615477e-02;
      c_[7] =  0.52305580273e-05;
      c_[8] = -0.10214454556e-07;
      c_[9] = -0.14526799362e-01;

      c_[10] =  0.64975356409e-01;
      c_[11] = -0.60304755494e-01;
      c_[12] = -0.14925537332e+00;
      c_[13] = -0.31664355868e-03;
      c_[14] =  0.28312781935e-01;

      c_[15] =  0.13039603845e-03;
      c_[16] =  0.10121435381e-01;
      c_[17] = -0.15425936014e-04;
      c_[18] = -0.61568007279e-01;
      c_[19] =  0.76001994423e-02;

      c_[20] = -0.18906040708e+00;
      c_[21] =  0.33141311846e+00;
      c_[22] = -0.25229604842e+00;
      c_[23] =  0.13145401812e+00;
      c_[24] = -0.48672350917e-01;

      c_[25] =  0.14756043863e-02;
      c_[26] = -0.85996667747e-02;
      c_[27] =  0.33880247915e-01;
      c_[28] =  0.69427495094e-02;
      c_[29] = -0.22271531045e-07;

      c_[30] = -0.22656880018e-03;
      c_[31] =  0.24056013779e-02;
      
      m_[0] = -2;
      m_[1] = -1;
      m_[2] = -1;
      m_[3] = -1 ;
      m_[4] = -0.5;
      
      m_[5] = -0.5;
      m_[6] = 0.5;
      m_[7] = 0.5;
      m_[8] = 1;
      m_[9] = -5;

      m_[10] = -4;
      m_[11] = -2;
      m_[12] = -2;
      m_[13] = -2;
      m_[14] = -1;

      m_[15] = -1;
      m_[16] = 0;
      m_[17] = 0;
      m_[18] = -5;
      m_[19] = -4;

      m_[20] = -3;
      m_[21] = -2;
      m_[22] = -2;
      m_[23] = -2;
      m_[24] = -1;

      m_[25] = -10;
      m_[26] = -6;
      m_[27] = -4;
      m_[28] = 0;
      m_[29] = -24;

      m_[30] = -10;
      m_[31] = -2;


      n_[0] = 9;
      n_[1] = 1;
      n_[2] = 2;
      n_[3] = 9;
      n_[4] = 8;
      
      n_[5] = 10;
      n_[6] = 1;
      n_[7] = 7;
      n_[8] = 10;
      n_[9] = 1;
      
      n_[10] = 1;
      n_[11] = 1;
      n_[12] = 2;
      n_[13] = 8;
      n_[14] = 1;

      n_[15] = 10;
      n_[16] = 4;
      n_[17] = 9;
      n_[18] = 2;
      n_[19] = 5;

      n_[20] = 1;
      n_[21] = 2;
      n_[22] = 3;
      n_[23] = 4;
      n_[24] = 2;

      n_[25] = 3;
      n_[26] = 4;
      n_[27] = 2;
      n_[28] = 2;
      n_[29] = 5;

      n_[30] = 2;
      n_[31] = 10;
      
      p_[0] = 0;
      p_[1] = 0;
      p_[2] = 0;
      p_[3] = 0;
      p_[4] = 0;

      p_[5] = 0;
      p_[6] = 0;
      p_[7] = 0;
      p_[8] = 0;
      p_[9] = -1;

      p_[10] = -1;
      p_[11] = -1;
      p_[12] = -1;
      p_[13] = -1;
      p_[14] = -1;

      p_[15] = -1;
      p_[16] = -1;
      p_[17] = -1;
      p_[18] = -1;
      p_[19] = -1;

      p_[20] = -1;
      p_[21] = -1;
      p_[22] = -1;
      p_[23] = -1;
      p_[24] = -1;

      p_[25] = -1;
      p_[26] = -1;
      p_[27] = -1;
      p_[28] = -1;
      p_[29] = -1;

      p_[30] = -1;
      p_[31] = -1;

      q_[0] = 0;
      q_[1] = 0;
      q_[2] = 0;
      q_[3] = 0;
      q_[4] = 0;

      q_[5] = 0;
      q_[6] = 0;
      q_[7] = 0;
      q_[8] = 0;
      q_[9] = 1;

      q_[10] = 1;
      q_[11] = 1;
      q_[12] = 1;
      q_[13] = 1;
      q_[14] = 1;

      q_[15] = 1;
      q_[16] = 1;
      q_[17] = 1;
      q_[18] = 2;
      q_[19] = 2;

      q_[20] = 2;
      q_[21] = 2;
      q_[22] = 2;
      q_[23] = 2;
      q_[24] = 2;

      q_[25] = 3;
      q_[26] = 3;
      q_[27] = 3;
      q_[28] = 3;
      q_[29] = 4;

      q_[30] = 4;
      q_[31] = 4;

      if(rc > 0.0)
	{
	  da_ = -(32.0/9.0)*M_PI*(pow(rc,-9.0)-1.5*pow(rc,-3.0));
	  if(noShift != 0) da_ += (8.0/3.0)*M_PI*(pow(rc,-9.0)-pow(rc,-3.0));
	}
    }

  virtual ~LJ_Mecke() {}

  virtual char const * const getName() const { return "LJMecke";}
  
  virtual double get_energy_per_atom(double rho)
  {
    double s = 0.1617*(rho/rhoc_)/(0.689+0.311*pow(kT_/kTc_,0.3674));
    double fid = log(rho)-1;
    double fhs = (4*s-3*s*s)*pow(1-s,-2);
    double fa  = 0.0;
    
    for(int i=0;i<32;i++)
      fa += c_[i]*pow(kT_/kTc_,m_[i])*pow(rho/rhoc_,n_[i])*exp(p_[i]*pow(rho/rhoc_,q_[i]));
    
    fa += da_*rho/kT_;

    return fid+fhs+fa;
  }
  // P/NkT
  virtual double get_pressure(double rho)
  {
    double s = 0.1617*(rho/rhoc_)/(0.689+0.311*pow(kT_/kTc_,0.3674));
    double pid = 1.0;
    double phs = (4*s-6*s*s)*pow(1-s,-2) + 2*s*(4*s-3*s*s)*pow(1-s,-3);
    double pa  = 0.0;
    
    for(int i=0;i<32;i++)
      pa += (n_[i]+p_[i]*q_[i]*pow(rho/rhoc_,q_[i]))
	*c_[i]*pow(kT_/kTc_,m_[i])*pow(rho/rhoc_,n_[i])*exp(p_[i]*pow(rho/rhoc_,q_[i]));

    pa += da_*rho/kT_;
    return pid+phs+pa;
  }

 private:
  double rhoc_;
  double kTc_;
  double da_;
  double c_[32];
  double m_[32];
  double n_[32];
  double p_[32];
  double q_[32];
};

class LJ_FMSA : public EOS
{
 public:
  LJ_FMSA(double kT, double rc = -1.0, int noShift = 0)
   : EOS(kT), da_(0.0)
    {
      if(rc > 0.0)
	{
	  da_ = -(32.0/9.0)*M_PI*(pow(rc,-9.0)-1.5*pow(rc,-3.0));
	  if(noShift != 0) da_ += (8.0/3.0)*M_PI*(pow(rc,-9.0)-pow(rc,-3.0));
	}
    }

  virtual ~LJ_FMSA() {}
  
  virtual char const * const getName() const { return "LJ_FMSA";}

  virtual double get_energy_per_atom(double rho)
  {

    double d = (1+0.2977*kT_)/(1+0.33163*kT_+1.0477*(1e-3)*kT_*kT_);

    double e = M_PI*rho*d*d*d/6;
    double phi_cs = log(rho)-1+e*(4-3*e)*pow(1-e,-2);

    double z1 = 2.9637;
    double z2 = 14.0167;

    double k1 = 2.1714*exp(z1*(1-d));
    double k2 = 2.1714*exp(z2*(1-d));

    double T1 = kT_*d/k1;
    double T2 = kT_*d/k2;

    double L1 = L_(z1*d,e);
    double L2 = L_(z2*d,e);
    
    double Q1 = Q_(z1*d,e);
    double Q2 = Q_(z2*d,e);
    
    double g0 = (1+e/2)*pow(1-e,-2);
    
    double phi_0 = -(k1/kT_)*L1*pow(1-e,-2)/(Q1*z1*z1);
    phi_0 -= -(k1/kT_)*(1+z1*d)/(z1*z1);
    phi_0 -=  -(k2/kT_)*L2*pow(1-e,-2)/(Q2*z2*z2);
    phi_0 += -(k2/kT_)*(1+z2*d)/(z2*z2);
    
    phi_0 += (4/kT_)*((1.0/9)*pow(d,-9)-(1.0/3)*pow(d,-3));
    phi_0 -= (4/kT_)*g0*((1.0/9)*pow(d,-9)-(1.0/3)*pow(d,-3)+(2.0/9));
        
    double phi_1 = -(k1*k1/(kT_*kT_))/(2*z1*Q1*Q1*Q1*Q1);
    phi_1 -= (k2*k2/(kT_*kT_))/(2*z2*Q2*Q2*Q2*Q2);
    phi_1 += 2*(k1*k2/(kT_*kT_))/((z2+z1)*Q1*Q1*Q2*Q2);
    
    phi_1 -= (4/kT_)*(d*d*d)*((1.0/(Q1*Q1*T1))-(1.0/(Q2*Q2*T2)))*((1.0/9)*pow(d,-9)-(1.0/3)*pow(d,-3)+(2.0/9));

    return phi_cs+4*M_PI*0.5*rho*(phi_0+0.5*phi_1) + da_*rho/kT_;
  }
  // P/NkT
  virtual double get_pressure(double rho)
  {
    double d = (1+0.2977*kT_)/(1+0.33163*kT_+1.0477*(1e-3)*kT_*kT_);
    
    double e = M_PI*rho*d*d*d/6;
    double P_cs = (1+e+e*e-e*e*e)*pow(1-e,-3);

    double z1 = 2.9637;
    double z2 = 14.0167;

    double k1 = 2.1714*exp(z1*(1-d));
    double k2 = 2.1714*exp(z2*(1-d));

    double T1 = kT_*d/k1;
    double T2 = kT_*d/k2;

    double L1 = L_(z1*d,e);
    double L2 = L_(z2*d,e);
    double L11 = L1_(z1*d,e);
    double L21 = L1_(z2*d,e);
    
    double Q1 = Q_(z1*d,e);
    double Q2 = Q_(z2*d,e);
    double Q11 = Q1_(z1*d,e);
    double Q21 = Q1_(z2*d,e);
    
    double g0 = (1+e/2)*pow(1-e,-2);
    double g01 = 0.5*pow(1-e,-2)+2*(1+e/2)*pow(1-e,-3);
    
    double phi_0 = 0.0;
    phi_0 -= (k1/kT_)*L1*pow(1-e,-2)/(Q1*z1*z1);
    phi_0 -= -(k1/kT_)*(1+z1*d)/(z1*z1);
    phi_0 -= -(k2/kT_)*L2*pow(1-e,-2)/(Q2*z2*z2);
    phi_0 += -(k2/kT_)*(1+z2*d)/(z2*z2);
    
    phi_0 += (4/kT_)*((1.0/9)*pow(d,-9)-(1.0/3)*pow(d,-3));
    phi_0 -= (4/kT_)*g0*((1.0/9)*pow(d,-9)-(1.0/3)*pow(d,-3)+(2.0/9));

    double phi_01 = 0;
    phi_01 -= (k1/kT_)*L11*pow(1-e,-2)/(Q1*z1*z1);
    phi_01 -= 2*(k1/kT_)*L1*pow(1-e,-3)/(Q1*z1*z1);
    phi_01 -= -Q11*(k1/kT_)*L1*pow(1-e,-2)/(Q1*Q1*z1*z1);

    phi_01 -= -(k2/kT_)*L21*pow(1-e,-2)/(Q2*z2*z2);
    phi_01 -= -2*(k2/kT_)*L2*pow(1-e,-3)/(Q2*z2*z2);
    phi_01 -= Q21*(k2/kT_)*L2*pow(1-e,-2)/(Q2*Q2*z2*z2);

    phi_01 -= (4/kT_)*g01*((1.0/9)*pow(d,-9)-(1.0/3)*pow(d,-3)+(2.0/9));
    
    double phi_1 = 0.0;
    phi_1 -= (k1*k1/(kT_*kT_))/(2*z1*Q1*Q1*Q1*Q1);
    phi_1 -= (k2*k2/(kT_*kT_))/(2*z2*Q2*Q2*Q2*Q2);
    phi_1 += 2*(k1*k2/(kT_*kT_))/((z2+z1)*Q1*Q1*Q2*Q2);
    
    phi_1 -= (4/kT_)*(d*d*d)*((1.0/(Q1*Q1*T1))-(1.0/(Q2*Q2*T2)))*((1.0/9)*pow(d,-9)-(1.0/3)*pow(d,-3)+(2.0/9));

    double phi_11 = 0.0;
    phi_11 -= -4*Q11*(k1*k1/(kT_*kT_))/(2*z1*Q1*Q1*Q1*Q1*Q1);
    phi_11 -= -4*Q21*(k2*k2/(kT_*kT_))/(2*z2*Q2*Q2*Q2*Q2*Q2);
    phi_11 += -2*Q11*2*(k1*k2/(kT_*kT_))/((z2+z1)*Q1*Q1*Q1*Q2*Q2);
    phi_11 += -2*Q21*2*(k1*k2/(kT_*kT_))/((z2+z1)*Q1*Q1*Q2*Q2*Q2);
    
    phi_11 -= (4/kT_)*(d*d*d)*((-2*Q11/(Q1*Q1*Q1*T1))-(-2*Q21/(Q2*Q2*Q2*T2)))*((1.0/9)*pow(d,-9)-(1.0/3)*pow(d,-3)+(2.0/9));

    return P_cs+4*M_PI*0.5*rho*(phi_0+0.5*phi_1)+4*M_PI*0.5*rho*e*(phi_01+0.5*phi_11) + da_*rho/kT_;
  }

 private:
  double L_(double z, double e) { return (1+e/2)*z+1+2*e;}
  double L1_(double z, double e) { return 0.5*z+2;}
  double S_(double z, double e) { return (1-e)*(1-e)*z*z*z+6*e*(1-e)*z*z+18*e*e*z-12*e*(1+2*e);}
  double S1_(double z, double e) { return -2*(1-e)*z*z*z+6*(1-2*e)*z*z+18*2*e*z-12*(1+4*e);}

  double Q_(double z, double e)
  {
    double L = L_(z,e);
    double S = S_(z,e);
    return (S+12*e*L*exp(-z))*pow(1-e,-2)/(z*z*z);
  }

  double Q1_(double z, double e)
  {
    double L = L_(z,e);
    double L1 = L1_(z,e);
    double S = S_(z,e);
    double S1 = S1_(z,e);
    return ((S1+12*L*exp(-z)+12*e*L1*exp(-z))*pow(1-e,-2)/(z*z*z))
      +2*((S+12*e*L*exp(-z))*pow(1-e,-3)/(z*z*z));
  }

 private:
  double da_;

};



#endif // __LUTSKO_EOS__
