#ifndef __LUTSKO_WEIGHTED_DENSITIES__
#define __LUTSKO_WEIGHTED_DENSITIES__

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <complex>

#include "Species.h"
#include "Density.h"
#include "Point.h"

/**
  *  @brief FMT Class
  *
  *    This class holds two 11-dimentional arrays called d0_ and d1_. These in turn hold the weighted densities as
  *             d0_ = {eta(N),s(N),V1(N),...,VD(N), T11(N), T12(N), ... T1D(N), T22(N), ..., TDD(N)}
  *             where for D=3 there are 1+1+3+(3+2+1) = 11 entries.
  *             Note that each of these, e.g. eta(N), is an N-dimensional vector holding the value of the weighted density at each point.
  *             The two copies d0_ and d1_ are to allow for two species but this is not fully implemented yet. 
  */  


class FMT
{
 public:
  /**
  *   @brief  Default  constructor for FMT 
  *  
  *   @param  lattice is the Lattice object 
  *   @return nothing 
  */  
 FMT()  : etaMax_(1e30){} 
  /**
  *   @brief  Default  destrctur for FMT 
  *  
  *   @return nothing 
  */  
  ~FMT(){}

  virtual double BulkMuex(const vector<double> &x, const vector<Species*> &allSpecies, int species) const = 0;
  virtual double BulkFex(const vector<double> &x, const vector<Species*> &allSpecies) const = 0;
  
  /**
  *   @brief  EtaMax is the maximum value of eta for which the f1,f2_,f3_ functions are calculated "honestly". For larger values of eta, 
  *           the original, divergent forms are replaced by non-divergent quadratic forms. This is intended to allow for calculations without backtracking. It is only implemented for some models. 
  *  
  *   @param  etaMax is the threshold beyond which quadratic forms are used
  *   @return nothing 
  */    
  void setEtaMax(double etaMax) { etaMax_ = etaMax;}

  /**
  *   @brief  Calculates total free energy and dOmega/dRho(i) for each lattice point using FFT convolutions
  *  
  *   @param  density is current density
  *   @param  mu is the chemical potential
  *   @return total free energy of system
  */  
  double calculateFreeEnergyAndDerivatives(vector<Species*> &allSpecies);

  double doFreeEnergyLoop(vector<Species*> &allSpecies);
  
  /**
  *   @brief  Calculate phi(i) and some derivative information. 
  *
  *    This function assumes the weighted densities at each lattice site are available.
  *           It calculates the function PHI(r) at lattice site i and also stores dPHI(i)/deta(i), dPHI/dv(i), etc. 
  *  
  *   @param  i is lattice position
  *   @return phi(i)
  */  
  double dPHI(long i, vector<Species*> &allSpecies);


 /**
  *   @brief  calculate f1 cofactor in FMT PHI function.  This is the same for all models (log(1-eta)) and so is instantiated here.
  *  
  *   @param  eta is the local volume-average density
  *   @return f1(eta).
  */    
 virtual double f1_(double eta)  const
 {
    if(eta > etaMax_)
      {
	double a0 = log(1.0-etaMax_);
	double a1 = -1.0/(1.0-etaMax_);
	double a2 = -a1*a1;
	return a0 + a1*(eta-etaMax_) + 0.5*a2*(eta-etaMax_)*(eta-etaMax_);
      }
    
    return log(1.0-eta);
 }
  
  /**
  *   @brief  calculate f2 cofactor in FMT PHI function.  Different derived classes give different models for this funciton.
  *  
  *   @param  eta is the local volume-average density
  *   @return f2(eta).
  */    
 virtual double f2_(double eta)  const = 0;

 /**
  *   @brief  calculate f3 cofactor in FMT PHI function.  Different derived classes give different models for this funciton.
  *  
  *   @param  eta is the local volume-average density
  *   @return f3(eta).
  */    
 virtual double f3_(double eta)  const= 0;

  /**
  *   @brief  derivative of f1_() with respect to its argument.
  *  
  *   @param  eta is the local volume-average density
  *   @return f1(eta).
  */     
 virtual double f1p_(double eta) const
 {
    if(eta > etaMax_)
      {
	double a1 = -1.0/(1.0-etaMax_);
	double a2 = -a1*a1;
	return a1 + a2*(eta-etaMax_);
      }
    
    return -1.0/(1.0-eta);
 }
 
 /**
  *   @brief  derivative of f2_() with respect to its argument.
  *  
  *   @param  eta is the local volume-average density
  *   @return f2(eta).
  */     
 virtual double f2p_(double eta) const = 0;

 /**
  *   @brief  derivative of f3_() with respect to its argument.
  *  
  *   @param  eta is the local volume-average density
  *   @return f3(eta).
  */     
 virtual double f3p_(double eta) const= 0;

  /**
  *   @brief  calculates third contribution to PHI = PHI1 + PHI2 + PHI3. This is model dependent
  *  
  *   @param  s2 one of the scalar weighted densities (shell-averaged density)
  *   @param  v2_v2 is v2 dot v2
  *   @param  vTv is v dot T dot v
  *   @param  T2 is trace(T dot T)
  *   @param  T3 is trace(T dot T dot T)
  *   @return PHI3
  */     
 virtual double Phi3(double s2, double v2_v2, double vTv, double T2, double T3) const = 0;

 /**
  *   @brief  calculates derivative of third contribution to PHI = PHI1 + PHI2 + PHI3 with respect to S2
  *  
  *   @param  s2 one of the scalar weighted densities (shell-averaged density)
  *   @param  v2_v2 is v2 dot v2
  *   @param  vTv is v dot T dot v
  *   @param  T2 is trace(T dot T)
  *   @param  T3 is trace(T dot T dot T)
  *   @return dPHI3/dS2
  */     
 virtual double dPhi3_dS2(double s2,double v2_v2,double vTv,double T2,double T3) const = 0;

 /**
  *   @brief  calculates derivative of third contribution to PHI = PHI1 + PHI2 + PHI3 with respect to v2(k)
  *  
  *   @param  k is the index of v2
  *   @param  s2 one of the scalar weighted densities (shell-averaged density)
  *   @param  v2_v2 is v2 dot v2
  *   @param  vTv is v dot T dot v
  *   @param  T2 is trace(T dot T)
  *   @param  T3 is trace(T dot T dot T)
  *   @return dPHI3/dv2(k)
  */      
 virtual double dPhi3_dV2(int k, double s2, double v2_v2, double v2[], double vT[]) const = 0;

 /**
  *   @brief  calculates derivative of third contribution to PHI = PHI1 + PHI2 + PHI3 with respect to T(j,k)
  *  
  *   @param  j is the first  index of T
  *   @param  k is the second index of T
  *   @param  s2 one of the scalar weighted densities (shell-averaged density)
  *   @param  v2_v2 is v2 dot v2
  *   @param  vTv is v dot T dot v
  *   @param  T2 is trace(T dot T)
  *   @param  T3 is trace(T dot T dot T)
  *   @return dPHI3/dT(j,k)
  */       
 virtual double dPhi3_dT(int j,int k,double s2, double v2[], double T0[3][3], double TT[3][3]) const = 0;


 virtual bool needsTensor() const { return true;}
 
 protected:

/**
  *   @brief  This function just sums up phi(i) over the lattice to 
  *           give the total free energy
  *
  *   @param  density is the Density object
  *   @return The total free energy
  */  
 double calculateFreeEnergy(vector<Species*> &allSpecies);

 /**
  *   @brief  name of model implemented by this class
  *
  *   @return the name as a string
  */        
 virtual string Name() const = 0;

 protected:
 
 double etaMax_; ///< cutoff used to control divergences
 
};

/**
  *  @brief  The original White Bear  FMT model 
  *
  *   White-bear FMT model with tensor densities and the CS equation of state
  */  
class WhiteBearI : public FMT
{
 public:
  WhiteBearI() : FMT(){};

   virtual bool needsTensor() const { return true;}


  virtual double f2_(double eta) const
  {
    return 1.0/(1.0-eta);
  }

  virtual double f2p_(double eta) const
  {
    double f = 1.0/(1.0-eta);
    return f*f;
  }

  virtual double f3_(double eta) const
  {
    if(eta < 1e-12)
      return 1.5+(8.0/3)*eta+3.75*eta*eta+4.8*eta*eta*eta;
    return (1.0/(eta*(1-eta)*(1-eta))) +  (log(1-eta)/(eta*eta));
  }

  virtual double f3p_(double eta) const
  {
    if(eta < 1e-12)
      return (8.0/3)+7.5*eta+14.4*eta*eta;
    return (1.0/(eta*eta))*((2/((1-eta)*(1-eta)*(1-eta)))-(3/((1-eta)*(1-eta)))-(1.0/(1-eta)))-(2/(eta*eta*eta))*log(1-eta);
  }

  virtual double Phi3(double s2, double v2_v2, double vTv, double T2, double T3) const
  {
    return (1.0/(8*M_PI))*(vTv-s2*v2_v2-T3+s2*T2);
  }

 virtual double dPhi3_dS2(double s2,double v2_v2,double vTv,double T2,double T3) const 
 { 
   return (1.0/(8*M_PI))*(-v2_v2+T2);
 }
 virtual double dPhi3_dV2(int k, double s2, double v2_v2, double v2[], double vT[]) const 
 { 
   return (1.0/(8*M_PI))*(2*vT[k]-2*s2*v2[k]);
 }

 virtual double dPhi3_dT(int j,int k,double s2, double v2[], double T0[3][3], double TT[3][3]) const 
 { 
   return (1.0/(8*M_PI))*(v2[j]*v2[k]-3*TT[j][k]+2*s2*T0[j][k]);
 }

 virtual double BulkMuex(const vector<double> &x, const vector<Species*> &allSpecies, int species) const
 {
   double n0 = 0.0;
   double n1 = 0.0;
   double n2 = 0.0;
   double n3 = 0.0;
   for(int i=0;i<x.size();i++)
     {
       double density = x[i];
       double hsd = allSpecies[i]->getHSD();
       n0 += density;
       n1 += 0.5*hsd*density;
       n2 += M_PI*hsd*hsd*density;
       n3 += (M_PI/6)*hsd*hsd*hsd*density;
     }

   double dPhi_dn0 = -log(1.0-n3);
   double dPhi_dn1 = n2/(1-n3);
   double dPhi_dn2 = (1.0/(12*M_PI))*(n2*n2/(n3*n3))*log(1-n3)+(n1/(1-n3))+n2*n2/(12*M_PI*n3*(1-n3)*(1-n3));
   double dPhi_dn3 = 0.0;

   dPhi_dn3  = -(1.0/(18*M_PI))*(n2*n2*n2/(n3*n3*n3))*log(1-n3);
   dPhi_dn3 += -((1/(36*M_PI))*(n2*n2*n2/(n3*n3))-n0)*(1.0/(1-n3));
   dPhi_dn3 += (n1*n2/((1-n3)*(1-n3)));
   dPhi_dn3 +=  (1.0/(36*M_PI))*((n2*n2*n2*(3*n3-1)/(n3*n3*(1-n3)*(1-n3)*(1-n3))));   


   double density = x[species];
   double hsd = allSpecies[species]->getHSD();

   return dPhi_dn0+dPhi_dn1*0.5*hsd+dPhi_dn2*M_PI*hsd*hsd+dPhi_dn3*(M_PI/6)*hsd*hsd*hsd; 
 }

 virtual double BulkFex(const vector<double> &x, const vector<Species*> &allSpecies) const
 {
   double n0 = 0.0;
   double n1 = 0.0;
   double n2 = 0.0;
   double n3 = 0.0;
   for(int i=0;i<x.size();i++)
     {
       double density = x[i];
       double hsd = allSpecies[i]->getHSD();
       n0 += density;
       n1 += 0.5*hsd*density;
       n2 += M_PI*hsd*hsd*density;
       n3 += (M_PI/6)*hsd*hsd*hsd*density;
     }

   double F = ((1.0/(36*M_PI))*(n2*n2*n2/(n3*n3))-n0)*log(1-n3);
   F += n1*n2/(1-n3);
   F += (1.0/(36*M_PI))*n2*n2*n2/(n3*(1-n3)*(1-n3));
		  
   return F;
 }

 
 virtual string Name() const { return string("WhiteBear I");}
};

/**
  *  @brief  The RSLT positive-definite FMT model (pre-White Bear)
  *
  *   This class implements the RSLT FMT model which is not as accurate as WhiteBear but which is positive definite and therefore stable.
  */  
class RSLT : public FMT
{
 public:
  RSLT() : FMT(){};


 virtual bool needsTensor() const { return false;}
  
  virtual double f2_(double eta) const
  {
    if(eta > etaMax_)
      {
	double a0 = 1.0/(1.0-etaMax_);
	double a1 = a0*a0;
	double a2 = 2*a1*a0;
	return a0 + a1*(eta-etaMax_) + 0.5*a2*(eta-etaMax_)*(eta-etaMax_);
      }
    
    return 1.0/(1.0-eta);
  }

  virtual double BulkMuex(const vector<double> &x, const vector<Species*> &allSpecies, int species) const
  {
    double n0 = 0.0;
    double n1 = 0.0;
    double n2 = 0.0;
    double n3 = 0.0;
    for(int i=0;i<x.size();i++)
      {
	double density = x[i];
	double hsd = allSpecies[i]->getHSD();
	n0 += density;
	n1 += 0.5*hsd*density;
	n2 += M_PI*hsd*hsd*density;
	n3 += (M_PI/6)*hsd*hsd*hsd*density;
      }

    double dPhi_dn0 = -log(1.0-n3);
    double dPhi_dn1 = n2/(1-n3);
    double dPhi_dn2 = (1.0/(12*M_PI))*(n2*n2/(n3*n3))*log(1-n3)+(n1/(1-n3))+n2*n2/(12*M_PI*n3*(1-n3)*(1-n3));
    double dPhi_dn3 = 0.0;

    dPhi_dn3  = -(1.0/(18*M_PI))*(n2*n2*n2/(n3*n3*n3))*log(1-n3);
    dPhi_dn3 += -((1/(36*M_PI))*(n2*n2*n2/(n3*n3))-n0)*(1.0/(1-n3));
    dPhi_dn3 += (n1*n2/((1-n3)*(1-n3)));
    dPhi_dn3 +=  (1.0/(36*M_PI))*((n2*n2*n2*(3*n3-1)/(n3*n3*(1-n3)*(1-n3)*(1-n3))));   

    double density = x[species];
    double hsd = allSpecies[species]->getHSD();

    return dPhi_dn0+dPhi_dn1*0.5*hsd+dPhi_dn2*M_PI*hsd*hsd+dPhi_dn3*(M_PI/6)*hsd*hsd*hsd; 
  }

  virtual double BulkFex(const vector<double> &x, const vector<Species*> &allSpecies) const
  {
    double n0 = 0.0;
    double n1 = 0.0;
    double n2 = 0.0;
    double n3 = 0.0;
    for(int i=0;i<x.size();i++)
      {
	double density = x[i];
	double hsd = allSpecies[i]->getHSD();
	n0 += density;
	n1 += 0.5*hsd*density;
	n2 += M_PI*hsd*hsd*density;
	n3 += (M_PI/6)*hsd*hsd*hsd*density;
      }

    double F = ((1.0/(36*M_PI))*(n2*n2*n2/(n3*n3))-n0)*log(1-n3);
    F += n1*n2/(1-n3);
    F += (1.0/(36*M_PI))*n2*n2*n2/(n3*(1-n3)*(1-n3));
		  
    return F;
  }
  

  virtual double f2p_(double eta) const
  {
    if(eta > etaMax_)
      {
	double a0 = 1.0/(1.0-etaMax_);	
	double a1 = a0*a0;
	double a2 = 2*a1*a0;
	return a1 + a2*(eta-etaMax_);
      }
    
    double f = 1.0/(1.0-eta);
    return f*f;
  }

  virtual double f3_(double eta) const
  {
    if(eta < 1e-12)
      return 1.5+(8.0/3)*eta+3.75*eta*eta+4.8*eta*eta*eta;

    if(eta > etaMax_)
      {
	double x = etaMax_;
	double a0 = (1.0/(x*(1-x)*(1-x))) +  (log(1-x)/(x*x));
	double a1 = -(1.0/(x*x))*(1.0/((1.0-x)*(1.0-x)*(1.0-x)))*(x*x-5*x+2)-(2.0/(x*x*x))*log(1-x);
	double a2 = -(1.0/(x*x*x))*(1.0/((1.0-x)*(1.0-x)*(1.0-x)*(1.0-x)))*(5*x*x*x-26*x*x+21*x-6)
	  +(6.0/(x*x*x*x))*log(1-x);
	return a0 + a1*(eta-etaMax_) + 0.5*a2*(eta-etaMax_)*(eta-etaMax_);
      }    
    return (1.0/(eta*(1-eta)*(1-eta))) +  (log(1-eta)/(eta*eta));
  }

  virtual double f3p_(double eta) const
  {
    if(eta < 1e-12)
      return (8.0/3)+7.5*eta+14.4*eta*eta;

    if(eta > etaMax_)
      {
	double x = etaMax_;
	double a1 = -(1.0/(x*x))*(1.0/((1.0-x)*(1.0-x)*(1.0-x)))*(x*x-5*x+2)-(2.0/(x*x*x))*log(1-x);
	double a2 = -(1.0/(x*x*x))*(1.0/((1.0-x)*(1.0-x)*(1.0-x)*(1.0-x)))*(5*x*x*x-26*x*x+21*x-6)
	  +(6.0/(x*x*x*x))*log(1-x);	
	return a1 + a2*(eta-etaMax_);
      }
    
    //    return (1.0/(eta*eta))*((2/((1-eta)*(1-eta)*(1-eta)))-(3/((1-eta)*(1-eta)))-(1.0/(1-eta)))-(2/(eta*eta*eta))*log(1-eta);
    return -(1.0/(eta*eta))*(1.0/((1.0-eta)*(1.0-eta)*(1.0-eta)))*(eta*eta-5*eta+2)-(2.0/(eta*eta*eta))*log(1-eta);
  }

  virtual double Phi3(double s2, double v2_v2, double vTv, double T2, double T3) const
  {
    if(s2 < 1e-20) return 0.0;
    double psi = v2_v2/(s2*s2);
    return (1.0/(36*M_PI))*s2*s2*s2*(1-psi)*(1-psi)*(1-psi);
  }

  virtual double dPhi3_dS2(double s2,double v2_v2,double vTv,double T2,double T3) const 
  { 
    if(s2 < 1e-20) return 0.0;
    double psi = v2_v2/(s2*s2);
    return (1.0/(12*M_PI))*s2*s2*(1+psi)*(1-psi)*(1-psi);
  }

  virtual double dPhi3_dV2(int k, double s2, double v2_v2, double v2[], double vT[]) const 
  { 
    if(s2 < 1e-20) return 0.0;
    double psi = v2_v2/(s2*s2);
    return -(1.0/(6*M_PI))*s2*v2[k]*(1-psi)*(1-psi);
  }

  virtual double dPhi3_dT(int j,int k,double s2, double v2[], double T0[3][3], double TT[3][3]) const 
  { 
    return 0;
  }

  virtual string Name() const { return string("RSLT");}
};

/**
  *  @brief  Modified RSLT positive-definite FMT model (pre-White Bear): for experimental purposes only!!
  *
  *   This class implements a modificaiton of the RSLT FMT model where I use the same basic idea to create a different type of model. In RSLT, the Rosenfeld numerator of s^3-3sv^2 (which is not positive definite) is replaced by s^3(1-v^2/s^2)^3 which agrees in the first two terms and is positive. This in fact works for  s^3(1-(3/n)v^2/s^2)^n for any n>=3. I used this class to experiment with different n but this was for research purposes only and should not be used without a thorough review ...  
  */  


class RSLT2: public RSLT
{
 public:
  RSLT2() : RSLT(){};

  virtual bool needsTensor() const { return false;}

  virtual double Phi3(double s2, double v2_v2, double vTv, double T2, double T3) const
  {
    double psi = v2_v2/(s2*s2);
    return (1.0/(36*M_PI))*s2*s2*s2*(1-1.5*psi)*(1-1.5*psi);
  }

 virtual double dPhi3_dS2(double s2,double v2_v2,double vTv,double T2,double T3) const 
 { 
   double psi = v2_v2/(s2*s2);
   return (1.0/(36*M_PI))*3*s2*s2*(1-1.5*psi)*(1-1.5*psi)
     +(1.0/(36*M_PI))*s2*s2*s2*2*(1-3*psi/2)*(-1.5)*(-2*psi/s2);
 }

 virtual double dPhi3_dV2(int k, double s2, double v2_v2, double v2[], double vT[]) const 
 { 
    double psi = v2_v2/(s2*s2);
    return (1.0/(36*M_PI))*s2*2*(1-3*psi/2)*(-1.5)*2*v2[k];
 }

 virtual double dPhi3_dT(int j,int k,double s2, double v2[], double T0[3][3], double TT[3][3]) const 
 { 
   return 0;
 }

 virtual string Name() const { return string("RSLT2");}

};
/**
  *  @brief  Rosenfeld modified to give CS EOS.
  *
  */  


class Rosenfeld: public RSLT
{
 public:
  Rosenfeld() : RSLT(){};
  /*  
  virtual double f3_(double eta) const
  {
    return 1.5/((1-eta)*(1-eta));
  }

  virtual double f3p_(double eta) const
  {
    return 1.5/((1-eta)*(1-eta)*(1-eta));
  }
  */

 virtual bool needsTensor() const { return false;}  
  
  virtual double Phi3(double s2, double v2_v2, double vTv, double T2, double T3) const
  {
    return (1.0/(36*M_PI))*(s2*s2*s2-3*s2*v2_v2);
  }

 virtual double dPhi3_dS2(double s2,double v2_v2,double vTv,double T2,double T3) const 
 { 
   return (1.0/(36*M_PI))*(3*s2*s2-3*v2_v2);
 }

 virtual double dPhi3_dV2(int k, double s2, double v2_v2, double v2[], double vT[]) const 
 { 
   return (1.0/(36*M_PI))*s2*(-6)*v2[k];
 }

 virtual double dPhi3_dT(int j,int k,double s2, double v2[], double T0[3][3], double TT[3][3]) const 
 { 
   return 0;
 }

 virtual string Name() const { return string("Rosenfeld");}

};


/**
  *  @brief  WhiteBear mark II
  *
  *   Modified WhiteBear which gives better vacancy densities for solid phase.
  */

class WhiteBearII : public WhiteBearI
{
 public:
  WhiteBearII() : WhiteBearI(){};

 virtual bool needsTensor() const { return true;}
  
  virtual double BulkMuex(const vector<double> &x, const vector<Species*> &allSpecies, int species) const
  {
    double n0 = 0.0;
    double n1 = 0.0;
    double n2 = 0.0;
    double n3 = 0.0;
    for(int i=0;i<x.size();i++)
      {
	double density = x[i];
	double hsd = allSpecies[i]->getHSD();
	n0 += density;
	n1 += 0.5*hsd*density;
	n2 += M_PI*hsd*hsd*density;
	n3 += (M_PI/6)*hsd*hsd*hsd*density;
      }


    double dPhi_dn0 = -log(1.0-n3);
    double dPhi_dn1 = (2.0/3.0)*(n2/n3)*log(1-n3)+((5.0-n3)/3)*n2/(1-n3);
    double dPhi_dn2 = ((2.0/3.0)*(n1/n3)-(1.0/(12*M_PI))*(n2*n2/(n3*n3)))*log(1-n3)+((5.0-n3)/3)*(n1/(1-n3))+n2*n2*(-n3*n3+3*n3-1)/(12*M_PI*n3*(1-n3)*(1-n3));
    double dPhi_dn3 = 0.0;

    dPhi_dn3  = (-(2.0/3.0)*(n1*n2/(n3*n3))+(1.0/(18*M_PI))*(n2*n2*n2/(n3*n3*n3)))*log(1-n3);
    dPhi_dn3 += -((2.0/3.0)*(n1*n2/n3)+(n1*n2/3)-(1/(36*M_PI))*(n2*n2*n2/(n3*n3))-n0)*(1.0/(1-n3));
    dPhi_dn3 += ((5.0-n3)/3)*(n1*n2/((1-n3)*(1-n3)));
    dPhi_dn3 +=  (1.0/(36*M_PI))*n2*n2*n2*(1-3*n3+5*n3*n3-n3*n3*n3)/(n3*n3*(1-n3)*(1-n3)*(1-n3));   

    double density = x[species];
    double hsd = allSpecies[species]->getHSD();

    return dPhi_dn0+dPhi_dn1*0.5*hsd+dPhi_dn2*M_PI*hsd*hsd+dPhi_dn3*(M_PI/6)*hsd*hsd*hsd;
      
  }

 virtual double BulkFex(const vector<double> &x, const vector<Species*> &allSpecies) const
 {
   double n0 = 0.0;
   double n1 = 0.0;
   double n2 = 0.0;
   double n3 = 0.0;
   for(int i=0;i<x.size();i++)
     {
       double density = x[i];
       double hsd = allSpecies[i]->getHSD();
       n0 += density;
       n1 += 0.5*hsd*density;
       n2 += M_PI*hsd*hsd*density;
       n3 += (M_PI/6)*hsd*hsd*hsd*density;
     }

   double F = ((2.0/3.0)*(n1*n2/n3)-(1.0/(36*M_PI))*(n2*n2*n2/(n3*n3))-n0)*log(1-n3);
   F += ((5-n3)/3)*n1*n2/(1-n3);
   F += (1.0/(36*M_PI))*n2*n2*n2*(-1+3*n3-n3*n3)/(n3*(1-n3)*(1-n3));
		  
   return F;
 }  

  
  virtual double f2_(double x) const
  {
    if(x < 1e-12)
      return 1+x*(1+x*((10.0/9)+(7.0/6)*x));
    return (1.0/3.0) + (4.0/3.0)*(1.0/(1.0-x)) + (2.0/(3.0*x))*log(1-x);
  }

  virtual double f2p_(double x) const 
  {
    if(x < 1e-12)
      return 1+x*((20.0/9)+3.5*x);
    return (2*(3*x-1)/(3*x*(1-x)*(1-x)))-(2.0/(3.0*x*x))*log(1-x);
  }

  virtual double f3_(double x) const
  {
    if(x < 1e-12)
      return 1.5+(7.0/3)*x+3.25*x*x+4.2*x*x*x;
    return -((1.0-3*x+x*x)/(x*(1-x)*(1-x))) - (1.0/(x*x))*log(1-x);
  }

  virtual double f3p_(double x) const 
  {
    if(x < 1e-12)
      return (7.0/3)+7.5*x+12.6*x*x;
    return ((2-5*x+6*x*x-x*x*x)/(x*x*(1-x)*(1-x)*(1-x)))+(2.0/(x*x*x))*log(1-x);
  }
  virtual string Name() const { return string("WhiteBear II");}
};


/**
  *  @brief  Exception used to isolate eta(r) > 1 errors
  *
  */  
class Eta_Too_Large_Exception : public runtime_error
{
 public:
 Eta_Too_Large_Exception() : runtime_error("Eta too large in FMT::dPHI"){}
  virtual const char* what() const throw() 
  {
    cnvt.str( "" );

    cnvt << runtime_error::what() << ": more info ... ";

    return cnvt.str().c_str();
  }
  private:
  static ostringstream cnvt;
};  



#endif // __LUTSKO_WEIGHTED_DENSITIES__
