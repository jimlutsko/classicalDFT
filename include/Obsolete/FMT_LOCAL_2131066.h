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
#include "Fundamental_Measures.h"


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
  virtual void get_dPdx_coeffs(vector<double> &num, vector<double> &denom) const = 0;
  virtual void get_P_coeffs(vector<double> &num, vector<double> &denom) const = 0;
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


  /**
  *   @brief  Calculates (d2F/dn_i dn_j)v_j
  *  
  *   @param  v: input vector
  *   @param  d2F: vector to be filled
  */    
  void add_second_derivative(vector<DFT_FFT> &v, vector<DFT_Vec> &d2F, vector<Species*> &allSpecies);

  // Brute force calculation: for checking only!
  double d2Phi_dn_dn(int I[3], int si, int J[3], int sj, vector<Species*> &allSpecies);

  
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
  *   @brief  Calculate phi given a set of fundamental measures
  *
  *   @param  fm is a FundamentalMeasures structure
  *   @return phi
  */    
  double Phi(const FundamentalMeasures &fm) const;

  /**
  *   @brief  Calculate dphi/dn_a for all fundamental measures, n_a
  *
  *   @param  fm are the fundamental measures
  *   @param  dPhi is the result
  */      
  void DPhi(const FundamentalMeasures& fm, FundamentalMeasures& dPhi) const;
  void D2Phi(const FundamentalMeasures& fm, vector<vector<double>>& d2Phi) const;
  
  /**
  *   @brief  Calculates the fundamental measures at lattice point i
  *
  *   @param  i is the lattice point
  *   @param allSpecies is the array of FMT_Species
  *   @return phi
  */      
  FundamentalMeasures getWeightedDensities(long i, vector<Species*> &allSpecies);

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
  *   @brief  secondderivative of f1_() with respect to its argument.
  *  
  *   @param  eta is the local volume-average density
  *   @return f1(eta).
  */     
 virtual double f1pp_(double eta) const
 {
    if(eta > etaMax_)
      {
	double a1 = -1.0/(1.0-etaMax_);
	double a2 = -a1*a1;
	return a2;
      }
    
    return -1.0/(1.0-2*eta+eta*eta);
 }
 
 /**
  *   @brief  derivative of f2_() with respect to its argument.
  *  
  *   @param  eta is the local volume-average density
  *   @return f2(eta).
  */     
  virtual double f2p_(double eta) const = 0;
  virtual double f2pp_(double eta) const = 0;

 /**
  *   @brief  derivative of f3_() with respect to its argument.
  *  
  *   @param  eta is the local volume-average density
  *   @return f3(eta).
  */     
  virtual double f3p_(double eta) const= 0;
  virtual double f3pp_(double eta) const= 0;

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
 virtual double Phi3(const FundamentalMeasures &fm) const = 0;

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
  virtual double dPhi3_dS2(const FundamentalMeasures &fm) const = 0;
  virtual double dPhi3_dS2_dS2(const FundamentalMeasures &fm) const { throw std::runtime_error("dPhi3_dS2_dS2 not implemented in FMT-derived class");;}
  virtual double dPhi3_dV2_dS2(int i, const FundamentalMeasures &fm) const { throw std::runtime_error("dPhi3_dV2_dS2 not implemented in FMT-derived class");;}
  virtual double dPhi3_dV2_dV2(int i, int j, const FundamentalMeasures &fm) const { throw std::runtime_error("dPhi3_dV2_dS2 not implemented in FMT-derived class");;}

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
 virtual double dPhi3_dV2(int k, const FundamentalMeasures &fm) const = 0;

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
  virtual double dPhi3_dT(int j,int k, const FundamentalMeasures &fm) const = 0;


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


  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & etaMax_;
  }
  
 protected:
 
 double etaMax_; ///< cutoff used to control divergences
 
};

/**
  *  @brief  The original White Bear  FMT model 
  *
  *   White-bear FMT model with tensor densities and the CS equation of state
  */  
class Stable : public FMT
{
 public:
  Stable(double A = 1, double B = 0) : FMT(), A_(A), B_(B){};

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

  virtual double f2pp_(double eta) const
  {
    double f = 1.0/(1.0-eta);
    return 2*f*f*f;
  }  

  virtual double f3_(double eta) const
  {
    double f = 1.0/(1.0-eta);
    return f*f;
  }

  virtual double f3p_(double eta) const
  {
    double f = 1.0/(1.0-eta);
    return 2*f*f*f;
  }
  virtual double f3pp_(double eta) const
  {
    double f = 1.0/(1.0-eta);
    return 6*f*f*f*f;    
  }  

  virtual double Phi3(const FundamentalMeasures &fm) const
  {
    double s2     = fm.s2;
    double vTv    = fm.vTv;
    double v2_v2  = fm.v2_v2;
    double T2     = fm.T2;
    double T3     = fm.T3;
    
    return (A_/(24*M_PI))*(s2*s2*s2-3*s2*v2_v2+3*vTv-T3)
      +(B_/(24*M_PI))*(s2*s2*s2-3*s2*T2+2*T3);
  }

 virtual double dPhi3_dS2(const FundamentalMeasures &fm) const 
 {
   double s2    = fm.s2;   
   double v2_v2 = fm.v2_v2;
   double T2    = fm.T2;
     
   return (A_/(24*M_PI))*(3*s2*s2-3*v2_v2) + (B_/(24*M_PI))*(3*s2*s2-3*T2);
 }

  virtual double dPhi3_dS2_dS2(const FundamentalMeasures &fm) const
  {
    double s2    = fm.s2;   
   
    return (A_/(24*M_PI))*(6*s2) + (B_/(24*M_PI))*(6*s2);
  }
  virtual double dPhi3_dV2_dS2(int i, const FundamentalMeasures &fm) const
  {
   return (A_/(24*M_PI))*(-6*fm.v2[i]);
  }
  virtual double dPhi3_dV2_dV2(int i, int j, const FundamentalMeasures &fm) const
  {
    return (A_/(24*M_PI))*(6*fm.T[i][j]-6*(i == j ? fm.s2 : 0.0));
  }
  
  virtual double dPhi3_dV2(int k, const FundamentalMeasures &fm) const 
 {
   double s2    = fm.s2;
   double v2_v2 = fm.v2_v2;
   double v2_k   = fm.v2[k];
   double vT_k   = fm.vT[k];
   
   return (A_/(24*M_PI))*(-6*s2*v2_k+6*vT_k);
 }

 virtual double dPhi3_dT(int j,int k, const FundamentalMeasures &fm) const 
 {
   double s2     = fm.s2;
   double v2_j   = fm.v2[j];
   double v2_k   = fm.v2[k];
   double T_T_jk = fm.TT[j][k];
   double T_jk   = fm.T[j][k];
   
    return (A_/(24*M_PI))*(3*v2_j*v2_k-3*T_T_jk)  + (B_/(24*M_PI))*(-6*s2*T_jk+6*T_T_jk);
 }
  //stable
  // phi = -n0*log(1-n3) +(n1*n2/(1-n3))+(1/(24*M_PI))*(8/0.9)*n2*n2*n2/(1-n3)^2
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
    double dPhi_dn2 = (n1/(1-n3))+(1.0/(27*M_PI))*3*n2*n2/((1-n3)*(1-n3));
    double dPhi_dn3 = 0.0;

    dPhi_dn3  = n0/(1-n3);
    dPhi_dn3 += (n1*n2/((1-n3)*(1-n3)));
    dPhi_dn3 += 2*(1.0/(27*M_PI))*(A_+B_/4)*n2*n2*n2/((1-n3)*(1-n3)*(1-n3));

    double density = x[species];
    double hsd = allSpecies[species]->getHSD();

    return dPhi_dn0+dPhi_dn1*0.5*hsd+dPhi_dn2*M_PI*hsd*hsd+dPhi_dn3*(M_PI/6)*hsd*hsd*hsd; 
  }

  virtual void get_P_coeffs(vector<double> &num, vector<double> &denom) const
  {
    num.clear();
    denom.clear();

    double C = (8*A_+2*B_-6)/3;
    
    num.push_back(1);
    num.push_back(1);
    num.push_back(C);

    denom.push_back( 1);
    denom.push_back(-3);
    denom.push_back( 3);
    denom.push_back(-1);
  }

  virtual void get_dPdx_coeffs(vector<double> &num, vector<double> &denom) const
  {
    num.clear();
    denom.clear();

    double C = (8*A_+2*B_-6)/3;    
    
    num.push_back( 1);
    num.push_back( 4);
    num.push_back( 1+3*C);

    denom.push_back( 1);
    denom.push_back(-4);
    denom.push_back( 6);
    denom.push_back(-4);
    denom.push_back( 1);
    
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

   double F = 0;
   F += -n0*log(1-n3);
   F += n1*n2/(1-n3);
   F += (1.0/(27*M_PI))*(A_+B_/4)*n2*n2*n2/((1-n3)*(1-n3));
   return F;
 }

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<FMT>(*this);
    boost::serialization::void_cast_register<Stable, FMT>(static_cast<Stable *>(NULL),static_cast<FMT *>(NULL));    
  }

  
  virtual string Name() const { return string("Stable with A = ") + to_string(A_) + string(" and B = ") + to_string(B_);}

protected:
  double A_ = 1;
  double B_ = 0;
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

  virtual double f2pp_(double eta) const
  {
    double f = 1.0/(1.0-eta);
    return 2*f*f*f;
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
    return ((-(1-3*eta)/((1-eta)*(1-eta)*(1-eta)))-(1.0/(1-eta))-(2.0/eta)*log(1-eta))/(eta*eta);
  }
  virtual double f3pp_(double eta) const { return 0.0;}  

  virtual double Phi3(const FundamentalMeasures &fm) const
  {
    double s2     = fm.s2;
    double vTv    = fm.vTv;
    double v2_v2  = fm.v2_v2;
    double T2     = fm.T2;
    double T3     = fm.T3;
    
    return (1.0/(8*M_PI))*(vTv-s2*v2_v2-T3+s2*T2);
  }

 virtual double dPhi3_dS2(const FundamentalMeasures &fm) const 
 {
   double v2_v2 = fm.v2_v2;
   double T2    = fm.T2;
   
   return (1.0/(8*M_PI))*(-v2_v2+T2);
 }

  virtual double dPhi3_dS2_dS2(const FundamentalMeasures &fm) const {return 0.0;}
  virtual double dPhi3_dV2_dS2(int i, const FundamentalMeasures &fm) const { return (1.0/(8*M_PI))*(-2*fm.v2[i]);}
  virtual double dPhi3_dV2_dV2(int i, int j, const FundamentalMeasures &fm) const
  {
    return (1.0/(8*M_PI))*2*(fm.T[i][j]-(i == j ? fm.s2 : 0.0));
  }
  
  virtual double dPhi3_dV2(int k, const FundamentalMeasures &fm) const 
 {
   double s2    = fm.s2;
   double v2_v2 = fm.v2_v2;
   double v2_k   = fm.v2[k];
   double vT_k   = fm.vT[k];
   
   return (1.0/(8*M_PI))*(2*vT_k-2*s2*v2_k);
 }

 virtual double dPhi3_dT(int j,int k, const FundamentalMeasures &fm) const 
 {
   double s2     = fm.s2;
   double v2_j   = fm.v2[j];
   double v2_k   = fm.v2[k];
   double T_T_jk = fm.TT[j][k];
   double T_jk   = fm.T[j][k];
   
   return (1.0/(8*M_PI))*(v2_j*v2_k-3*T_T_jk+2*s2*T_jk);
 }


  virtual void get_P_coeffs(vector<double> &num, vector<double> &denom) const
  {
    num.clear();
    denom.clear();

    num.push_back(1);
    num.push_back(1);
    num.push_back(1);
    num.push_back(-1);

    denom.push_back( 1);
    denom.push_back(-3);
    denom.push_back( 3);
    denom.push_back(-1);
  }

  virtual void get_dPdx_coeffs(vector<double> &num, vector<double> &denom) const
  {
    num.clear();
    denom.clear();
    
    num.push_back( 1);
    num.push_back( 4);
    num.push_back( 4);
    num.push_back(-4);
    num.push_back( 1);

    denom.push_back( 1);
    denom.push_back(-4);
    denom.push_back( 6);
    denom.push_back(-4);
    denom.push_back( 1);
    
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

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<FMT>(*this);
    boost::serialization::void_cast_register<WhiteBearI, FMT>(static_cast<WhiteBearI *>(NULL),static_cast<FMT *>(NULL));    
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

  virtual double f2pp_(double eta) const
  {
    if(eta > etaMax_)
      {
	double a0 = 1.0/(1.0-etaMax_);	
	double a1 = a0*a0;
	double a2 = 2*a1*a0;
	return a2;
      }
    
    double f = 1.0/(1.0-eta);
    return 2*f*f*f;
  }

  virtual void get_P_coeffs(vector<double> &num, vector<double> &denom) const
  {
    num.clear();
    denom.clear();

    num.push_back(1);
    num.push_back(1);
    num.push_back(1);
    num.push_back(-1);

    denom.push_back( 1);
    denom.push_back(-3);
    denom.push_back( 3);
    denom.push_back(-1);
  }

  virtual void get_dPdx_coeffs(vector<double> &num, vector<double> &denom) const
  {
    num.clear();
    denom.clear();
    
    num.push_back( 1);
    num.push_back( 4);
    num.push_back( 4);
    num.push_back(-4);
    num.push_back( 1);

    denom.push_back( 1);
    denom.push_back(-4);
    denom.push_back( 6);
    denom.push_back(-4);
    denom.push_back( 1);
    
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


  virtual double f3pp_(double eta) const
  {
    if(eta < 1e-12)
      return 7.5+28.8*eta;

    if(eta > etaMax_)
      {
	double x = etaMax_;
	double a1 = -(1.0/(x*x))*(1.0/((1.0-x)*(1.0-x)*(1.0-x)))*(x*x-5*x+2)-(2.0/(x*x*x))*log(1-x);
	double a2 = -(1.0/(x*x*x))*(1.0/((1.0-x)*(1.0-x)*(1.0-x)*(1.0-x)))*(5*x*x*x-26*x*x+21*x-6)
	  +(6.0/(x*x*x*x))*log(1-x);	
	return a2;
      }
    
    return -(1.0/(eta*eta*eta))*(1.0/((1.0-eta)*(1.0-eta)*(1.0-eta)*(1-eta)))*(5*eta*eta*eta-26*eta*eta+21*eta-6)+(6.0/(eta*eta*eta*eta))*log(1-eta);
  }

  
  virtual double Phi3(const FundamentalMeasures &fm) const
  {
    double s2     = fm.s2;
    double v2_v2  = fm.v2_v2;
    
    if(s2 < 1e-20) return 0.0;
    double psi = v2_v2/(s2*s2);
    return (1.0/(36*M_PI))*s2*s2*s2*(1-psi)*(1-psi)*(1-psi);
  }

  virtual double dPhi3_dS2(const FundamentalMeasures &fm) const 
  {
    double s2    = fm.s2;
    double v2_v2 = fm.v2_v2;
    
    if(s2 < 1e-20) return 0.0;
    double psi = v2_v2/(s2*s2);
    return (1.0/(12*M_PI))*s2*s2*(1+psi)*(1-psi)*(1-psi);
  }

  virtual double dPhi3_dS2_dS2(const FundamentalMeasures &fm) const 
  {
    double s2    = fm.s2;
    double v2_v2 = fm.v2_v2;
    
    if(s2 < 1e-20) return 0.0;
    double psi = v2_v2/(s2*s2);
    return (1.0/(6*M_PI))*s2*(1-psi)*(1+psi+2*psi*psi);
  }  

  virtual double dPhi3_dV2_dS2(int k, const FundamentalMeasures &fm) const
  {
    double s2    = fm.s2;
    double v2_v2 = fm.v2_v2;
    double v2_k  = fm.v2[k];
    
    if(s2 < 1e-20) return 0.0;
    double psi = v2_v2/(s2*s2);
    return (1.0/(36*M_PI))*(-6*v2_k)*(1-psi)*(1+3*psi);
  }
  virtual double dPhi3_dV2_dV2(int i, int j, const FundamentalMeasures &fm) const
  {
    double s2    = fm.s2;
    double v2_v2 = fm.v2_v2;
    double v2_i  = fm.v2[i];
    double v2_j  = fm.v2[j];

    if(s2 < 1e-20) return 0.0;
    double psi = v2_v2/(s2*s2);
    
    if(i == j) return (1.0/(36*M_PI))*(-6*s2)*(1-psi)*(1-psi-4*v2_i*v2_i/(s2*s2));

    return (1.0/(36*M_PI))*(24/s2)*v2_i*v2_j*(1-psi);
  }  

  
  virtual double dPhi3_dV2(int k, const FundamentalMeasures &fm) const
  {
    double s2    = fm.s2;
    double v2_v2 = fm.v2_v2;
    double v2_k  = fm.v2[k];
    
    if(s2 < 1e-20) return 0.0;
    double psi = v2_v2/(s2*s2);
    return -(1.0/(6*M_PI))*s2*v2_k*(1-psi)*(1-psi);
  }

  virtual double dPhi3_dT(int j,int k, const FundamentalMeasures &fm) const 
  { 
    return 0;
  }

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<FMT>(*this);
    boost::serialization::void_cast_register<RSLT, FMT>(static_cast<RSLT *>(NULL),static_cast<FMT *>(NULL));    
  }

  
  virtual string Name() const { return string("RSLT");}
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

  virtual void get_d2Phi(vector< vector<double> > &d2Phi, double eta, double s0, double s1, double s2, double v1[3], double v2[3])
  { throw std::runtime_error("get_d2Phi not implemented in clas Rosenfeld");}  

 virtual bool needsTensor() const { return false;}  
  
  virtual double Phi3(const FundamentalMeasures &fm) const
  {
    double s2     = fm.s2;
    double v2_v2  = fm.v2_v2;
    return (1.0/(36*M_PI))*(s2*s2*s2-3*s2*v2_v2);
  }

 virtual double dPhi3_dS2(const FundamentalMeasures &fm) const 
 {
   double s2    = fm.s2;
   double v2_v2 = fm.v2_v2;
   
   return (1.0/(36*M_PI))*(3*s2*s2-3*v2_v2);
 }

 virtual double dPhi3_dS2_dS2(const FundamentalMeasures &fm) const 
 {
   double s2    = fm.s2;
   double v2_v2 = fm.v2_v2;
   
   return (1.0/(6*M_PI))*s2;
 }  

 virtual double dPhi3_dV2(int k, const FundamentalMeasures &fm) const 
 {
   double s2  = fm.s2;
   double v2_k = fm.v2[k];
   
   return (1.0/(36*M_PI))*s2*(-6)*v2_k;
 }

  virtual double dPhi3_dV2_dS2(int k, const FundamentalMeasures &fm) const
  {
    double v2_k  = fm.v2[k];

    return -v2_k/(6*M_PI);
  }

  virtual double dPhi3_dV2_dV2(int i, int j, const FundamentalMeasures &fm) const
  {
    double s2    = fm.s2;
    return (i == j ? -s2/(6*M_PI) : 0.0);
  }  
  
 virtual double dPhi3_dT(int j,int k, const FundamentalMeasures &fm) const 
 { 
   return 0;
 }

 virtual string Name() const { return string("Rosenfeld");}

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<RSLT>(*this);
    boost::serialization::void_cast_register<Rosenfeld, RSLT>(static_cast<Rosenfeld *>(NULL),static_cast<RSLT *>(NULL));    
  }

  
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
  virtual double f3pp_(double eta) const { return 0.0;}
  
  virtual string Name() const { return string("WhiteBear II");}

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<WhiteBearI>(*this);
    boost::serialization::void_cast_register<WhiteBearII, WhiteBearI>(static_cast<WhiteBearII *>(NULL),static_cast<WhiteBearI *>(NULL));    
  }
  
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
