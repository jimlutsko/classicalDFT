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
  FMT(){} 
  /**
   *   @brief  Default  destrctur for FMT 
   *  
   *   @return nothing 
   */  
  ~FMT(){}

  /**
   *   @brief  This function calculates total FMT contribution to the chemical potential of the requested species.
   *           It allows for the complications of extended models (like AO model). 

   *   @param  density is an array of densities (one for each species).
   *   @param  allSpecies is the array of species (and obviously matching the density array).
   *   @param  species is the species for which the chemical potential is returned. 
   *   @return Total FMT contribution to excess beta times Helmholtz free energy per unit volume. 
   */   
  double BulkMuex(const vector<double> &x, const vector<Species*> &allSpecies, int species) const
  {
    FundamentalMeasures fm(0.0,1.0);

    for(int s=0; s<allSpecies.size(); s++)
      {
	FMT_Species *sfmt = dynamic_cast<FMT_Species*>(allSpecies[s]);
	if(sfmt)
	  {
	    FundamentalMeasures fm1(x[s],sfmt->getHSD());
	    fm.add(fm1);
	  }
      }
    FundamentalMeasures dPhi;
    calculate_dPhi_wrt_fundamental_measures(fm, dPhi);
    vector<double> dPhi_v = dPhi.get_as_vector();
    
    FundamentalMeasures weights(1.0, allSpecies[species]->getHSD()); // a non-fmt species gives zero hsd so OK.
    vector<double> weights_v = weights.get_as_vector();

    double mu = 0.0;
    for(int i=0;i<dPhi_v.size();i++) mu += dPhi_v[i]*weights_v[i];

    // I now assume there is only one species as I do not (now, yet) have the generalization
    // of these expressions to the case of multiple species.
    FMT_AO_Species *sao = dynamic_cast<FMT_AO_Species*>(allSpecies[0]);
    if(sao)
      {
	double hsdp = sao->getHSDP();
	double rhop = sao->getReservoirDensity();
	FundamentalMeasures n(x[0],sao->getHSD());
	FundamentalMeasures dPhi;  
	calculate_dPhi_wrt_fundamental_measures(n,dPhi);
	vector<double> v_dPhi = dPhi.get_as_vector();
	
	vector<double> v_wp = (FundamentalMeasures(1.0,hsdp)).get_as_vector();

	FundamentalMeasures w(1.0,sao->getHSD());  
	FundamentalMeasures d2Phi_V;	    
	calculate_d2Phi_dot_V(n,w,d2Phi_V);
	vector<double> result = d2Phi_V.get_as_vector();
	
	double arg  = 0;
	double pref = 0;
	for(int i=0;i<v_dPhi.size();i++)
	  {
	    arg  += v_dPhi[i]*v_wp[i];
	    pref += result[i]*v_wp[i];
	  }
	    
	mu += rhop*pref*exp(-arg);
      }    

    return mu;
  }

  /**
   *   @brief  This function calculates total FMT contribution to the helmholtz free energy.
   *           It allows for the complications of extended models (like AO model). 

   *   @param  density is an array of densities (one for each species).
   *   @param  allSpecies is the array of species (and obviously matching the density array).
   *   @return Total FMT contribution to excess beta times Helmholtz free energy per unit volume. 
   */      
  double BulkFex(const vector<double> &x, const vector<Species*> &allSpecies) const
  {
    FundamentalMeasures fm(0.0,1.0);

    for(int s=0; s<allSpecies.size(); s++)
      {
	FMT_Species *sfmt = dynamic_cast<FMT_Species*>(allSpecies[s]);
	if(sfmt)
	  {
	    FundamentalMeasures fm1(x[s],sfmt->getHSD());
	    fm.add(fm1);
	  }
      }
    double f = calculate_Phi(fm);

    // I now assume there is only one species as I do not (now, yet) have the generalization
    // of these expressions to the case of multiple species.
    FMT_AO_Species *sao = dynamic_cast<FMT_AO_Species*>(allSpecies[0]);
    if(sao)
      {
	double hsdp = sao->getHSDP();
	double rhop = sao->getReservoirDensity();
	FundamentalMeasures n(x[0],sao->getHSD());
	FundamentalMeasures dPhi;  
	calculate_dPhi_wrt_fundamental_measures(n,dPhi);
	vector<double> v_dPhi = dPhi.get_as_vector();
       
	vector<double> v_wp = (FundamentalMeasures(1.0,hsdp)).get_as_vector();
	
	double arg = 0;
	for(int i=0;i<v_dPhi.size();i++) arg += v_dPhi[i]*v_wp[i];
	f -= rhop*exp(-arg);

	// This is the "standard" shift of the free energy density. 
	f += rhop;	
      }
    return f;
  }
  
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
  double calculate_Phi(const FundamentalMeasures &fm) const;

  /**
   *   @brief  Calculate dphi/dn_a for all fundamental measures, n_a
   *
   *   @param  fm are the fundamental measures
   *   @param  dPhi is the result
   */      
  void calculate_dPhi_wrt_fundamental_measures(const FundamentalMeasures& fm, FundamentalMeasures& dPhi) const;

  void calculate_d2Phi_dot_V(const FundamentalMeasures& n, const FundamentalMeasures &v, FundamentalMeasures &result) const;  
  
  /**
   *   @brief  Calculates the fundamental measures at lattice point i
   *
   *   @param  i is the lattice point
   *   @param allSpecies is the array of FMT_Species
   *   @return phi
   */      
  FundamentalMeasures getWeightedDensities(long i, vector<Species*> &allSpecies);

  // The FMT functional is coded as PHI1 = f1(eta)*Phi1(eta,s,v,T) etc.
  // Phi1 and Phi2 are the same in all models (so far) so they are directly coded in the implementation, FMT.cpp. Phi3 must be specified.
  virtual double f1_(double eta)  const = 0;
  virtual double f2_(double eta)  const = 0;
  virtual double f3_(double eta)  const= 0;

  //derivatives
  virtual double f1p_(double eta) const = 0;
  virtual double f1pp_(double eta) const = 0;

  virtual double f2p_(double eta) const = 0;
  virtual double f2pp_(double eta) const = 0;

  virtual double f3p_(double eta) const= 0;
  virtual double f3pp_(double eta) const= 0;

  /**
   *   @brief  calculates third contribution to PHI = PHI1 + PHI2 + PHI3. This is model dependent
   */     
  virtual double Phi3(const FundamentalMeasures &fm) const = 0;

  /**
   *   @brief  calculates derivative of third contribution to PHI = PHI1 + PHI2 + PHI3 with respect to S2
   */     
  virtual double dPhi3_dS2(const FundamentalMeasures &fm) const = 0;
  virtual double dPhi3_dS2_dS2(const FundamentalMeasures &fm) const { throw std::runtime_error("dPhi3_dS2_dS2 not implemented in FMT-derived class");;}
  virtual double dPhi3_dV2_dS2(int i, const FundamentalMeasures &fm) const { throw std::runtime_error("dPhi3_dV2_dS2 not implemented in FMT-derived class");;}
  virtual double dPhi3_dV2_dV2(int i, int j, const FundamentalMeasures &fm) const { throw std::runtime_error("dPhi3_dV2_dS2 not implemented in FMT-derived class");;}

  /**
   *   @brief  calculates derivative of third contribution to PHI = PHI1 + PHI2 + PHI3 with respect to v2(k)
   */      
  virtual double dPhi3_dV2(int k, const FundamentalMeasures &fm) const = 0;

  /**
   *   @brief  calculates derivative of third contribution to PHI = PHI1 + PHI2 + PHI3 with respect to T(j,k)
   */       
  virtual double dPhi3_dT(int j,int k, const FundamentalMeasures &fm) const  { return 0;}
  virtual double dPhi3_dS2_dT(int j,int k, const FundamentalMeasures &fm) const { return 0;}
  virtual double dPhi3_dV2_dT(int i, int j,int k, const FundamentalMeasures &fm) const { return 0;}
  virtual double dPhi3_dT_dT(int i, int j,int k, int l, const FundamentalMeasures &fm) const{ return 0;}

  virtual bool needsTensor() const  = 0;
  
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
  }
  
};

/**
 *  @brief  Original Rosenfeld 
 *
 */  

class Rosenfeld: public FMT
{
 public:
  Rosenfeld() : FMT(){};

  virtual bool needsTensor() const { return false;}

  virtual double f1_(double eta)  const
  {
    return log(1.0-eta);
  }
  
 virtual double f1p_(double eta) const
 {
    return -1.0/(1.0-eta);
 }

 virtual double f1pp_(double eta) const
 {
   return -1.0/(1.0-2*eta+eta*eta);
 }

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
    double v2_v2  = fm.v2_v2;

    return (1.0/(24*M_PI))*(s2*s2*s2-3*s2*v2_v2);
  }

 virtual double dPhi3_dS2(const FundamentalMeasures &fm) const 
 {
   double s2    = fm.s2;   
   double v2_v2 = fm.v2_v2;
     
   return (1.0/(24*M_PI))*(3*s2*s2-3*v2_v2);
 }

  virtual double dPhi3_dS2_dS2(const FundamentalMeasures &fm) const
  {
    double s2    = fm.s2;   
   
    return (1.0/(24*M_PI))*(6*s2);
  }
  virtual double dPhi3_dV2_dS2(int i, const FundamentalMeasures &fm) const
  {
    return (1.0/(24*M_PI))*(-6*fm.v2[i]);
  }
  virtual double dPhi3_dV2_dV2(int i, int j, const FundamentalMeasures &fm) const
  {
    return 0.0;
  }
  
  virtual double dPhi3_dV2(int k, const FundamentalMeasures &fm) const 
  {
    double s2    = fm.s2;
    double v2_v2 = fm.v2_v2;
    double v2_k   = fm.v2[k];
   
    return (1.0/(24*M_PI))*(-6*s2*v2_k);
  }

  virtual void get_d2Phi(vector< vector<double> > &d2Phi, double eta, double s0, double s1, double s2, double v1[3], double v2[3])
  { throw std::runtime_error("get_d2Phi not implemented in clas Rosenfeld");}  


  virtual string Name() const { return string("Rosenfeld");}

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<FMT>(*this);
    boost::serialization::void_cast_register<Rosenfeld, FMT>(static_cast<Rosenfeld *>(NULL),static_cast<FMT *>(NULL));    
  }  
};

/**
  *  @brief  The RSLT positive-definite FMT model (pre-White Bear)
  *
  *   This class implements the RSLT FMT model which is not as accurate as WhiteBear but which is positive definite and therefore stable.
  *   Differences from Rosenfeld are all in the Phi3 term. Note that (unlike original RSLT) this gives CS equation of state.
  */  
class RSLT : public Rosenfeld
{
 public:
  RSLT() : Rosenfeld(){};

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

    return -(1.0/(eta*eta))*(1.0/((1.0-eta)*(1.0-eta)*(1.0-eta)))*(eta*eta-5*eta+2)-(2.0/(eta*eta*eta))*log(1-eta);
  }


  virtual double f3pp_(double eta) const
  {
    if(eta < 1e-12)
      return 7.5+28.8*eta;

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

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Rosenfeld>(*this);
    boost::serialization::void_cast_register<RSLT, Rosenfeld>(static_cast<RSLT *>(NULL),static_cast<Rosenfeld *>(NULL));    
  }
  
  virtual string Name() const { return string("RSLT");}
};


/**
  *  @brief  "Explicitly Stable" FMT
  *
  *   Note that the bulk eos for the fluid is neither PY nor CS.
  */  
class esFMT : public Rosenfeld
{
 public:
  esFMT(double A = 1, double B = 0) : Rosenfeld(), A_(A), B_(B){};

  virtual bool needsTensor() const { return true;}

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
    return (A_/(24*M_PI))*(3*fm.T[i][j]+3*fm.T[j][i]-6*(i == j ? fm.s2 : 0.0));
  }
  
  virtual double dPhi3_dV2(int k, const FundamentalMeasures &fm) const 
  {
    double s2    = fm.s2;
    double v2_v2 = fm.v2_v2;
    double v2_k   = fm.v2[k];
    double vT_k   = fm.vT[k];
    double Tv_k   = fm.Tv[k];
   
    return (A_/(24*M_PI))*(-6*s2*v2_k+3*vT_k+3*Tv_k);
  }

  virtual double dPhi3_dT(int j,int k, const FundamentalMeasures &fm) const 
  {
   double s2     = fm.s2;
   double v2_j   = fm.v2[j];
   double v2_k   = fm.v2[k];
   double T_T_jk = fm.TT[k][j];
   double T_jk   = fm.T[k][j];
   
   return (A_/(8*M_PI))*(v2_j*v2_k-T_T_jk)  + (B_/(4*M_PI))*(-s2*T_jk+T_T_jk);
  }
  virtual double dPhi3_dS2_dT(int j,int k, const FundamentalMeasures &fm) const
  {
    double s2     = fm.s2;
    double v2_j   = fm.v2[j];
    double v2_k   = fm.v2[k];
    double T_T_jk = fm.TT[k][j];
    double T_jk   = fm.T[k][j];
    
    return (B_/(4*M_PI))*(-T_jk);
  }    
  virtual double dPhi3_dV2_dT(int i, int j,int k, const FundamentalMeasures &fm) const
  {
    double s2     = fm.s2;
    double v2_j   = fm.v2[j];
    double v2_k   = fm.v2[k];
    double T_T_jk = fm.TT[j][k];
    double T_jk   = fm.T[j][k];

    double val = 0;
    if(i == j) val += A_*v2_k;
    if(i == k) val += A_*v2_j;
    
    return val/(8*M_PI);
  }        
  virtual double dPhi3_dT_dT(int j, int k,int l, int m, const FundamentalMeasures &fm) const
  {
    double s2     = fm.s2;

    double val = 0;
    if(l == k) val += -A_*fm.T[m][j] + 2*B_*fm.T[m][j];
    if(m == j) val += -A_*fm.T[k][l] + 2*B_*fm.T[k][l];
    if(l == k && m == j) val += -2*B_*s2;
    return val/(8*M_PI);
  }        
  

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Rosenfeld>(*this);
    boost::serialization::void_cast_register<esFMT, Rosenfeld>(static_cast<esFMT *>(NULL),static_cast<Rosenfeld *>(NULL));    
  }

  
  virtual string Name() const { return string("esFMT with A = ") + to_string(A_) + string(" and B = ") + to_string(B_);}

protected:
  double A_ = 1;
  double B_ = 0;
};



/**
  *  @brief  Variant of the "Explicitly Stable" FMT to enforce high vacancies
  *
  *   A parameter "lambda" is added to the numerator of Phi3 so that there
  *   is a repulsion term "lambda/(1-eta)^2" that prevent the minimum from
  *   being too close to the eta>=1 region which is unphysical and causes 
  *   backtracking in calculations. This is an artificial functional that
  *   should be used at the early stage of the minimisation.
  */  
class esFMT_lambda : public esFMT
{
 public:
  esFMT_lambda(double A = 1, double B = 0, double lambda = 0) 
  : esFMT(A,B), lambda_(lambda){};

  virtual bool needsTensor() const { return true;}

  virtual double Phi3(const FundamentalMeasures &fm) const
  {
    double s2     = fm.s2;
    double vTv    = fm.vTv;
    double v2_v2  = fm.v2_v2;
    double T2     = fm.T2;
    double T3     = fm.T3;
    
    return (A_/(24*M_PI))*(s2*s2*s2-3*s2*v2_v2+3*vTv-T3)
      +(B_/(24*M_PI))*(s2*s2*s2-3*s2*T2+2*T3) + lambda_;
  }
  
 protected:
  double lambda_;
};



/**
  *  @brief  The original White Bear  FMT model 
  *
  *   White-bear FMT model with tensor densities and the CS equation of state
  */  
class WhiteBearI : public esFMT
{
 public:

  // In the paper, it says that this should be esFMT(3/2, -3/2): here the 3/2 has been moved into f3. 
  WhiteBearI() : esFMT(1,-1){};

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

  virtual double f3_(double x) const
  {
    if(x < 1e-12)
      return 1.5+(8.0/3)*x+3.75*x*x+4.8*x*x*x;
    
    return (1.0/(x*(1-x)*(1-x))) +  (log(1-x)/(x*x));
  }

  virtual double f3p_(double x) const
  {
    if(x < 1e-12)
      return (8.0/3)+7.5*x+14.4*x*x;
    
    return ((-2+5*x-x*x)/(x*x*(1-x)*(1-x)*(1-x)))-(2/(x*x*x))*log(1-x);
  }
  virtual double f3pp_(double x) const
  {
    if(x < 1e-12)
      return 7.5+28.8*x;
    
    return ((6-21*x+26*x*x-5*x*x*x)/(x*x*x*(1-x)*(1-x)*(1-x)*(1-x)))+(6*log(1-x))/(x*x*x*x);
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
  *  @brief  WhiteBear mark II
  *
  *   Modified WhiteBear which gives better vacancy densities for solid phase.
  */

class WhiteBearII : public esFMT //WhiteBearI
{
 public:

  // In the paper, it says that this should be esFMT(3/2, -3/2): here the 3/2 has been moved into f3. 
  WhiteBearII() : esFMT(1,-1){};
  
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

  virtual double f2pp_(double x) const 
  {
    if(x < 1e-12)
      return (20.0/9)+7*x;
    return (2*(7*x*x-5*x+2)/(3*x*x*(1-x)*(1-x)*(1-x)))+(4.0/(3.0*x*x*x))*log(1-x);
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
  virtual double f3pp_(double x) const
  {
    if(x < 1e-12)
      return 7.5+25.2*x;
    
    return (-(6-21*x+26*x*x-19*x*x*x+2*x*x*x*x)/(x*x*x*(1-x)*(1-x)*(1-x)*(1-x)))-(6.0/(x*x*x*x))*log(1-x);    
  }
  
  virtual string Name() const { return string("WhiteBear II");}

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<WhiteBearI>(*this);
    boost::serialization::void_cast_register<WhiteBearII, WhiteBearI>(static_cast<WhiteBearII *>(NULL),static_cast<WhiteBearI *>(NULL));    
  }
  
};



// lambda variant

class WhiteBearII_lambda : public esFMT_lambda
{
 public:
  WhiteBearII_lambda(double lambda) : esFMT_lambda(1,-1, lambda){};
  
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

  virtual double f2pp_(double x) const 
  {
    if(x < 1e-12)
      return (20.0/9)+7*x;
    return (2*(7*x*x-5*x+2)/(3*x*x*(1-x)*(1-x)*(1-x)))+(4.0/(3.0*x*x*x))*log(1-x);
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
  virtual double f3pp_(double x) const
  {
    if(x < 1e-12)
      return 7.5+25.2*x;
    
    return (-(6-21*x+26*x*x-19*x*x*x+2*x*x*x*x)/(x*x*x*(1-x)*(1-x)*(1-x)*(1-x)))-(6.0/(x*x*x*x))*log(1-x);    
  }
  
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
