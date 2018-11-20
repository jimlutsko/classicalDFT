#ifndef __LUTSKO_WEIGHTED_DENSITIES__
#define __LUTSKO_WEIGHTED_DENSITIES__

//#include "spherical.h"

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <complex>

#include "FMT_Weighted_Density.h"

#include "Density.h"
//#include "perturbative/PotentialSplitter.h"
//#include "perturbative/HardSphereMapper.h"

/**
  *  @brief Point Class: holds x,y,z coordinates
  *
  *  @detailed  This class is just a triple with members named X,Y and Z
  */  


class Point
{
 public:
  /**
  *   @brief  Default  constructor for Point: all corrdinates are set to zero
  *  
  *   @return nothing 
  */  
 Point() : x_(0.0),y_(0.0),z_(0.0) {};

  /**
  *   @brief  Constructor for Point: 
  *  
  *   @param  x is a coordinate
  *   @param  y is a coordinate
  *   @param  z is a coordinate
  *   @return nothing 
  */  
 Point(double x, double y, double z) : x_(x),y_(y),z_(z) {};

  /**
  *   @brief  Accessor for x coordinate
  *  
  *   @return x_ 
  */  
  double X() const {return x_;}

  /**
  *   @brief  Accessor for y coordinate
  *  
  *   @return y_
  */  
  double Y() const {return y_;}

  /**
  *   @brief  Accessor for z coordinate
  *  
  *   @return z_ 
  */  
  double Z() const {return z_;}

 protected:
  double x_;  ///< x coordinate
  double y_;  ///< y coordinate
  double z_;  ///< z coordinate
};

/**
  *  @brief FMT Class
  *
  *  @detailed  This class holds two 11-dimentional arrays called d0_ and d1_. These in turn hold the weighted densities as
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
  *   @param  dx,dy,dz are lattice spacings
  *   @param  density is the Density object
  *   @param  hsd is the hard sphere diameter
  *   @param  pointsFile contains the points for spherical integration
  *   @param  hsd1 is the hard-sphere diameter for the second species
  *   @return nothing 
  */  
  FMT(Lattice &lattice, double hsd, string& pointsFile, double hsd1 = -1);
  /**
  *   @brief  Default  destrctur for FMT 
  *  
  *   @return nothing 
  */  
  ~FMT(){}

  /**
  *   @brief  EtaMax is the maximum value of eta for which the f1,f2_,f3_ functions are calculated "honestly". For larger values of eta, 
  *           the original, divergent forms are replaced by non-divergent quadratic forms. This is intended to allow for calculations without backtracking. It is only implemented for some models. 
  *  
  *   @param  etaMax is the threshold beyond which quadratic forms are used
  *   @return nothing 
  */    
  void setEtaMax(double etaMax) { etaMax_ = etaMax;}

  /**
  *   @brief  Performs trilinear interpolation: this is a conienience function that is not necessary internally. 
  *  
  *   @param  p is the point for which we interpolate
  *   @param  p0 is one corner of the cube containing p (the one with the minimum values of x,y,z)
  *   @param  p1 is one corner of the cube containing p (the one with the maximal values of x,y,z)
  *   @param  coeffs recieves the weights given to each of the 8 lattice points surrounding p
  *   @return nothing
  */  
  void interpolate(Point &p, Point &p0, Point &p1, double coeffs[8]);

  /**
  *   @brief  Calculates total free energy and dOmega/dRho(i) for each lattice point using FFT convolutions
  *  
  *   @param  density is current density
  *   @param  mu is the chemical potential
  *   @return total free energy of system
  */  
  double calculateFreeEnergyAndDerivatives_fourier_space1(Density& density, DFT_Vec &dF0);
  //  double calculateFreeEnergyAndDerivatives_fourier_space1(Density& density, vec &dF0, vec &dF1);

  /**
  *   @brief  Calculate phi(i) and some derivative information. 
  *
  *   @detailed This function assumes the weighted densities at each lattice site are available.
  *           It calculates the function PHI(r) at lattice site i and also stores dPHI(i)/deta(i), dPHI/dv(i), etc. 
  *  
  *   @param  i is lattice position
  *   @return phi(i)
  */  
  double dPHI(long i);

  /**
  *   @brief  Accessor for weighted density eta
  *  
  *   @param  pos is lattice position
  *   @param  species is the particle species (either 0 or 1).
  *   @return eta(pos)
  */  

  double getEta(long pos, int species) const { return (species == 0 ? d0_[0].r(pos) : d1_[0].r(pos));}

  /**
  *   @brief  Accessor the array holding the weight corresponding to eta, in fourier space. 
  *  
  *   @param  species is either 0 or 1.
  *   @return wek for the requested species. 
  */  
  const DFT_Vec_Complex& getWEK(int species) const { return (species == 0 ? d0_[0].wk() : d1_[0].wk());}


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
  *   @brief  derivative of @f1 with respect to its argument.
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
  *   @brief  derivative of @f2 with respect to its argument.
  *  
  *   @param  eta is the local volume-average density
  *   @return f2(eta).
  */     
 virtual double f2p_(double eta) const = 0;

 /**
  *   @brief  derivative of @f3 with respect to its argument.
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

 protected:
 /**
  *   @brief  convenience accessor for weighted density Eta
  *  
  *   @param  d is the array of weighted densities
  *   @return d[0]
  */      
 FMT_Weighted_Density & Eta(vector<FMT_Weighted_Density>  &d) { return d[0];}

 /**
  *   @brief  convenience accessor for weighted density S
  *  
  *   @param  d is the array of weighted densities
  *   @return d[1]
  */       
 FMT_Weighted_Density & S(vector<FMT_Weighted_Density>  &d)   { return d[1];}

 /**
  *   @brief  convenience accessor for weighted density v(i)
  *  
  *   @param  i is the index of v
  *   @param  d is the array of weighted densities
  *   @return d[2+i]
  */       
 FMT_Weighted_Density & V(int i, vector<FMT_Weighted_Density>  &d) { return d[2+i];}

 /**
  *   @brief  convenience accessor for weighted density T(i,j). 
  *
  *   @detailed elements 0 and 1 are the scalar densities; 2,3,4 are v; 5,6,7 are T(0,k); 8,9 are T(1,k>=1) and d(10) is T(2,2). Note that T is symmetric.
  *  
  *   @param  i is the first  index of T
  *   @param  j is the second index of T
  *   @param  d is the array of weighted densities
  *   @return (i == 0 ? d[5+j] : d[7+j])
  */       
 FMT_Weighted_Density & T(int i, int j, vector<FMT_Weighted_Density>  &d)
   {
     if(i > j) swap(i,j);
     if(i == 0) return d[5+j];
     else if (i == 1) return d[7+j];
     return d[10];
   }

 /**
  *   @brief  This is a one-time-only evaluation of the numerical approximation to the FMT weight functions. These are all 
  *           functions w_{alpha}(i,j) = w_{alpha}(abs(i-j)). Their evaluation involves real-space integrations for which the 
  *           integration points are given in the file pointsFile. Most of the work occurs via a call to initializeWeightedDensities@
  *
  *   @param  densities: the array of weighted densities
  *   @param  hsd is the hard-sphere diameter
  *   @param  pointsFile is the file holding the integration points for integrating a spherical shell
  *   @param  Nx is the number of lattice points in the x-direction
  *   @param  Ny is the number of lattice points in the y-direction
  *   @param  Nz is the number of lattice points in the z-direction
  */        
 void generateWeights(vector<FMT_Weighted_Density> &densities, double hsd, string& pointsFile, long Nx, long Ny, long Nz);

/**
  *   @brief  This is a one-time-only evaluation of the numerical approximation to the FMT weight functions. These are all 
  *           functions w_{alpha}(i,j) = w_{alpha}(abs(i-j)). Their evaluation involves real-space integrations for which the 
  *           integration points are given in the file pointsFile. 
  *
  *   @param  dd: the array of weighted densities
  *   @param  hsd is the hard-sphere diameter
  *   @param  Nx is the number of lattice points in the x-direction
  *   @param  Ny is the number of lattice points in the y-direction
  *   @param  Nz is the number of lattice points in the z-direction
  *   @param  pointsFile is the file holding the integration points for integrating a spherical shell
  */  
 void initializeWeightedDensities(vector<FMT_Weighted_Density> &dd, double hsd, long Nx, long Ny, long Nz, string& pointsFile);


/**
  *   @brief  This function determines the final weighted densities by summing
  *    over the contributions from each species with the appropriate factors of the hard-sphere diameter.
  *
  *   @param  dd: the array of weighted densities
  *   @param  hsd is the hard-sphere diameter
  *   @param  eta is the accumulated weighted density eta
  *   @param  s0  is the accumulated weighted density s/(hsd*hsd)
  *   @param  s1  is the accumulated weighted density s/(hsd)
  *   @param  s2  is the accumulated weighted density s
  *   @param  v1  is the accumulated weighted density v/(hsd)
  *   @param  v2  is the accumulated weighted density v
  *   @param  T  is the accumulated weighted density T
  */  
 void addWeightedDensityContributions(int i, vector<FMT_Weighted_Density> &dd, double hsd,double &eta,double &s0, double &s1, double &s2, double v1[], double v2[], double T0[3][3]);


/**
  *   @brief  This function calculates dPhi/dEta(i), etc where i labes the species.
  *
  *   @param  dd: the array of weighted densities
  *   @param  hsd is the hard-sphere diameter
  *   @param  eta is the accumulated weighted density eta for species i
  *   @param  s0  is the weighted density s/(hsd*hsd) for species i
  *   @param  s1  is the weighted density s/(hsd) for species i
  *   @param  s2  is the weighted density s for species i
  *   @param  v1  is the weighted density v/(hsd) for species i
  *   @param  v2  is the weighted density v for species i
  *   @param  T  is the weighted density T for species i
  *   @param  v1_v2  is v1 dot v2 for species i
  *   @param  v2_v2  is v2 dot v2 for species i
  *   @param  vTV  is v dot T dot v for species i
  *   @param  T2  is trace(T dot T) for species i
  *   @param  T3  is trace(T dot T dot T) for species i
  *   @param  f2  is result of function f2_
  *   @param  f3  is result of function f3_
  *   @param  vT  is v dot T for species i
  *   @param  TT  is T dot T for species i
  */  
 void add_dPhi_d_WeightedDensity(int i, vector<FMT_Weighted_Density> &dd, double hsd, double eta, double s0, double s1, double s2, double v1[], double v2[], double T0[3][3], double v1_v2, double v2_v2, double vTv, double T2, double T3, double f2, double f3, double vT[], double TT[3][3]);


/**
  *   @brief  This function just sums up phi(i) over the lattice to 
  *           give the total free energy
  *
  *   @param  density is the Density object
  *   @return The total free energy
  */  
 double calculateFreeEnergy(Density& density);


  /**
  *   @brief  Calculate dF(i) = dV * dPHI/drho(i)
  *
  *   @param  dd is the array of weighted densities 
  *   @param  dV is the volume element
  *   @param  dF is the vector to be filled.
  */       
 void calculateFreeEnergyDerivatives(vector<FMT_Weighted_Density> &dd, double dV, DFT_Vec &dF);

 /**
  *   @brief  name of model implemented by this class
  *
  *   @return the name as a string
  */        
 virtual string Name() const = 0;

 protected:

 double dx_; ///< spacing of lattice in x direction
 double dy_; ///< spacing of lattice in y direction
 double dz_; ///< spacing of lattice in z direction
 
 double hsd0_; ///< hard sphere diameter for first species
 double hsd1_; ///< hard sphere diameter for second species
 
 vector<FMT_Weighted_Density>  d0_; ///< all weighted densities in real & fourier space for first species
 vector<FMT_Weighted_Density>  d1_; ///< all weighted densities in real & fourier space for second species
 
 DFT_FFT dPhi_; ///< dPHI/drho(i)

 double etaMax_; ///< cutoff used to control divergences
 
};

/**
  *  @brief  The original White Bear  FMT model 
  *
  *  @detailed White-bear FMT model with tensor densities and the CS equation of state
  */  
class WhiteBearI : public FMT
{
 public:
 WhiteBearI(Lattice &lattice, double hsd, string& pointsFile, double hsd1 = -1) 
   : FMT(lattice, hsd, pointsFile, hsd1){};


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

 virtual string Name() const { return string("WhiteBear I");}
};

/**
  *  @brief  The RSLT positive-definite FMT model (pre-White Bear)
  *
  *  @detailed This class implements the RSLT FMT model which is not as accurate as WhiteBear but which is positive definite and therefore stable.
  */  
class RSLT : public FMT
{
 public:
 RSLT(Lattice &lattice, double hsd, string& pointsFile, double hsd1 = -1) 
   : FMT(lattice, hsd, pointsFile, hsd1){};


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
  *  @detailed This class implements a modificaiton of the RSLT FMT model where I use the same basic idea to create a different type of model. In RSLT, the Rosenfeld numerator of s^3-3sv^2 (which is not positive definite) is replaced by s^3(1-v^2/s^2)^3 which agrees in the first two terms and is positive. This in fact works for  s^3(1-(3/n)v^2/s^2)^n for any n>=3. I used this class to experiment with different n but this was for research purposes only and should not be used without a thorough review ...  
  */  


class RSLT2: public RSLT
{
 public:
 RSLT2(Lattice &lattice, double hsd, string& pointsFile, double hsd1 = -1) 
   : RSLT(lattice, hsd, pointsFile, hsd1){};


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
  *  @brief  WhiteBear mark II
  *
  *  @detailed Modified WhiteBear which gives better vacancy densities for solid phase.
  */

class WhiteBearII : public WhiteBearI
{
 public:
 WhiteBearII(Lattice &lattice, double hsd, string& pointsFile, double hsd1 = -1) 
   : WhiteBearI(lattice, hsd, pointsFile, hsd1){};

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
