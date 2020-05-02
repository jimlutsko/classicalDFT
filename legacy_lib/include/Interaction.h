/* This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
 * To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
 *
 * Author: James F. Lutsko
 * www.lutsko.com
 */

#ifndef __LUTSKO__INTERACTION__
#define __LUTSKO__INTERACTION__

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <complex>
#include <complex.h>

#include "Species.h"
#include "Log.h"
#include "myColor.h"

/**
  *  @brief This class encapsulates the interaction between two species (or one species with itself)
  *
  *  @detailed Its main job is to compute an array $w({\mathbf R})$ with which the mean-field contribution to the 
  *         enegy can be calulated as $\sum_{\mathbf R} \sum_{\mathbf R'} \rho({\mathbf R})\rho({\mathbf R'})w({\mathbf R}-{\mathbf R'})$
  *         The array w is computed from the attractive part of the potential and its exact form depends on what approximations are being made.
  *         There are two main classes: The "Energy" class in which w is determined starting with the energy and the "Force" class in which it is 
  *         determined starting with the forces. 
  */  
class Interaction_Base
{
 public:
  /**
   *   @brief  Constructor 
   *  
   *   @param s1: the first species.
   *   @param s2: the second species.
   *   @param kT: the temperature
   *   @param log: the log object
   */  
  Interaction_Base(Species &s1, Species &s2, Potential1 &v, double kT, Log &log);

  /**
   *   @brief  Initialization function: must be called before using the object.
   *  
   */  
  virtual void initialize();

  /**
   *   @brief  This executes the basic functionality of computing energies and forces
   *  
   *   @returns The mean-field contribution to the (total) free energy divided by kT.
   */    
  double getInteractionEnergyAndForces();

  /**
   *   @brief  An internal, debugging function. It should be suppressed at some point.
   *  
   *   @returns 
   */        
  double checkCalc(int jx, int jy, int jz);

  /**
   *   @brief  Returns the requested entry in the array of weights
   *  
   *   @param  s: position of entry requested
   *
   *   @returns  w(s)
   */        
  double getW(long s)  { return w_att_.Real().get(s);}

  /**
   *   @brief  Calculates the mean-field contribution (divided by kT) to the bulk chemical potential for a given species
   *  
   *   @param  x: array holding the density of each species
   *   @param  species: the species whose chemical potential is being calculated
   *
   *   @return the chemical potential divided by kT
   */        
  double Mu(const vector<double> &x, int species) const
  {
    double mu = 0.0;

    if(s1_.getSequenceNumber() == species)
      mu += 0.5*a_*a_vdw_*x[s2_.getSequenceNumber()];

    if(s2_.getSequenceNumber() == species)
      mu += 0.5*a_*a_vdw_*x[s1_.getSequenceNumber()];

    mu += 3*b_*a_vdw_*a_vdw_*x[s1_.getSequenceNumber()]*x[s1_.getSequenceNumber()];

    
    return mu;
  }

  /**
   *   @brief  Calculates the mean-field contribution (divided by kT) to the bulk free energy per unit volume for this object
   *  
   *   @param  x: array holding the density of each species
   *
   *   @return the mean-field contribution to the free energy per unit volume divided by kT
   */  
  double Fhelmholtz(const vector<double> &x) const {return 0.5*a_*a_vdw_*x[s1_.getSequenceNumber()]*x[s2_.getSequenceNumber()] + b_*a_vdw_*a_vdw_*pow(x[s1_.getSequenceNumber()],3);}  

  /**
   *   @brief  Returns the vDW parameter calculated from the weights
   *
   *   @return sum_{\mathbf R} w({\mathbf R})/kT
   */    
  double getVDWParameter() const { if(!initialized_) throw std::runtime_error("Interaction object must be initialized before calling getVDWParameter()"); return a_*a_vdw_;}
  double a_ = 1.0;
  double b_ = 0.0;
  


  
 protected:
  
  /**
   *   @brief  This calculates the array w. 
   *  
   *   @param  v: the interatomic potential
   *   @param  log: the log object for output
   */      
  virtual void generateWeights(Potential1 &v, Log& log);  

  /**
   *   @brief  Calculates the weight w(Sx,Sy,Sz)
   *  
   *   @param Sx,Sy,Sz: the cell 
   *   @param v: the potential
   *   @param dx: lattice spacing
   *
   *   @returns the value of the kernel in cell (Sx,Sy,Sz) without the global dV*dV
   */          
  virtual double generateWeight(int Sx, int Sy, int Sz, double dx, Potential1& v) = 0;


 protected:
  Species &s1_; ///< First of the interacting species
  Species &s2_; ///< Second of the interacting species

  double a_vdw_;  ///< vDW parameter calculted by summing the weights and dividing by temperature

  DFT_FFT w_att_; ///< The weights

  bool initialized_ = false; ///< Flag to make sure the object has been initialized
  Potential1 &v_;            ///< The interatomic potential
  double kT_;                ///< The temperature
  Log& log_;                 ///< Log object
};


/**
  *  @brief This is an old form of the Interaction object that calculates the weights based on spherical integration using points supplied in a file. 
  *         This is only maintained for backward compatability and will be eliminated at some point.
  */  
class Interaction : public Interaction_Base
{
 public:
  /**
   *   @brief  Constructor 
   *  
   *   @param s1: the first species.
   *   @param s2: the second species.
   *   @param kT: the temperature
   *   @param log: the log object
   *   @param pointsFile: the name of the file from which the spherical integration points are to be read.
   */  
 Interaction(Species &s1, Species &s2, Potential1 &v, double kT, Log &log, string& pointsFile) :
  Interaction_Base(s1,s2,v,kT,log), pointsFile_(pointsFile) {};

  /**
   *   @brief  Initialization function: must be called before using the object.
   *  
   */  
  virtual void initialize();

  
 protected:
  /**
   *   @brief  Reads the weights fromt the file with the supplied name prefix (plus stuff added in this function)
   *  
   *   @returns  true, if successful
   */        
  bool readWeights();

  /**
   *   @brief  Checks whether the input file corresponds to this object or not.
   *  
   *   @param in: the ifstream attached to the weight file
   *  
   *   @returns false if the file cannot be used.
   */      
  virtual bool checkWeightsFile(ifstream &in);

  /**
   *   @brief  This calculates the array w. 
   *  
   *   @param  v: the interatomic potential
   *   @param  log: the log object for output
   */        
  virtual void generateWeights(Potential1 &v, Log& log);
  
 protected:
  string pointsFile_; ///< The name of the file with the spherical integation points.
  stringstream weightsFile_; ///< The name of the file with the spherical integation points.
};

/**
  *  @brief Calculates the weight array using Gauss-Legendre integration. Typically, Ngauss = 5 is sufficient. This is an abstract base class.
  */  

class Interaction_Gauss : public Interaction_Base
{
 public:
  /**
   *   @brief  Constructor 
   *  
   *   @param s1: the first species.
   *   @param s2: the second species.
   *   @param kT: the temperature
   *   @param log: the log object
   *   @param Ngauss: number of gauss-lengendre points in each direction
   */    
 Interaction_Gauss(Species &s1, Species &s2, Potential1 &v, double kT, Log &log, int Ngauss) :
  Interaction_Base(s1,s2,v,kT,log)
  {
    gsl_integration_glfixed_table *tr = gsl_integration_glfixed_table_alloc(Ngauss);
    for(int i=0;i<Ngauss;i++)
      {
	double p,w;
	gsl_integration_glfixed_point(0, 1, i, &p, &w, tr);
	gauss_p.push_back(p); gauss_w.push_back(w);
      }
  }

 protected:
  /**
   *   @brief  Calculates the weight w(Sx,Sy,Sz)
   *  
   *   @param Sx,Sy,Sz: the cell 
   *   @param v: the potential
   *   @param dx: lattice spacing
   *
   *   @returns the value of the kernel in cell (Sx,Sy,Sz) without the global dV*dV
   */        
  virtual double generateWeight(int Sx, int Sy, int Sz, double dx, Potential1& v);
  /**
   *   @brief  Returns the kernel of the integral being evaluated
   *  
   *   @param Sx,Sy,Sz: the cell 
   *   @param v: the potential
   *   @param x,y,z: the position within the cell
   *
   *   @returns the value of the kernel in cell (Sx,Sy,Sz) at position (x,y,z).
   */        
  virtual double getKernel(int Sx, int Sy, int Sz, double dx, Potential1 &v, double x, double y, double z) = 0;
  
 protected:
  vector<double> gauss_p;
  vector<double> gauss_w;
};

/**
  *  @brief Calculates the weight array using Gauss-Legendre integration. Typically, Ngauss = 5 is sufficient. This uses the "Energy route".
  */  
class Interaction_Gauss_E : public Interaction_Gauss
{
 public:
  /**
   *   @brief  Constructor 
   *  
   *   @param s1: the first species.
   *   @param s2: the second species.
   *   @param kT: the temperature
   *   @param log: the log object
   *   @param Ngauss: number of gauss-lengendre points in each direction
   */    
 Interaction_Gauss_E(Species &s1, Species &s2, Potential1 &v, double kT, Log &log, int Ngauss) :
  Interaction_Gauss(s1,s2,v,kT,log,Ngauss){}

 protected:
  /**
   *   @brief  Returns the kernel of the integral being evaluated
   *  
   *   @param Sx,Sy,Sz: the cell 
   *   @param v: the potential
   *   @param x,y,z: the position within the cell
   *
   *   @returns the value of the kernel in cell (Sx,Sy,Sz) at position (x,y,z).
   */        
  double getKernel(int Sx, int Sy, int Sz, double dx, Potential1 &p, double x, double y, double z);
};

/**
  *  @brief Calculates the weight array using Gauss-Legendre integration. Typically, Ngauss = 5 is sufficient. This uses the "Force route".
  */  
class Interaction_Gauss_F : public Interaction_Gauss
{
 public:
  /**
   *   @brief  Constructor 
   *  
   *   @param s1: the first species.
   *   @param s2: the second species.
   *   @param kT: the temperature
   *   @param log: the log object
   *   @param Ngauss: number of gauss-lengendre points in each direction
   */    
 Interaction_Gauss_F(Species &s1, Species &s2, Potential1 &v, double kT, Log &log, int Ngauss) :
  Interaction_Gauss(s1,s2,v,kT,log,Ngauss){}

 protected:
  /**
   *   @brief  Returns the kernel of the integral being evaluated
   *  
   *   @param Sx,Sy,Sz: the cell 
   *   @param v: the potential
   *   @param x,y,z: the position within the cell
   *
   *   @returns the value of the kernel in cell (Sx,Sy,Sz) at position (x,y,z).
   */      
  double getKernel(int Sx, int Sy, int Sz, double dx, Potential1 &v, double x, double y, double z);
};

/**
  *  @brief Calculates the weight array using polynomial interpolation of the potential. 
  */  
class Interaction_Interpolation : public Interaction_Base
{
 public:
  /**
   *   @brief  Constructor 
   *  
   *   @param s1: the first species.
   *   @param s2: the second species.
   *   @param kT: the temperature
   *   @param log: the log object
   */      
 Interaction_Interpolation(Species &s1, Species &s2, Potential1 &v, double kT, Log &log): 
  Interaction_Base(s1,s2,v,kT,log) {}

 protected:
  /**
   *   @brief  Calculates the weight w(Sx,Sy,Sz)
   *  
   *   @param Sx,Sy,Sz: the cell 
   *   @param v: the potential
   *   @param dx: lattice spacing
   *
   *   @returns the value of the kernel in cell (Sx,Sy,Sz) without the global dV*dV
   */        
  virtual double generateWeight(int Sx, int Sy, int Sz, double dx, Potential1& v);
  
 protected:
  vector<double> vv_; ///< Weights of the interpolant
  vector<double> pt_; ///< Positions at which the potential must be evaluated
};

/**
  *  @brief Calculates the weight array as the potential evaluated at the lattice points.
  */  
class Interaction_Interpolation_Zero : public Interaction_Interpolation
{
 public:
  /**
   *   @brief  Constructor 
   *  
   *   @param s1: the first species.
   *   @param s2: the second species.
   *   @param kT: the temperature
   *   @param log: the log object
   */    
 Interaction_Interpolation_Zero(Species &s1, Species &s2, Potential1 &v, double kT, Log &log) :
  Interaction_Interpolation(s1,s2,v,kT,log) { vv_.push_back(1.0); pt_.push_back(0.0); }
};

/**
  *  @brief Calculates the weight array using linear interpolation of the potential and the "Energy route". 
  */  
class Interaction_Interpolation_LE : public Interaction_Interpolation
{
 public:
  /**
   *   @brief  Constructor 
   *  
   *   @param s1: the first species.
   *   @param s2: the second species.
   *   @param kT: the temperature
   *   @param log: the log object
   */      
 Interaction_Interpolation_LE(Species &s1, Species &s2, Potential1 &v, double kT, Log &log) :
  Interaction_Interpolation(s1,s2,v,kT,log)
    { vv_.push_back( 1.0); vv_.push_back(26.0); vv_.push_back(66.0); vv_.push_back(26.0); vv_.push_back(1.0);
      pt_.push_back(-2.0); pt_.push_back(-1.0); pt_.push_back( 0.0); pt_.push_back( 1.0); pt_.push_back(2.0);
      for(auto &x: vv_) x /= 120.0;
    }
};

/**
 *  @brief Calculates the weight array using quadratic interpolation of the potential and the "Energy route". 
 */  
class Interaction_Interpolation_QE : public Interaction_Interpolation
{
 public:
  /**
   *   @brief  Constructor 
   *  
   *   @param s1: the first species.
   *   @param s2: the second species.
   *   @param kT: the temperature
   *   @param log: the log object
   */      
 Interaction_Interpolation_QE(Species &s1, Species &s2, Potential1 &v, double kT, Log &log) :
  Interaction_Interpolation(s1,s2,v,kT,log)
    { vv_.push_back(-1.0); vv_.push_back(8.0);  vv_.push_back(18.0); vv_.push_back(112.0); vv_.push_back(86.0); vv_.push_back(112.0); vv_.push_back(18.0); vv_.push_back(8.0); vv_.push_back(-1.0);
      pt_.push_back(-2.0); pt_.push_back(-1.5); pt_.push_back(-1.0); pt_.push_back(-0.5);  pt_.push_back(0.0);  pt_.push_back(0.5);   pt_.push_back(1.0);  pt_.push_back(1.5); pt_.push_back(2.0);
      for(auto &x: vv_) x /= 360.0;
    }
};

/**
 *  @brief Calculates the weight array using linear interpolation of the potential and the "Force route". 
 */  
class Interaction_Interpolation_LF : public Interaction_Interpolation
{
 public:
  /**
   *   @brief  Constructor 
   *  
   *   @param s1: the first species.
   *   @param s2: the second species.
   *   @param kT: the temperature
   *   @param log: the log object
   */      
 Interaction_Interpolation_LF(Species &s1, Species &s2, Potential1 &v, double kT, Log &log) :
  Interaction_Interpolation(s1,s2,v,kT,log)
    { vv_.push_back(1.0);  vv_.push_back(4.0); vv_.push_back(1.0);
      pt_.push_back(-1.0); pt_.push_back(0.0); pt_.push_back(1.0);
      for(auto &x: vv_) x /= 6;
    }
};

/**
 *  @brief Calculates the weight array using guadratic interpolation of the potential and the "Force route". 
 */  
class Interaction_Interpolation_QF : public Interaction_Interpolation
{
 public:
  /**
   *   @brief  Constructor 
   *  
   *   @param s1: the first species.
   *   @param s2: the second species.
   *   @param kT: the temperature
   *   @param log: the log object
   */      
 Interaction_Interpolation_QF(Species &s1, Species &s2, Potential1 &v, double kT, Log &log) :
  Interaction_Interpolation(s1,s2,v,kT,log)
    { vv_.push_back(1.0);  vv_.push_back(1.0); vv_.push_back(1.0);
      pt_.push_back(-0.5); pt_.push_back(0.0); pt_.push_back(0.5);
      for(auto &x: vv_) x /= 3;      
    }
};



#endif // __LUTSKO__INTERACTION__
