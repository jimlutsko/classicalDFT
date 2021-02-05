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

#include "Species.h"

// There is some sort of a conflict between Boost headers and <complex> so this needs to come last in order to compile.
#include <complex>


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

  Interaction_Base(Species *s1, Species *s2, Potential1 *v, double kT);
  Interaction_Base() {}

  // Initialization function: must be called before using the object.
  virtual void initialize();

  // Scale the attractive interaction: i.e. same as multiplying the potential by a factor.
  virtual void scale_interaction(double scale_fac)
  {
    a_vdw_ *= scale_fac;
    w_att_.MultBy(scale_fac);  
  }
  
  // This executes the basic functionality of computing energies and forces
  virtual double getInteractionEnergyAndForces();

  // Calculates (d2F/dn_i dn_j)v_j
  void add_second_derivative(vector<DFT_FFT> &v, vector<DFT_Vec> &d2F);

  // This is only used for testing add_second_derivative.
  double second_derivative_brute_force(int I[3], int J[3], vector<DFT_FFT> &v);

  // An internal, debugging function. It should be suppressed at some point.
  double checkCalc(int jx, int jy, int jz);

  // Returns the requested entry in the array of weights
  double getW(long s)  { return w_att_.Real().get(s);}

  // Calculates the mean-field contribution (divided by kT) to the bulk chemical potential for a given species
  double Mu(const vector<double> &x, int species) const
  {
    double mu = 0.0;

    if(s1_->getSequenceNumber() == species)
      mu += 0.5*a_vdw_*x[s2_->getSequenceNumber()];

    if(s2_->getSequenceNumber() == species)
      mu += 0.5*a_vdw_*x[s1_->getSequenceNumber()];

    return mu;
  }

  // Calculates the mean-field contribution (divided by kT) to the bulk free energy per unit volume for this object
  double Fhelmholtz(const vector<double> &x) const {return 0.5*a_vdw_*x[s1_->getSequenceNumber()]*x[s2_->getSequenceNumber()];}  

  // Returns the vDW parameter calculated from the weights
  double getVDWParameter() const { if(!initialized_) throw std::runtime_error("Interaction object must be initialized before calling getVDWParameter()"); return a_vdw_;}
  
 protected:

  virtual void generateWeights();  

  // Calculates the weight w(Sx,Sy,Sz)
  virtual double generateWeight(int Sx, int Sy, int Sz, double dx) = 0;

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & s1_;
    ar & s2_;
    ar & a_vdw_;
    ar & w_att_;
    ar & initialized_;
    ar & v_;
    ar & kT_;
  }

 protected:
  Species *s1_; ///< First of the interacting species
  Species *s2_; ///< Second of the interacting species

  double a_vdw_;  ///< vDW parameter calculted by summing the weights and dividing by temperature

  DFT_FFT w_att_; ///< The weights

  bool initialized_ = false; ///< Flag to make sure the object has been initialized
  Potential1 *v_;            ///< The interatomic potential
  double kT_;                ///< The temperature
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
   *   @param pointsFile: the name of the file from which the spherical integration points are to be read.
   */  
 Interaction(Species *s1, Species *s2, Potential1 *v, double kT, string pointsFile) :
  Interaction_Base(s1,s2,v,kT), pointsFile_(pointsFile) {};

  Interaction() :
    Interaction_Base(){};

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
   */        
  virtual void generateWeights();

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Interaction_Base>(*this);
    ar & weightsFile_;
    ar & pointsFile_;
    
    boost::serialization::void_cast_register<Interaction, Interaction_Base>(static_cast<Interaction *>(NULL),static_cast<Interaction_Base *>(NULL));        
  }
  
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
   *   @param Ngauss: number of gauss-lengendre points in each direction
   */    
 Interaction_Gauss(Species *s1, Species *s2, Potential1 *v, double kT, int Ngauss) :
  Interaction_Base(s1,s2,v,kT)
  {
    gsl_integration_glfixed_table *tr = gsl_integration_glfixed_table_alloc(Ngauss);
    for(int i=0;i<Ngauss;i++)
      {
	double p,w;
	gsl_integration_glfixed_point(0, 1, i, &p, &w, tr);
	gauss_p.push_back(p); gauss_w.push_back(w);
      }
  }

  Interaction_Gauss():
    Interaction_Base() {}  

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
  virtual double generateWeight(int Sx, int Sy, int Sz, double dx);
  /**
   *   @brief  Returns the kernel of the integral being evaluated
   *  
   *   @param Sx,Sy,Sz: the cell 
   *   @param v: the potential
   *   @param x,y,z: the position within the cell
   *
   *   @returns the value of the kernel in cell (Sx,Sy,Sz) at position (x,y,z).
   */        
  virtual double getKernel(int Sx, int Sy, int Sz, double dx, double x, double y, double z) = 0;

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Interaction_Base>(*this);
    ar & gauss_p;
    ar & gauss_w;

    boost::serialization::void_cast_register<Interaction_Gauss, Interaction_Base>(static_cast<Interaction_Gauss *>(NULL),static_cast<Interaction_Base *>(NULL));            
  }
  
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
   *   @param Ngauss: number of gauss-lengendre points in each direction
   */    
 Interaction_Gauss_E(Species *s1, Species *s2, Potential1 *v, double kT, int Ngauss) :
  Interaction_Gauss(s1,s2,v,kT,Ngauss){}

  Interaction_Gauss_E():
    Interaction_Gauss(){}  

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Interaction_Gauss>(*this);
    boost::serialization::void_cast_register<Interaction_Gauss_E, Interaction_Gauss>(static_cast<Interaction_Gauss_E *>(NULL),static_cast<Interaction_Gauss *>(NULL));                
  }
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
  double getKernel(int Sx, int Sy, int Sz, double dx, double x, double y, double z);
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
   *   @param Ngauss: number of gauss-lengendre points in each direction
   */    
 Interaction_Gauss_F(Species *s1, Species *s2, Potential1 *v, double kT, int Ngauss) :
  Interaction_Gauss(s1,s2,v,kT,Ngauss){}

  Interaction_Gauss_F():
    Interaction_Gauss(){}

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Interaction_Gauss>(*this);
    boost::serialization::void_cast_register<Interaction_Gauss_F, Interaction_Gauss>(static_cast<Interaction_Gauss_F *>(NULL),static_cast<Interaction_Gauss *>(NULL));                    
  }

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
  double getKernel(int Sx, int Sy, int Sz, double dx, double x, double y, double z);
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
   */      
 Interaction_Interpolation(Species *s1, Species *s2, Potential1 *v, double kT): 
  Interaction_Base(s1,s2,v,kT) {}

  Interaction_Interpolation():
  Interaction_Base() {}

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
  virtual double generateWeight(int Sx, int Sy, int Sz, double dx);

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Interaction_Base>(*this);
    ar & vv_;
    ar & pt_;

    boost::serialization::void_cast_register<Interaction_Interpolation, Interaction_Base>(static_cast<Interaction_Interpolation *>(NULL),static_cast<Interaction_Base *>(NULL));                    
  }
  
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
 Interaction_Interpolation_Zero(Species *s1, Species *s2, Potential1 *v, double kT) :
  Interaction_Interpolation(s1,s2,v,kT) { vv_.push_back(1.0); pt_.push_back(0.0); }

 Interaction_Interpolation_Zero() :
   Interaction_Interpolation() {};

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Interaction_Interpolation>(*this);
    boost::serialization::void_cast_register<Interaction_Interpolation_Zero, Interaction_Interpolation>(static_cast<Interaction_Interpolation_Zero *>(NULL),static_cast<Interaction_Interpolation *>(NULL));                        
  }
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
   */      
 Interaction_Interpolation_LE(Species *s1, Species *s2, Potential1 *v, double kT) :
  Interaction_Interpolation(s1,s2,v,kT)
    { vv_.push_back( 1.0); vv_.push_back(26.0); vv_.push_back(66.0); vv_.push_back(26.0); vv_.push_back(1.0);
      pt_.push_back(-2.0); pt_.push_back(-1.0); pt_.push_back( 0.0); pt_.push_back( 1.0); pt_.push_back(2.0);
      for(auto &x: vv_) x /= 120.0;
    }

  Interaction_Interpolation_LE() : 
    Interaction_Interpolation(){}


  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Interaction_Interpolation>(*this);
    boost::serialization::void_cast_register<Interaction_Interpolation_LE, Interaction_Interpolation>(static_cast<Interaction_Interpolation_LE *>(NULL),static_cast<Interaction_Interpolation *>(NULL));                            
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
   */      
 Interaction_Interpolation_QE(Species *s1, Species *s2, Potential1 *v, double kT) :
  Interaction_Interpolation(s1,s2,v,kT)
    { vv_.push_back(-1.0); vv_.push_back(8.0);  vv_.push_back(18.0); vv_.push_back(112.0); vv_.push_back(86.0); vv_.push_back(112.0); vv_.push_back(18.0); vv_.push_back(8.0); vv_.push_back(-1.0);
      pt_.push_back(-2.0); pt_.push_back(-1.5); pt_.push_back(-1.0); pt_.push_back(-0.5);  pt_.push_back(0.0);  pt_.push_back(0.5);   pt_.push_back(1.0);  pt_.push_back(1.5); pt_.push_back(2.0);
      for(auto &x: vv_) x /= 360.0;
    }

  Interaction_Interpolation_QE() : 
    Interaction_Interpolation(){}  

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Interaction_Interpolation>(*this);
    boost::serialization::void_cast_register<Interaction_Interpolation_QE, Interaction_Interpolation>(static_cast<Interaction_Interpolation_QE *>(NULL),static_cast<Interaction_Interpolation *>(NULL));                            
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
   */      
 Interaction_Interpolation_LF(Species *s1, Species *s2, Potential1 *v, double kT) :
  Interaction_Interpolation(s1,s2,v,kT)
    { vv_.push_back(1.0);  vv_.push_back(4.0); vv_.push_back(1.0);
      pt_.push_back(-1.0); pt_.push_back(0.0); pt_.push_back(1.0);
      for(auto &x: vv_) x /= 6;
    }

  Interaction_Interpolation_LF() : 
    Interaction_Interpolation(){}
  
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Interaction_Interpolation>(*this);
    boost::serialization::void_cast_register<Interaction_Interpolation_LF, Interaction_Interpolation>(static_cast<Interaction_Interpolation_LF *>(NULL),static_cast<Interaction_Interpolation *>(NULL));
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
   */      
 Interaction_Interpolation_QF(Species *s1, Species *s2, Potential1 *v, double kT) :
  Interaction_Interpolation(s1,s2,v,kT)
    { vv_.push_back(1.0/3);  vv_.push_back(1.0/3); vv_.push_back(1.0/3);
      pt_.push_back(-0.5); pt_.push_back(0.0); pt_.push_back(0.5);
    }  
  Interaction_Interpolation_QF() : 
    Interaction_Interpolation(){}
  
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Interaction_Interpolation>(*this);
    boost::serialization::void_cast_register<Interaction_Interpolation_QF, Interaction_Interpolation>(static_cast<Interaction_Interpolation_QF *>(NULL),static_cast<Interaction_Interpolation *>(NULL));    
  }  
};



class Interaction_Gaussian_Density : public Interaction_Interpolation_Zero
{
 public:
  Interaction_Gaussian_Density(FMT_Gaussian_Species *s1, FMT_Gaussian_Species *s2, Potential1 *v, double kT, int n_gauss_legendre = 10);

 Interaction_Gaussian_Density() :
   Interaction_Interpolation_Zero() {};  

  // Initialization function: must be called before using the object.
  virtual void initialize();
  
  // Scale the attractive interaction: i.e. same as multiplying the potential by a factor.
  virtual void scale_interaction(double scale_fac)
  {
    Interaction_Interpolation_Zero::scale_interaction(scale_fac);
    throw std::runtime_error("scaling not implemented in Interaction_Gaussian_Density");
  }
  
  //  executes the basic functionality of computing energies and forces
  virtual double getInteractionEnergyAndForces();


  double Klm(double al, double am, double A, double r) const;
  double dKlm_dal(double al, double am, double A, double r) const;
  double dKlm_dr(double al, double am, double A, double r) const;

  double Ell(double al, double A) const;
  double dEll_dal(double al, double A) const;
  
  double Elm(double al, double am, double Rlm, double A) const;
  double dElm_dal(double al, double am, double Rlm, double A) const;
  double dElm_dRlm(double al, double am, double Rlm, double A) const;  
    
  
protected:
  vector<double> r_;  // Gauss-Legendre points
  vector<double> w_;  // Gauss-Legendre weights
};

#endif // __LUTSKO__INTERACTION__
