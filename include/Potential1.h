#ifndef __LUTSKO_POTENTIAL1__
#define __LUTSKO_POTENTIAL1__

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/export.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

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
  Potential1(double sigma, double eps, double rcut);
  Potential1();
    
  void set_WCA_limit(double r); ///< Distance from origin that attractive continuate extends
  void setBH();                 ///< Use the BH split of the potential
  
  double V(double r)                 const;    ///< The cut and shifted potential at point r
  double V2(double r2)               const;  ///< The cut and shifted potential at point r calculated from r^2

  double getRcut()                   const;
  double getR0()                     const;  
  double Watt(double r)              const; ///< The attractive tail 
  double Watt2(double r2)            const; ///< Attractive part calcuated from r2
  double V0(double r)                const; ///< The repulsive part of the potential
  double getHSD(double kT)           const; ///< Calculates HSD by numeric integration. 
  double getVDW_Parameter(double kT) const; ///< Calculates the vDW parameter by numeric integration

  virtual double getHardCore()   const = 0; ///< Hard core diameter of the original potential (e.g. zero for LJ).   
  virtual double getRmin()       const = 0; ///< min of potential  
  virtual string getIdentifier() const = 0; ///< Unique name of this potential
  
 protected:
  virtual double vr(double r)   const = 0; ///< The underlying potential
  virtual double vr2(double r2) const = 0; ///< The underlying potential

  virtual double dBH_Kernal(double r) const {return (1.0-exp(-V0(r)/kT_));}  ///< Kernal for calcuating hsd
  virtual double vdw_Kernal(double r) const {return r*r*Watt(r);}  ///< Kernel for calculating vDW parameter

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & sigma_;
    ar & eps_;
    ar & rcut_;
    ar & rmin_;
    ar & shift_;
    ar & Vmin_;
    ar & r_att_min_;
    ar & r0_;
    ar & bhFlag_;
    ar & kT_;    
  }
  
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

// Truncated Lennard-Jones potential
class LJ : public Potential1
{
 public:
  LJ(double sigma, double eps, double rcut);
  LJ();
  
  virtual double getRmin()       const;
  virtual double getHardCore()   const;
  virtual string getIdentifier() const;

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Potential1>(*this);
    boost::serialization::void_cast_register<LJ, Potential1>(static_cast<LJ *>(NULL),static_cast<Potential1 *>(NULL));            
  }  
  
 protected:
  virtual double vr(double r)    const;
  virtual double vr2(double r2)  const;
};

// Truncated ten Wolde-Frenkel potential
class tWF : public Potential1
{
 public:
  tWF(double sigma, double eps, double rcut, double alpha = 50.0);
  tWF();
  
  virtual double getRmin()       const;
  virtual double getHardCore()   const;
  virtual string getIdentifier() const;

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Potential1>(*this);
    ar & alpha_;

    boost::serialization::void_cast_register<tWF, Potential1>(static_cast<tWF *>(NULL),static_cast<Potential1 *>(NULL));            
  }
  
 protected:
  virtual double vr(double r)   const;
  virtual double vr2(double r2) const;

  double alpha_;
};

// Frenkel's potential from https://arxiv.org/pdf/1910.05746.pdf
class WHDF : public Potential1
{
 public:
  WHDF(double sigma, double eps, double rcut);
  WHDF();
  
  virtual double getRmin()       const;
  virtual double getHardCore()   const;
  virtual string getIdentifier() const;
  
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Potential1>(*this);
    boost::serialization::void_cast_register<WHDF, Potential1>(static_cast<WHDF *>(NULL),static_cast<Potential1 *>(NULL));            
  }  

 protected:
  virtual double vr(double r)   const;
  virtual double vr2(double r2) const;

};

#endif // __POTENTIAL
