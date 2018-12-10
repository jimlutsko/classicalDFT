#ifndef __LUTSKO_DFT__
#define __LUTSKO_DFT__

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <complex>




#include "Enskog.h"

#include "Potential1.h"
#include "VDW1.h"
#include "Density.h"
#include "FMT.h"


/*! \mainpage classicalDFT: Finite Temperature Density Functional Theory in 3 dimensions
 *
 * \section intro_sec Introduction
 *
 * This code is a library that can be used to build practical applications of ftDFT. It supports the most common Fundamental Measure Theory DFT models for hard spheres (including White Bear) and can be used for crystallization. At present, it only handles single-component systems: generalization to multicomponent systems is in the works. 
 *
 * \section install_sec Installation
 *
 *  The library just needs to be deposited into a directory. Applications should be built in sub-directories of this main directory.
 *
 * etc...
 */


/**
  *  @brief DFT base Class. This class enapsulates the DFT models (e.g. Fundamental measure Theory - FMT; Mean Field Theories; etc.)
  *
  *  @detailed This is the base class for all of the DFT objects.
  */  

class DFT
{
 public:
  /**
  *   @brief  Default  constructor for DFT 
  *  
  *   @return nothing 
  */  
  DFT(){}
  /**
  *   @brief  Default  destructor for DFT
  *  
  *   @return nothing 
  */  
  ~DFT(){}


  double getEta(int i, int species) const { return 0;}

/**
  *   @brief  Calculates the ideal gas contribution to the free energy. Also calculates forces if dF.size >= Density.Ntot.
  *  
  *   @param  density is current density
  *   @param  force vector
  *   @return total free energy of system
  */  
  double  F_IdealGas(Density& density, DFT_Vec& dF);


/**
  *   @brief  Calculates the external potential to the free energy. Also calculates forces if dF.size >= Density.Ntot.
  *  
  *   @param  density is current density
  *   @param  mu is the chemical potential
  *   @param  force vector
  *   @return total free energy of system
  */  
  double F_External(Density& density, double mu, DFT_Vec& dF);

  
/**
  *   @brief  Calculates total grand canonical free energy and dOmega/dRho(i) for each lattice point using FFT convolutions
  *  
  *   @param  density is current density
  *   @param  mu is the chemical potential
  *   @return total free energy of system
  */  
  virtual double calculateFreeEnergyAndDerivatives(Density& density, double mu, DFT_Vec& dF, bool onlyFex) = 0;

  /**
  *   @brief  Compute chemical potential/kT for a uniform system with the given density
  *  
  *   @param  x is the density
  *   @return mu/kT
  */   
  virtual double Mu(double x) const  = 0;

/**
  *   @brief  Compute grand potenial/kT/V for a uniform system with the given density
  *  
  *   @param  x is the density
  *   @return Omega/(kT * V)
  */   
  virtual double Omega(double x) const = 0;

  /**
  *   @brief  Compute Helmholtz free energy/kT/V for a uniform system with the given density
  *  
  *   @param  x is the density
  *   @return F/(kT * V)
  */   
  virtual double Fhelmholtz(double x) const = 0;

  /**
   *   @brief  Returns the name of the DFT object - i.e. identifies the model being used.
   *   @return the name as a string object.
   */     
  virtual string Name() const = 0;

  /**
   *   @brief  This conveys information to the FMT object for models that use it.
   */     
  virtual void setEtaMax(double etaMax) {}

  /**
   *   @brief  Only implemented for DFT_FMT object.
   */     
  virtual double Xliq_From_Mu(double mu) const {throw std::runtime_error("Not implemented");}

};


/**
  *  @brief This represents an ideal gas. 
  *
  *  @detailed It is for testing purposes only.
  */  

class DFT_IdealGas : public DFT
{
 public:
  /**
  *   @brief  Default  constructor for DFT_IdealGas 
  *  
  *   @return nothing 
  */  
  DFT_IdealGas(){}
  /**
  *   @brief  Default  destructor for DFT_IdealGas
  *  
  *   @return nothing 
  */  
  ~DFT_IdealGas(){}


  double getEta(int i, int species) const { return 0;}

/**
  *   @brief  Calculates total grand canonical free energy and dOmega/dRho(i) for each lattice point using FFT convolutions
  *  
  *   @param  density is current density
  *   @param  mu is the chemical potential
  *   @return total free energy of system
  */  
  virtual double calculateFreeEnergyAndDerivatives(Density& density, double mu, DFT_Vec& dF, bool onlyFex = false);
/**
  *   @brief  Compute chemical potential/kT for given density
  *  
  *   @param  x is the density
  *   @return mu/kT
  */   
  virtual double Mu(double x) const {return log(x);}

/**
  *   @brief  Compute grand potential/kT/V for given density
  *  
  *   @param  x is the density
  *   @return Omega/(kT * V)
  */   
  virtual double Omega(double x) const {return -x;}

  /**
  *   @brief  Compute Helmholtz free energy/kT/V for given density
  *  
  *   @param  x is the density
  *   @return Omega/(kT * V)
  */   
  virtual double Fhelmholtz(double x) const {return x*log(x)-x;}

  virtual string Name() const { return string("DFT_IdealGas");}


  double getDF_DRHO(Density &density, double mu, DFT_Vec &dF);
};



/**
  *  @brief DFT_FMT Class: A single-species hard-sphere object.
  *
  *   @detailed A single-species hard-sphere object.
  *  
  */  

template <class T> class DFT_FMT : public DFT
{
 public:
  /**
  *   @brief  Default  constructor for DFT_FMT
  *
  *   @param  density is the Density object
  *   @param  pointsFile gives the location of the file containing the spherical-integration points
  *   @param  hsd is hard sphere diameter
  *   @return nothing 
  */  
  DFT_FMT(Lattice &lattice, FMT_Species  *species);
  /**
  *   @brief  Default  destructor for DFT_FMT
  *  
  *   @return nothing 
  */  
  ~DFT_FMT(){}

/**
  *   @brief  Calculates total grand canonical free energy and dOmega/dRho(i) for each lattice point using FFT convolutions
  *  
  *   @param  density is current density
  *   @param  mu is the chemical potential
  *   @return total free energy of system
  */  
  virtual double calculateFreeEnergyAndDerivatives(Density& density, double mu, DFT_Vec& dF, bool onlyFex = false);
/**
  *   @brief  Compute chemical potential/kT for given density
  *  
  *   @param  x is the density `
  *   @return mu/kT
  */   
  virtual double Mu(double x) const { Enskog en(x); return en.chemPotentialCS();}

/**
  *   @brief  Find liquid with given chem potential/kT
  *  
  *   @param  mu is chem potential
  *   @return mu/kT
  */   
  virtual double Xliq_From_Mu(double mu) const;

/**
  *   @brief  Compute grand potenial/kT/V for given density
  *  
  *   @param  x is the density
  *   @return Omega/(kT * V)
  */   
  virtual double Omega(double x) const {Enskog en(x); return x*(en.freeEnergyCS() - en.chemPotentialCS());}

/**
  *   @brief  Compute Helmhltz Free Energy/kT/V for given density
  *  
  *   @param  x is the density
  *   @return F/(kT * V)
  */   
  virtual double Fhelmholtz(double x) const {Enskog en(x); return x*en.freeEnergyCS();}

/**
  *   @brief  Accessor for hard sphere diameter
  *  
  *   @param  none
  *   @return hard sphere diameter
  */   
  double HSD() const { return hsd_;}

/**
  *   @brief  Accessor for weighted density eta: local packing fraction
  *  
  *   @param  position in array
  *   @return eta
  */   
  double getEta(int i, int species) const { return fmt_.getEta(i,species);}

  /**
  *   @brief  Accessor for wek on FMT object
  *  
  *   @param  pos is lattice position
  *   @return eta(pos)
  */  
 const DFT_Vec_Complex& getWEK(int species) const { return fmt_.getWEK(species);}

 virtual string Name() const { return string("DFT_FMT : ") + fmt_.Name();}

 /**
  *   @brief  This conveys information to the FMT object for models that use it.
  */     
 virtual void setEtaMax(double etaMax) {fmt_.setEtaMax(etaMax);}

 /**
  *   @brief  Accessor for real-space FMT vector density
  *  
  *   @param  J is cartesian component
  *   @return v - real space
  */  
 const DFT_Vec &getV_Real(int J, int species = 0) const { return fmt_.getV_Real(J,species);}

   /**
  *   @brief  Accessor for kernal for FMT vector density
  *  
  *   @param  J is cartesian component
  *   @return wv - infourier space
  */  
 const DFT_Vec_Complex &getVweight_Four(int J, int species = 0) const { return fmt_.getVweight_Four(J,species);}

 protected:
  T         fmt_;   ///< Hard-sphere FMT object
  double    hsd_;       ///< Hard sphere diameter
};



/**
  *  @brief FMT plus van der Waals mean field tail.
  *
  *   @detailed A single-species mean-field DFT model. It holds a DFT_FMT object to do the hard-sphere stuff.
  *  
  */  

template <class T> class DFT_VDW : public DFT
{
 public:
  /**
  *   @brief  Default  constructor for DFT 
  *  
  *   @param  density is the Density object
  *   @param  potential is the particle-particle potential
  *   @param  pointsFile gives the location of the file containing the spherical-integration points
  *   @return nothing 
  */  
  DFT_VDW(Lattice &lattice, Potential1 &potential,  string& pointsFile, double kT);
  /**
  *   @brief  Default  destructor for DFT : deletes  fmt_ object. 
  *  
  *   @return nothing 
  */  
  ~DFT_VDW(){if(dft_fmt_) delete dft_fmt_; if(species_) delete species_;}

/**
  *   @brief  Calculates total grand canonical free energy and dOmega/dRho(i) for each lattice point using FFT convolutions
  *  
  *   @param  density is current density
  *   @param  mu is the chemical potential
  *   @return total free energy of system
  */  
  virtual double calculateFreeEnergyAndDerivatives(Density& density, double mu, DFT_Vec& dF, bool onlyFex = false);

/**
  *   @brief  Compute chemical potential/kT for given density
  *  
  *   @param  x is the density
  *   @return mu/kT
  */   
  virtual double Mu(double x) const { return vdw_.chemPotential(x);}

/**
  *   @brief  Compute grand potenial/kT/V for given density
  *  
  *   @param  x is the density
  *   @return Omega/(kT * V) = -Pressure/kT
  */   
  virtual double Omega(double x) const {return -vdw_.pressure(x);}

  /**
  *   @brief  Compute helmholtz free energy/kT/V for given density
  *  
  *   @param  x is the density
  *   @return F/(kT * V)
  */   
  
  virtual double Fhelmholtz(double x) const {return vdw_.helmholtzPerUnitVolume(x);}

/**
  *   @brief  Get coexisting liquid and gas
  *  
  *   @param  references to liq and gas densities to be filled in. 
  *   @return gas density in program units
  */   
  void coexistence(double& xliq, double &xgas) const {vdw_.findCoexistence(xliq,xgas);}

/**
  *   @brief  Find spinodals
  *  
  *   @param  references to spinodal densities to be filled in. 
  *   @return nothing
  */   
  void spinodal(double& xs1, double &xs2) const {vdw_.spinodal(xs1,xs2);}

  //  void spinodals(double& x1, double& x2) const ;

/**
  *   @brief  find liquid density with given chemical potential/kT
  *  
  *   @param  mu is chem potential / kT
  *   @param  mu_coex is chem potential / kT for coexistence
  *   @param  xliq_coex  is coexisting liquid density
  *   @return gas density in program units
  */   
  double findLiquidFromMu(double mu, double mu_coex, double xliq_coex) const
  { return vdw_.findLiquidFromMu(mu,mu_coex,xliq_coex);}

/**
  *   @brief  Accessor for coexisting hsd
  *  
  *   @param  none
  *   @return hsd
  */   
  double HSD() const { return dft_fmt_->HSD();}

  /**
  *   @brief  Accessor for coexisting vdw parameter
  *  
  *   @param  none
  *   @return vdw_.a_
  */   
  double get_VDW_Param() const { return vdw_.get_VDW_Parameter();}
  
  virtual string Name() const { return string("DFT_VDW : ") + dft_fmt_->Name();}
  
  /**
   *   @brief  This conveys information to the FMT object for models that use it.
   */     
  virtual void setEtaMax(double etaMax) {dft_fmt_->setEtaMax(etaMax);}

  
 protected:
  DFT_FFT w_att_;
  DFT_FFT v_mean_field_;

  FMT_Species *species_;
  DFT_FMT<T> *dft_fmt_; ///< The hard-sphere dft 
  VDW1 vdw_;
};

/**
  *  @brief DFT_VDW Class
  *
  *   @detailed A single-species mean-field DFT model. It holds a DFT_FMT object to do the hard-sphere stuff.
  *  
  */  

template <class T> class DFT_VDW_Surfactant : public DFT_VDW<T>
{
 public:
  /**
  *   @brief  Default  constructor for DFT 
  *  
  *   @param  density is the Density object
  *   @param  potential is the particle-particle potential
  *   @param  pointsFile gives the location of the file containing the spherical-integration points
  *   @return nothing 
  */  
  DFT_VDW_Surfactant(Lattice &lattice, Potential1 &potential,  string& pointsFile, double kT, double Asurf, double rhosurf);
  /**
  *   @brief  Default  destructor for DFT : deletes  fmt_ object. 
  *  
  *   @return nothing 
  */  
  ~DFT_VDW_Surfactant(){}

  virtual double calculateFreeEnergyAndDerivatives(Density& density, double mu, DFT_Vec& dF, bool onlyFex = false);

  void setSurfactant(double rhos, double A) {rho_surf_ = rhos; Asurf_ = A;}
  void setFixedN(bool flag) {bFixedN_ = flag;}
  
  double getSurfactant(long i, Density &density) const {return surfactant_density_.cReal().get(i);}
  
  virtual double Mu(double x) const;
  virtual double Omega(double x) const;
  virtual double Fhelmholtz(double x) const; 
  void coexistence(double& xliq, double &xgas) const;
  void spinodal(double& xs1, double &xs2) const; 
  double findLiquidFromMu(double mu, double mu_coex, double xliq_coex) const;


  virtual string Name() const { return string("DFT_VDW_Surfactant : ") +  DFT_VDW<T>::dft_fmt_->Name();}
  
 protected:
  DFT_FFT surfactant_density_;   //< Arrays holding actual surfactant density and its FFT
  DFT_FFT surfactant_potential_; //< Arrays holding surfactant assymetric potential
  double Asurf_;                     //< strength of assymetric (v^2 term) surfactant interaction
  double rho_surf_;              //< parameter rho_0
  double ax_;                    //< VDW parameter for symmetric part of surfactant interaction
  bool bFixedN_;                 ///< Fixed particle number flag
};

#endif // __LUTSKO_DFT__
