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


  void addSpecies(Species* s) {allSpecies_.push_back(s);}
  int getNumberOfSpecies() const {return allSpecies_.size();}

  double getNumberAtoms(int i) const {return allSpecies_[i]->getDensity().getNumberAtoms();}

  const Lattice& lattice() const {return allSpecies_.front()->getLattice();}

  double get_convergence_monitor() const { double d = 0; for(auto& s: allSpecies_) d += s->get_convergence_monitor(); return d;}
  
  const Density& getDensity(int i) const {return allSpecies_[i]->getDensity();}

  DFT_Vec &getDF(int i) {return allSpecies_[i]->getDF();}
  
  void doDisplay(int i, string &title, string &file){ allSpecies_[i]->doDisplay(title,file);}
  void set_density_from_amplitude(int i,DFT_Vec &x) {allSpecies_[i]->set_density_from_amplitude(x);}
  void writeDensity(int i, string &of) const {allSpecies_[i]->getDensity().writeDensity(of);}

  /**
  *   @brief  Calculates the ideal gas contribution to the free energy. Also calculates forces if flag is true.
  *  
  *   @param  bCalcForce is a flag telling whether or not to calculate the forces
  *   @return ideal gas contribution to the free energy
  */  
  double  F_IdealGas(bool bCalcForce = true);


/**
  *   @brief  Calculates the external potential to the free energy. Also calculates forces if dF.size >= Density.Ntot.
  *  
  *   @param  bCalcForce is a flag telling whether or not to calculate the forces
  *   @return total external force contribution to the free energy
  */  
  double F_External(bool bCalcForce = true);

  
/**
  *   @brief  Calculates total grand canonical free energy and dOmega/dRho(i) for each lattice point using FFT convolutions
  *  
  *   @param  onlyFex : if true, the ideal and external contributions are not added
  *   @return total free energy of system
  */  
  virtual double calculateFreeEnergyAndDerivatives(bool onlyFex) = 0;

  /**
  *   @brief  Compute chemical potential/kT for a uniform system with the given density
  *  
  *   @param  x is the density
  *   @return mu/kT
  */   
  virtual double Mu(double x) const  = 0;

/**
  *   @Brief  Compute grand potenial/kT/V for a uniform system with the given density
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

 protected:
  vector<Species*> allSpecies_;

  
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


/**
  *   @brief  Calculates total grand canonical free energy and dOmega/dRho(i) for each lattice point using FFT convolutions
  *  
  *   @param  density is current density
  *   @param  mu is the chemical potential
  *   @return total free energy of system
  */  
  virtual double calculateFreeEnergyAndDerivatives(bool onlyFex = false);
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
  *   @return nothing 
  */  
  
  // External field
  DFT_FMT(int Nx, int Ny, int Nz)
    : fmt_(Nx,Ny,Nz) {}


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
  virtual double calculateFreeEnergyAndDerivatives(bool onlyFex = false);
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
  *   @brief  The name of this DFT object
  *  
  *   @return name as string
  */   
  virtual string Name() const { return string("DFT_FMT : ") + fmt_.Name();}

 /**
  *   @brief  This conveys information to the FMT object for models that use it.
  */     
 virtual void setEtaMax(double etaMax) {fmt_.setEtaMax(etaMax);}

 protected:
  T         fmt_;   ///< Hard-sphere FMT object
};



/**
  *  @brief FMT plus van der Waals mean field tail.
  *
  *   @detailed A single-species mean-field DFT model. It holds a DFT_FMT object to do the hard-sphere stuff.
  *  
  */  

template <class T> class DFT_VDW : public DFT_FMT<T>
{
 public:
  /**
  *   @brief  Default  constructor for DFT 
  *  
  *   @param  density is the Density object
  *   @return nothing 
  */  
  DFT_VDW(int Nx, int Ny, int Nw);
  /**
  *   @brief  Default  destructor for DFT : deletes  fmt_ object. 
  *  
  *   @return nothing 
  */  
  ~DFT_VDW(){} //if(dft_fmt_) delete dft_fmt_;}


  void addSpecies(VDW_Species* species,  double kT);
  
/**
  *   @brief  Calculates total grand canonical free energy and dOmega/dRho(i) for each lattice point using FFT convolutions
  *  
  *   @param  density is current density
  *   @param  mu is the chemical potential
  *   @return total free energy of system
  */  
  virtual double calculateFreeEnergyAndDerivatives(bool onlyFex = false);

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
  double HSD() const { return species_->getHSD();}

  /**
  *   @brief  Accessor for coexisting vdw parameter
  *  
  *   @param  none
  *   @return vdw_.a_
  */   
  double get_VDW_Param() const { return vdw_.get_VDW_Parameter();}
  
  virtual string Name() const { return string("DFT_VDW : ") + DFT_FMT<T>::Name();} //dft_fmt_->Name();}
  
  /**
   *   @brief  This conveys information to the FMT object for models that use it.
   */     
  //  virtual void setEtaMax(double etaMax) {dft_fmt_->setEtaMax(etaMax);}

  
 protected:
  DFT_FFT v_mean_field_;

  VDW_Species *species_;
  //  DFT_FMT<T> *dft_fmt_; ///< The hard-sphere dft 
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
  *   @return nothing 
  */  
  DFT_VDW_Surfactant(int Nx, int Ny, int Nz, double Asurf, double rhosurf);
  /**
  *   @brief  Default  destructor for DFT : deletes  fmt_ object. 
  *  
  *   @return nothing 
  */  
  ~DFT_VDW_Surfactant(){}

  void addSpecies(VDW_Species *species, double kT);

  virtual double calculateFreeEnergyAndDerivatives(bool onlyFex = false);

  void setSurfactant(double rhos, double A) {rho_surf_ = rhos; Asurf_ = A;}
  void setFixedN(bool flag) {bFixedN_ = flag;}
  
  double getSurfactant(long i) const {return surfactant_density_.cReal().get(i);}
  
  virtual double Mu(double x) const;
  virtual double Omega(double x) const;
  virtual double Fhelmholtz(double x) const; 
  void coexistence(double& xliq, double &xgas) const;
  void spinodal(double& xs1, double &xs2) const; 
  double findLiquidFromMu(double mu, double mu_coex, double xliq_coex) const;


  virtual string Name() const { return string("DFT_VDW_Surfactant : ") +  DFT_VDW<T>::Name();} //dft_fmt_->Name();}
  
 protected:
  DFT_FFT surfactant_density_;   //< Arrays holding actual surfactant density and its FFT
  DFT_FFT surfactant_potential_; //< Arrays holding surfactant assymetric potential
  double Asurf_;                     //< strength of assymetric (v^2 term) surfactant interaction
  double rho_surf_;              //< parameter rho_0
  double ax_;                    //< VDW parameter for symmetric part of surfactant interaction
  bool bFixedN_;                 ///< Fixed particle number flag
};

#endif // __LUTSKO_DFT__
