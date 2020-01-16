#ifndef __LUTSKO_DFT__
#define __LUTSKO_DFT__

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <complex>

//TEST


#include "Enskog.h"

#include "Potential1.h"
#include "Density.h"
#include "FMT.h"

#include "Interaction.h"


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
   *   @param s: the first species. It makes no sense to create a DFT object without at least one species.
   */  
 DFT(Species *s) : fmt_(NULL) {allSpecies_.push_back(s);}
  /**
   *   @brief  Default  destructor for DFT
   *  
   *   @return nothing 
   */  
  ~DFT(){}

  /**
   *   @brief  Tells the DFT that there is another species
   *  
   *   @param  s is the Species object
   */ 
  void addSpecies(Species* s) {allSpecies_.push_back(s);}

  /**
   *   @brief  Tells the DFT that there is another interaction
   *  
   *   @param  I is the Interction object
   */ 
  void addInteraction(Interaction* Interaction) {Interactions_.push_back(Interaction);}

  /**
   *   @brief  Specify a hard-sphere free energy functional
   *  
   *   @param  fmt is the FMT object
   */ 
  void addHardCoreContribution(FMT *fmt) {fmt_ = fmt;}  

  /**
   *   @brief  Requests the number of species
   *  
   *   @returns the size of allSpecies_
   */ 
  int getNumberOfSpecies() const {return allSpecies_.size();}

  /**
   *   @brief  Requests the number of atoms of species s
   *  
   *   @param  species 
   *   @returns the Density-object's calculation of the mass
   */ 
  double getNumberAtoms(int species) const {return allSpecies_[species]->getDensity().getNumberAtoms();}

  /**
   *   @brief  Requests a read-only access to the Lattice object
   *  
   *   @returns a read-only reference to the Lattice object
   */ 
  const Lattice& lattice() const {return allSpecies_.front()->getLattice();}

  /**
   *   @brief  Requests the value of the convergence criterion for free energy minimization
   *  
   *   @returns the value of the convergence criterion for free energy minimization
   */   
  double get_convergence_monitor() const { double d = 0; for(auto& s: allSpecies_) d += s->get_convergence_monitor(); return d;}

  /**
   *   @brief  Requests the Density object associated with one of the species
   *  
   *   @param   species: the species
   *   @returns a read-only reference to the density object
   */     
  const Density& getDensity(int species) const {return allSpecies_[species]->getDensity();}

  
  DFT_Vec &getDF(int i) {return allSpecies_[i]->getDF();}
  void setDF(int i, DFT_Vec &df) {return allSpecies_[i]->setDF(df);}

  
  void doDisplay(string &title, string &file)
  {
    for(auto &x: allSpecies_) x->doDisplay(title,file);
  }

  /**
   *   @brief  Sets the densities based on amplitudes by passing to the species object
   *  
   *   @param  i is the species
   *   @param x is the amplitde
   */  
  void set_density_from_amplitude(int i,DFT_Vec &x) {allSpecies_[i]->set_density_from_amplitude(x);}
  void set_density(int i,DFT_Vec &x) {allSpecies_[i]->set_density(x);}
  void set_density(int i, long j, double x) {allSpecies_[i]->set_density(j,x);}

  void set_density_from_eta(int Jspecies) {((FMT_Species*) allSpecies_[Jspecies])->set_density_from_eta();}
  
  void writeDensity(int i, string &of) const {allSpecies_[i]->getDensity().writeDensity(of);}

  /**
   *   @brief  Calculates total grand canonical free energy and dOmega/dRho(i) for each lattice point using FFT convolutions. 
   *            This is just a wrapper that calls calculateFreeEnergyAndDerivatives_internal_ which does the real work. 
   *  
   *   @param  onlyFex : if true, the ideal and external contributions are not added: IS THIS PARAMETER NEEDED???
   *   @return total free energy of system
   */  
  double calculateFreeEnergyAndDerivatives(bool onlyFex);

  /**
   *   @brief  Compute chemical potential/kT for a uniform system with the given density
   *  
   *   @param  x is the array of densities
   *   @param  species is the species for which we calculate the chemical potential
   *   @return mu/kT
   */   
  virtual double Mu(const vector<double> &x, int species) const;

  /**
   *   @Brief  Compute grand potenial/kT/V for a uniform system with the given density
   *  
   *   @param  x is the array of densities
   *   @return Omega/(kT * V)
   */   
  virtual double Omega(const vector<double> &x) const;

  /**
   *   @brief  Compute Helmholtz free energy/kT/V for a uniform system with the given density
   *  
   *   @param  x is the array of densities
   *   @return F/(kT * V)
   */   
  virtual double Fhelmholtz(const vector<double> &x) const;

  /**
   *   @brief  Returns the name of the DFT object - i.e. identifies the model being used.
   *   @return the name as a string object.
   */     
  virtual string Name() const { return string("DFT_Ideal_Gas");}

  /**
   *   @brief  This conveys information to the FMT object for models that use it.
   */     
  virtual void setEtaMax(double etaMax) {}

  /**
   *   @brief  Only implemented for DFT_FMT object.
   */     
  virtual double Xliq_From_Mu(double mu) const {throw std::runtime_error("Not implemented");}

 protected:
  
  /**
   *   @brief  Calculates total grand canonical free energy and dOmega/dRho(i) for each lattice point using FFT convolutions. Here, this is just the ideal gas contribution.
   *  
   *   @param  onlyFex : if true, the ideal and external contributions are not added
   *   @return total free energy of system
   */  
  virtual double calculateFreeEnergyAndDerivatives_internal_(bool onlyFex);


 protected:
  vector<Species*> allSpecies_; ///< array holding the species objects
  vector<Interaction*> Interactions_; ///< array holding the interactions
  FMT *fmt_;
};


/**
  *  @brief DFT_VDW_Surfactant class
  *
  *   @detailed This extends DFT_VDW with functionality to model a surfactant. The additional contributions depend only on gradients
  *              of the other species, so the bulk thermodynamics is just the same as that inherited from DFT_VDW. 
  *  
  */  
/*
template <class T> class DFT_VDW_Surfactant : public DFT_VDW<T>
{
 public:
  **
  *   @brief  Default  constructor for DFT 
  *  
  *   @param s: the first species. It makes no sense to create a DFT object without at least one species. We assume that this is the "water".
  *   @param pot: the surfactant potential.
  *   @param kT: the temperature
  *   @return nothing 
  *  
  DFT_VDW_Surfactant(Species *species, Potential1 &pot, double kT);

  **
  *   @brief  Default  destructor
  *  
  ~DFT_VDW_Surfactant(){}

  **
   *   @brief  The name of this DFT object
   *  
   *   @return name as string
   *   
  virtual string Name() const { return string("DFT_VDW_Surfactant : ") +  DFT_VDW<T>::Name();} 

 protected:
  **
   *   @brief  Calculates total grand canonical free energy and dOmega/dRho(i) for each lattice point using FFT convolutions
   *  
   *   @param  density is current density
   *   @param  mu is the chemical potential
   *   @return total free energy of system
   *  
  virtual double calculateFreeEnergyAndDerivatives_internal_(bool onlyFex = false);
  
 protected:
  DFT_FFT surfactant_potential_; //< Arrays holding surfactant assymetric potential
};
*/
#endif // __LUTSKO_DFT__
