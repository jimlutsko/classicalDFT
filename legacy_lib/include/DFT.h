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
  DFT(Species *s = NULL) : fmt_(NULL) {if(s) allSpecies_.push_back(s);}
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
  void addInteraction(Interaction_Base* Interaction) {Interactions_.push_back(Interaction);}

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

  //  void set_density_from_eta(int Jspecies) {((FMT_Species*) allSpecies_[Jspecies])->set_density_from_eta();}
  
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
   *   @brief  Only implemented for single species and one interaction
   */     
  virtual double XLiq_From_Mu(double mu, double high_density) const;

  /**
   *   @brief  Only implemented for single species and one interaction
   */     
  virtual double XVap_From_Mu(double mu, double high_density) const;

  /**
   *   @brief  Determines xliq (largest value less than close packing) from beta P. Only implemented for single species and one interaction
   */     
  virtual double XLiq_from_P(double P) const;
  
  /**
   *   @brief  Only implemented for single species and one interaction
   */     
  virtual void spinodal(double &xs1, double &xs2) const;

  /**
   *   @brief  Only implemented for single species and one interaction
   */     
  virtual void liq_vap_coex(double &xs1, double &xs2, double &x1, double &x2) const;

  /**
   *   @brief  Only implemented for single species and one interaction
   */     
  virtual void getCriticalPoint(Potential1& p, double &xc, double &Tc)
  {
    //    if(Interactions_.size() != 1 || allSpecies_.size() != 1) throw std::runtime_error("DFT::getCriticalPoint only implemented for exactly 1 species and 1 interaction");
    double T1 = Tc;
    
    do{
      T1 = Tc;
      double hsd = p.getHSD(Tc);
      xc = 0.13044*(6.0/M_PI)*pow(hsd,-3.0);
      Tc = -0.090082*2*p.getVDW_Parameter(Tc)*pow(hsd,-3.0)*Tc;
    } while(fabs(T1-Tc) > 1e-8*(T1+Tc));
  }

  // protected:
  
  /**
   *   @brief  Calculates total grand canonical free energy and dOmega/dRho(i) for each lattice point using FFT convolutions. Here, this is just the ideal gas contribution.
   *  
   *   @param  onlyFex : if true, the ideal and external contributions are not added
   *   @return total free energy of system
   */  
  virtual double calculateFreeEnergyAndDerivatives_internal_(bool onlyFex);

  /**
   *   @brief  Returns ideal gas contribution to free energy
   *  
   *   @return Ideal part of free energy
   */  
  double get_f_id() const { return F_id_;}

    /**
   *   @brief  Returns external field contribution to free energy including chemical potential
   *  
   *   @return External field contribution to free energy (including chemical potential)
   */  
  double get_f_ext() const { return F_ext_;}

    /**
   *   @brief  Returns hard-sphere contribution to free energy
   *  
   *   @return Hard-sphere contribution to free energy
   */  
  double get_f_hs() const { return F_hs_;}

  /**
   *   @brief  Returns mean field contribution to free energy
   *  
   *   @return Mean-field contribution to free energy
   */  
  double get_f_mf() const { return F_mf_;}
  

 protected:
  vector<Species*> allSpecies_; ///< array holding the species objects
  vector<Interaction_Base*> Interactions_; ///< array holding the interactions
  FMT *fmt_;
  double F_id_  = 0.0; ///< Ideal part of free energy
  double F_ext_ = 0.0; ///< External field contribution to free energy (including chemical potential)
  double F_hs_  = 0.0; ///< Hard-sphere contribution to free energy
  double F_mf_  = 0.0; ///< Mean-field contribution to free energy
};


namespace dft_util {
/**
 *   @brief  This will fit the value of thethe vDW coefficient to the supplied data for a single species. The input densities are supplied in dimensionless form, 
 *           so they imply some length scale  l. The output is the hsd/l and the vdw parmeter in the form  beta a/l^3. Normally, for a potential with length scale sigma and energy scale epsilon, the
 *           calculated vdw parameter is some number times beta epsilon sigma^3 so one could match the two either by adjusting epsilon or sigma or both. 
 *  
 *   @param  data is the array of pairs of densities and z_factor = P/(n kT) for which the pressure is being supplied
 *   @param hsd: this is the hard-sphere diameter.
 *   @returns the vDW parameter in the dimensionless form beta avdw/l^3. 
 */   
  double fit_avdw_to_data(vector<pair<double, double> > &data, double hsd);

/**
 *   @brief  This will fit the value of the hsd and the vDW coefficient to the supplied data for a single species. The input densities are supplied in dimensionless form, 
 *           so they imply some length scale  l. The output is the hsd/l and the vdw parmeter in the form  beta a/l^3. Normally, for a potential with length scale sigma and energy scale epsilon, the
 *           calculated vdw parmaeter is some number times beta epsilon sigma^3 so one could match the two either by adjusting epsilon or sigma or both. 
 *  
 *   @param  data is the array of pairs of densities and z_factor = P/(n kT) for which the pressure is being supplied
 *   @param hsd: this is an in/out parameter. The supplied value is used as an initial guess in the numerical search. Afterwards, it contains the determined value.
 *   @param aVDW this is an output parameter. It is actually  beta aVDW/l^3 . 
 *   @param tol is the tolerence of the fit. 
 */   

  double fit_to_data(std::vector<std::pair<double, double> > &data, double &avdw, double tol = 1e-8);
}  


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
