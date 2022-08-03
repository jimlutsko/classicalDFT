#ifndef __LUTSKO_DFT__
#define __LUTSKO_DFT__

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <complex>

#include <boost/serialization/vector.hpp>



#include "Enskog.h"

#include "Potential1.h"
#include "Density.h"
#include "FMT.h"

#include "Interaction.h"
#include "Dynamical_Matrix.h"

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

class DFT : public Dynamical_Matrix
{
 public:
  // Xtors
  DFT(Species *s = NULL) : fmt_(NULL) {if(s) addSpecies(s);}
  ~DFT(){}

  // Add pieces
  void addSpecies(Species* s) {if(s == NULL) return; allSpecies_.push_back(s); s->setIndex(allSpecies_.size()-1);}
  void addInteraction(Interaction_Base* Interaction) {Interactions_.push_back(Interaction);}
  void addHardCoreContribution(FMT *fmt) {fmt_ = fmt;}  

  // Accessors
  const Lattice& get_lattice() const {return allSpecies_.front()->getLattice();}
  const Density& getDensity(int species) const {return allSpecies_[species]->getDensity();}
  Species *getSpecies(int s) const { return allSpecies_[s];}
  
  int    getNumberOfSpecies() const {return allSpecies_.size();}
  double getNumberAtoms(int species) const {return allSpecies_[species]->getDensity().getNumberAtoms();}  
  double get_convergence_monitor() const { double d = 0; for(auto& s: allSpecies_) d += s->get_convergence_monitor(); return d;}
  
  DFT_Vec &getDF(int i) {return allSpecies_[i]->getDF();}

  // Set
  void setDF(int i, DFT_Vec &df) {return allSpecies_[i]->setDF(df);}
  void set_density(int i,DFT_Vec &x) {allSpecies_[i]->set_density(x);}
  void set_density(int i, long j, double x) {allSpecies_[i]->set_density(j,x);}
  void set_densities_from_aliases(vector<DFT_Vec> &x_);
  void convert_dF_to_alias_derivs(vector<DFT_Vec> &x_);
  
  // A few actions  
  void doDisplay(string &title, string &file, void *param = NULL) { for(auto &x: allSpecies_) x->doDisplay(title,file, param);}
  void writeDensity(int i, string &of) const {allSpecies_[i]->getDensity().writeDensity(of);}
  double calculateFreeEnergyAndDerivatives(bool onlyFex);

  // Bulk Thermodynamics

  virtual double mu_times_beta(const vector<double> &densities, int species) const;
  virtual double omega_times_beta_over_volume(const vector<double> &x) const;
  virtual double fhelmholtz_times_beta_over_volume(const vector<double> &x) const;

  virtual double mu_times_beta(double densities) const;
  virtual double omega_times_beta_over_volume(double density) const;
  virtual double fhelmholtz_times_beta_over_volume(double density) const;  

  void   findSpinodal(double xmax, double dx, double &xs1, double &xs2, double tol) const;
  double find_density_from_mu(double mu, double xmin, double xmax, double tol) const;
  void findCoex(double xmax, double dx, double &x1, double &x2, double tol) const;
  virtual void getCriticalPoint(Potential1& p, double &xc, double &Tc, double HSD = -1) const;

  // Bulk structural properties
  double real_space_dcf(double r, double x) const;
  double fourier_space_dcf(double k, double x) const;
  
  // Identifiers
  virtual string Name() const { return string("DFT_Ideal_Gas");}
  string get_fmt_name() const
  {
    if(fmt_ == NULL) return string("none");
    return fmt_->get_name();
  }
  string get_potential_name(int which = 0) const
  {
    if(Interactions_.size() < which+1) return string("none");
    return Interactions_[which]->get_name();
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


  // Implement Dynamical_Matrix interface.
  // Second derivatives contracted into arbitrary vector
  virtual void     matrix_dot_v(const vector<DFT_FFT> &v, vector<DFT_Vec> &result, void *param = NULL) const;

  virtual unsigned get_dimension(int direction) const {return allSpecies_[0]->getLattice().get_dimension(direction);}
  virtual long     get_Nboundary()              const {return allSpecies_[0]->getDensity().get_Nboundary();}
  virtual long     boundary_pos_2_pos(int p)    const {return allSpecies_[0]->getDensity().boundary_pos_2_pos(p);}
  virtual bool     get_next_boundary_point(int &ix, int &iy, int &iz) const {return allSpecies_[0]->getLattice().get_next_boundary_point(ix,iy,iz);}
  virtual bool     get_next_boundary_point(long &pos) const {return allSpecies_[0]->getLattice().get_next_boundary_point(pos);}
  virtual bool     is_fixed_boundary() const { return allSpecies_[0]->is_fixed_boundary();}
  virtual bool is_boundary_point(long p) const {return allSpecies_[0]->getDensity().is_boundary_point(p);}
  
  void set_full_hessian(bool full) { full_hessian_ = full;}
  
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & allSpecies_;
    ar & Interactions_;
    ar & fmt_;
    ar & F_id_;
    ar & F_ext_;
    ar & F_hs_;
    ar & F_mf_;	
  }

 protected:
  vector<Species*> allSpecies_; ///< array holding the species objects
  vector<Interaction_Base*> Interactions_; ///< array holding the interactions
  FMT *fmt_;
  double F_id_  = 0.0; ///< Ideal part of free energy
  double F_ext_ = 0.0; ///< External field contribution to free energy (including chemical potential)
  double F_hs_  = 0.0; ///< Hard-sphere contribution to free energy
  double F_mf_  = 0.0; ///< Mean-field contribution to free energy

  mutable bool full_hessian_ = true; // used in matrix_holder to distinguish full hessian from excess hessian.
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
