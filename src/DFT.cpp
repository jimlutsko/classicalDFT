#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

using namespace std;

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_gamma.h>

#ifdef USE_OMP
#include <omp.h>
#endif

#include "DFT.h"


// note that here, r = |r1-r2|
// NEEDS UPDATE FOR FMT_SPECIES_EOS
double DFT::real_space_dcf(double r, double x) const
{
  if(allSpecies_.size() != 1)
    throw std:: runtime_error("DFT::real_space_dcf only implemented for single species");
  
  double dcf = 0;
  double hsd = allSpecies_[0]->getHSD();

  if(fmt_) dcf += fmt_->get_real_space_dcf(r,x,hsd);

  for(auto &x: Interactions_)
    dcf -= x->getW(r);
    
  return dcf;
}

// NEEDS UPDATE FOR FMT_SPECIES_EOS
double DFT::fourier_space_dcf(double k, double x) const
{
  if(allSpecies_.size() != 1)
    throw std:: runtime_error("DFT::real_space_dcf only implemented for single species");
  
  double dcf = 0;
  double hsd = allSpecies_[0]->getHSD();

  if(fmt_) dcf += fmt_->get_fourier_space_dcf(k,x,hsd);
  
  return dcf;
}

// NEEDS UPDATE FOR FMT_SPECIES_EOS
double DFT::mu_times_beta(double density) const
{
  return mu_times_beta(vector<double>(1,density),0);
}

// NEEDS UPDATE FOR FMT_SPECIES_EOS
double DFT::mu_times_beta(const vector<double> &x, int species) const
{
  Summation mu;

  mu += log(x[species]);

  if(fmt_)
    mu += fmt_->BulkMuex(x, allSpecies_, species);
  
  for(auto &interaction: Interactions_)
    mu += interaction->Mu(x,species);
  
  return mu.sum();
}

// NEEDS UPDATE FOR FMT_SPECIES_EOS
double DFT::omega_times_beta_over_volume(double density) const
{
  return omega_times_beta_over_volume(vector<double>(1,density));
}

// NEEDS UPDATE FOR FMT_SPECIES_EOS
double DFT::omega_times_beta_over_volume(const vector<double> &x) const
{
  double omega = fhelmholtz_times_beta_over_volume(x);
  for(int i=0;i<allSpecies_.size();i++)
    omega -= x[i]*mu_times_beta(x,i);

  return omega;
}

// NEEDS UPDATE FOR FMT_SPECIES_EOS
double DFT::fhelmholtz_times_beta_over_volume(double density) const
{
  return fhelmholtz_times_beta_over_volume(vector<double>(1,density));
}

// NEEDS UPDATE FOR FMT_SPECIES_EOS
double DFT::fhelmholtz_times_beta_over_volume(const vector<double> &x) const
{
  double F = 0.0;

  double V = get_lattice().getVolume();

  for(auto &y: x)
    F += y*log(y)-y;

  double Fhs = 0.0;
  if(fmt_)
    {
      Fhs += fmt_->BulkFex(x, allSpecies_);
      F += Fhs;
    }

  double Fmf = 0.0;
  for(auto &interaction: Interactions_)
    Fmf += interaction->Fhelmholtz(x);

  F += Fmf;

  return F;  
}  

// NEEDS UPDATE FOR FMT_SPECIES_EOS
void DFT::set_densities_from_aliases(vector<DFT_Vec> &x_)
{
  for(int s=0; s<allSpecies_.size();s++)
    allSpecies_[s]->set_density_from_alias(x_[s]);
}

// NEEDS UPDATE FOR FMT_SPECIES_EOS
void DFT::convert_dF_to_alias_derivs(vector<DFT_Vec> &x_)
{
  for(int s = 0; s<allSpecies_.size(); s++)
    allSpecies_[s]->convert_to_alias_deriv(x_[s],getDF(s));
}

// NEEDS UPDATE FOR FMT_SPECIES_EOS
void DFT::convert_dF_to_alias_derivs()
{
  for(int s = 0; s<allSpecies_.size(); s++)
    allSpecies_[s]->convert_to_alias_deriv(getDF(s));
}

// NEEDS UPDATE FOR FMT_SPECIES_EOS
double DFT::calculateFreeEnergyAndDerivatives(bool onlyFex)
{
  for(auto &species : allSpecies_)  
    species->zeroForce();
  
  for(auto &s: allSpecies_) s->beginForceCalculation();

  double F = calculateFreeEnergyAndDerivatives_internal_(onlyFex);

  for(auto &s: allSpecies_)
    s->endForceCalculation();
  
  return F;
}

#ifdef USE_OMP    
#pragma omp declare reduction(SummationPlus: Summation: omp_out += omp_in) 
#endif
// NEEDS UPDATE FOR FMT_SPECIES_EOS
double DFT::calculateFreeEnergyAndDerivatives_internal_(bool onlyFex)
{
  Summation F;
  // Ideal gas contribution  
  if(!onlyFex) 
    for(auto &species : allSpecies_)
      {
	const Density& density = species->getDensity();
	double dV = density.dV();
	long Ntot = density.Ntot();
	long pos;
#ifdef USE_OMP    
#pragma omp parallel for shared(species, dV) private(pos) schedule(static) reduction(SummationPlus:F)
#endif
	for(pos=0;pos<Ntot;pos++)
	  {
	    double d0 = density.get(pos);
	    F += (d0*log(d0)-d0)*dV;
	    species->addToForce(pos,log(d0)*dV); //HERE
	  }
      }
  F_id_ = F;

  // Hard-sphere contribution
  if(fmt_)
    {    
      try{
	F_hs_ = fmt_->calculateFreeEnergyAndDerivatives(allSpecies_);
	F += F_hs_;
      } catch( Eta_Too_Large_Exception &e) {
	throw e;
      }
    } else {   // Need the following only if the fmt object is not called
    for(auto &species : allSpecies_)
      species->doFFT();
  }

  //< Mean field contribution to F and dF
  F_mf_ = 0;
  for(auto &interaction: DFT::Interactions_)    
    F_mf_ += interaction->getInteractionEnergyAndForces();
  F += F_mf_;

  // External field + chemical potential
  F_ext_ = 0;


  // add in the chemical potential ...
  
  for(auto &species: allSpecies_)
    F_ext_ += species->evaluate_contribution_chemical_potential();

  
  for(auto &field : external_fields_)
    F_ext_ += allSpecies_[field->get_species()]->evaluate_external_field(*field);
  
  //  for(auto &species : allSpecies_)
  //    F_ext_ += species->externalField(true); // bCalcForces = true: obsolete?
  F += F_ext_;

  if(offset_)
    {
      double Foff = 0;
      for(auto &species : allSpecies_)
	{
	  const Density& density = species->getDensity();
	  double dV = density.dV();
	  long Ntot = density.Ntot();
	  long pos;
#ifdef USE_OMP    
#pragma omp parallel for shared(species, dV) private(pos) schedule(static) reduction(SummationPlus:F)
#endif
	  for(pos=0;pos<Ntot;pos++)
	    {
	      double d0 = density.get(pos);
	      Foff -= log(d0);
	      species->addToForce(pos,-1.0/d0); //HERE
	    }
	}
      F += Foff;      

    }

  return F.sum();  
}

// computes (d2F/dn_i dn_j) v_j: 
//
//  Uses a cheap fix for fixed boundaries: set v_j=0 for j on boundary and F_{ij}v_j=0 for i on boundary

// NEEDS UPDATE FOR FMT_SPECIES_EOS
void DFT::matrix_dot_v_intern(const vector<DFT_FFT> &v, vector<DFT_Vec> &result, void *param, bool only_d2F) const
{
  // I would like to do this but it violates the const declaration of v
  //  for(int i=0;i<v.size();i++)
  //    v[i].do_real_2_fourier(); // make sure this is done! An internal flag should prevent needless FFT's
  
  // Boundary terms must be zero if the boundary is fixed
  for(int s=0;s<allSpecies_.size();s++)  
    if(allSpecies_[s]->is_fixed_boundary())
      for(unsigned p=0;p<allSpecies_[s]->getLattice().get_Nboundary();p++)
      {
        unsigned pp = allSpecies_[s]->getLattice().boundary_pos_2_pos(p);
        if(fabs(v[s].cReal().get(pp)) > 0.0)
          throw std::runtime_error("Input vector v must have zero boundary entries in DFT::hessian_dot_v when the species has fixed boundaries");
      }

  // ideal gas contribution: v_i/n_i  
  if(full_hessian_)
  {
    double dV = allSpecies_[0]->getDensity().dV();

    for(int s=0;s<allSpecies_.size();s++)
      #ifdef USE_OMP
      #pragma omp parallel for
      #endif
      for(unsigned pos=0;pos<v[s].cReal().size();pos++)
        result[s].set(pos, dV*v[s].cReal().get(pos)/allSpecies_[s]->getDensity().get(pos));
  }

  // Hard-sphere
  if(fmt_)
  {
    try {fmt_->add_second_derivative(v,result, allSpecies_);}
    catch( Eta_Too_Large_Exception &e) {throw e;}
  } else for(auto &species : allSpecies_) // needed for interactions - its done automatically if there is an fmt evaluation
	   species->doFFT(); 
  
  // Mean field
  for(auto &interaction: DFT::Interactions_)    
    interaction->add_second_derivative(v,result);

  // Remove boundary terms if the boundary is fixed
  for(int s=0;s<allSpecies_.size();s++)  
    if(allSpecies_[s]->is_fixed_boundary())
    {
      long p = 0;
      do{result[s].set(p,0.0);} while(get_next_boundary_point(p));
    }
}

// computes result_I = F_{I I+J} for fixed J
// NEEDS UPDATE FOR FMT_SPECIES_EOS
void DFT::diagonal_matrix_elements(int jx, int jy, int jz, vector<DFT_Vec> &result) const
{
  // ideal gas contribution: delta_J0/n_i  

  int Nx = allSpecies_[0]->getDensity().Nx();
  int Ny = allSpecies_[0]->getDensity().Ny();
  int Nz = allSpecies_[0]->getDensity().Nz();
  long Ntot = allSpecies_[0]->getDensity().Ntot();
  long J = allSpecies_[0]->getDensity().get_PBC_Pos(jx,jy,jz);
  double dV = allSpecies_[0]->getDensity().dV();

  if(full_hessian_)
    if((jx%Nx == 0) && (jy%Ny == 0) && (jz%Nz == 0))
      {
	for(int s=0;s<allSpecies_.size();s++)
#ifdef USE_OMP
#pragma omp parallel for
#endif
	  for(unsigned pos=0;pos<Ntot;pos++)
	    result[s].set(pos, dV/allSpecies_[s]->get_density(pos));
      }

  for(auto &species : allSpecies_)
    species->doFFT(); 
  
  // Hard-sphere
  if(fmt_)
  {    
    try {fmt_->add_second_derivative(jx,jy,jz, allSpecies_,result);}
    catch( Eta_Too_Large_Exception &e) {throw e;}
  } else for(auto &species : allSpecies_)
	   species->doFFT(); 
  
  // Mean field: just shift all entries by w(0)*dV*dV
  // N.B. w already contains one factor of dV.
  for(auto &interaction: DFT::Interactions_)
    if(interaction->get_s1() == interaction->get_s2())
      result[0].add(interaction->getW(J)*dV);

  // Remove boundary terms if the boundary is fixed
  
  for(int s=0;s<allSpecies_.size();s++)  
    if(allSpecies_[s]->is_fixed_boundary())
    {
      long p = 0;
      do{result[s].set(p,0.0);} while(get_next_boundary_point(p));
    }
  
}


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
double fit_avdw_to_data(vector<pair<double, double> > &data, double hsd)
{
  double s1 = 0.0;
  double s2 = 0.0;

  for(auto &p: data)
    {
      double e = p.first*M_PI*hsd*hsd*hsd/6;
      double Z = (1+e*(1+e*(1-e)))/pow(1-e,3);
      s1 += (p.second - Z)*p.first;
      s2 += p.first*p.first;
    }
  return 2*s1/s2;
}

static double objective(vector<pair<double, double> > &data, double hsd)
{
  double sum = 0.0;
  double a   = fit_avdw_to_data(data,hsd);

  for(auto &p: data)
    {      
      double e = p.first*M_PI*hsd*hsd*hsd/6;
      double Z = (1+e*(1+e*(1-e)))/pow(1-e,3);
      Z += 0.5*a*p.first;
      Z -= p.second;

      sum += Z*p.first*(2+2*e-e*e)*pow(1-e,-4);
    }
  return sum; 
}

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

double fit_to_data(std::vector<std::pair<double, double> > &data, double &avdw, double tol)
{
  double a = 0;
  double b = 0;
  double fa = 0;
  double fb = 0;

  for(double d = 0.5; d < 2.0; d += 0.1)
    {
      a = b;
      b = d;
      fa = fb;
      fb = objective(data, d);
      if(fa*fb < 0) break;
    }

  if(fa*fb > 0) throw std::runtime_error("No crossing found in fit_to_data");

  while(fabs(a-b) > tol*(a+b))
    {
      double c = (a+b)/2;
      double fc = objective(data,c);
      if(fa*fc > 0) {a = c; fa = fc;}
      else { b = c; fb = fc;}
    }

  avdw = fit_avdw_to_data(data,(a+b)/2);  
  return (a+b)/2;
}
}  
