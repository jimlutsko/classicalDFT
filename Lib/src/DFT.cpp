#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

#include <gsl/gsl_integration.h>

using namespace std;

#ifdef USE_OMP
#include <omp.h>
#endif

#include "DFT.h"

double DFT::Mu(const vector<double> &x, int species) const
{
  double mu = log(x[species]);

  if(fmt_)
    mu += fmt_->BulkMuex(x, allSpecies_, species);
  
  for(auto &interaction: Interactions_)
    mu += interaction->Mu(x,species);
  return mu;
}

double DFT::Omega(const vector<double> &x) const
{
  double omega = Fhelmholtz(x);
  for(int i=0;i<allSpecies_.size();i++)
    if(fabs(x[i]) > SMALL_VALUE) omega -= x[i]*Mu(x,i);
  return omega;
}

double DFT::Fhelmholtz(const vector<double> &x) const
{
  double F = 0.0;
  for(auto &y: x)
    if(fabs(y) > SMALL_VALUE)
      F += y*log(y)-y;
  if(fmt_)
    F += fmt_->BulkFex(x, allSpecies_);
       
  for(auto &interaction: Interactions_)
    F += interaction->Fhelmholtz(x);  

  return F;  
}  


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


#pragma omp declare reduction(SummationPlus: Summation: omp_out += omp_in) 

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
	
	int chunk = Ntot/20;
	long i;
  
#pragma omp parallel for				\
  shared( chunk, species, dV)				\
  private(i)						\
  schedule(static,chunk)				\
  reduction(SummationPlus:F)
	for(i=0;i<Ntot;i++)
	  {
	    double d0 = density.getDensity(i);
	    if(d0 <= -SMALL_VALUE) d0 = SMALL_VALUE; // should never happen!
	    F += (d0*log(d0)-d0)*dV;
	    species->addToForce(i,log(d0)*dV);
	  }
      }
  
  for(auto &species : allSpecies_)
    F += species->externalField(true); // bCalcForces = true: obsolete?  

  if(fmt_)
    {    
      try{
	F += fmt_->calculateFreeEnergyAndDerivatives(allSpecies_);
      } catch( Eta_Too_Large_Exception &e) {
	throw e;
      }
    }
  //< Mean field contribution to F and dF
  // Need the following only if the fmt object is not called
  if(!fmt_)
    for(auto &species : allSpecies_)
      {
	species->doFFT();
      }
  
  for(auto &interaction: DFT::Interactions_)
    F += interaction->getInteractionEnergyAndForces();

  return F.sum();  
}
