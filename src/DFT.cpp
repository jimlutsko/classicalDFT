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


double DFT::calculateFreeEnergyAndDerivatives_internal_(bool onlyFex)
{
  double F = 0.0;
  
  if(onlyFex) return F; // do nothing.

  for(auto &species : allSpecies_)
    {
      const Density& density = species->getDensity();
      double dV = density.dV();
      
      for(long i=0;i<density.Ntot();i++)
	{
	  double d0 = density.getDensity(i);
	  if(d0 > -SMALL_VALUE)
	    {
	      F += (d0*log(d0)-d0)*dV;
	      species->addToForce(i,log(d0)*dV);
	    } else {
	    species->addToForce(i,log(SMALL_VALUE)*dV);
	  }
	}
    }

  for(auto &species : allSpecies_)
    F += species->externalField(true); // bCalcForces = true: obsolete?  

  return F;
  
}

template <class T>
double DFT_FMT<T>::calculateFreeEnergyAndDerivatives_internal_(bool onlyFex)
{
  double F = DFT::calculateFreeEnergyAndDerivatives_internal_(onlyFex);
  
  try{
    F += fmt_.calculateFreeEnergyAndDerivatives(allSpecies_);
  } catch( Eta_Too_Large_Exception &e) {
    throw e;
  }

  return F;
}


template <class T>
double DFT_FMT<T>::Xliq_From_Mu(double mu) const 
{ 
  double xa = 0.5;
  double xb = 0.5;

  while(Enskog(xa).chemPotentialCS() > mu) { xa /= 2;}

  double xmax = 6/M_PI;

  while(Enskog(xb).chemPotentialCS() < mu) { xb = (xmax+xb)/2; }

  while(fabs(xa-xb) > 1e-8)
    {
      double x = (xa+xb)/2;
      if(Enskog(x).chemPotentialCS() > mu) xb = x;
      else xa = x;
    }
  return (xa+xb)/2;
}

#include "Potential1.h"

template <class T>
DFT_VDW<T>::DFT_VDW(Species *s) : DFT_FMT<T>(s) {}

template <class T>
double DFT_VDW<T>::Mu(const vector<double> &x, int species) const
{
  VDW_Species *s = (VDW_Species*) DFT_FMT<T>::allSpecies_[species];
  
  return DFT_FMT<T>::Mu(x,species) + 2*s->get_VDW_Constant()*x[species];
}

template <class T>
double DFT_VDW<T>::Fhelmholtz(const vector<double> &x) const
{
  double f = DFT_FMT<T>::Fhelmholtz(x);
  for(int s=0; s < DFT_FMT<T>::allSpecies_.size(); s++)
    {
      VDW_Species *sp = (VDW_Species*) (DFT_FMT<T>::allSpecies_[s]);
      f += sp->get_VDW_Constant()*x[s]*x[s];
    }
  return f;
}

template <class T>
double DFT_VDW<T>::calculateFreeEnergyAndDerivatives_internal_(bool onlyFex)
{
  double F = 0;

  // Hard sphere contributions to F and dF
  try {
    F = DFT_FMT<T>::calculateFreeEnergyAndDerivatives_internal_(onlyFex);
  } catch( Eta_Too_Large_Exception &e) {
    throw e;
  }
  // Mean field contribution to F and dF
  for(auto &s: DFT::allSpecies_)
      F += ((VDW_Species *) s)->getInteractionEnergyAndForces();

  // Mean field contribution to F and dF
  //  for(auto &I: DFT::Interactions_)
  //      F += I->getInteractionEnergyAndForces();  

  return F;
}



// These lines are needed to force instantiation of the relevant classes. 


template class DFT_FMT<WhiteBearI>;
template class DFT_FMT<WhiteBearII>;
template class DFT_FMT<RSLT>;
template class DFT_FMT<RSLT2>;
template class DFT_FMT<Rosenfeld>;

template class DFT_VDW<WhiteBearI>;
template class DFT_VDW<WhiteBearII>;
template class DFT_VDW<RSLT>;
template class DFT_VDW<RSLT2>;
template class DFT_VDW<Rosenfeld>;
