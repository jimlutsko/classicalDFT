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

double DFT::F_IdealGas(bool bCalcForces)
{
  double Fideal = 0.0;
  
  for(auto &species : allSpecies_)
    {
      const Density& density = species->getDensity();
      double dV = density.dV();
      
      for(long i=0;i<density.Ntot();i++)
	{
	  double d0 = density.getDensity(i);
	  if(d0 > -SMALL_VALUE)
	    {
	      Fideal += (d0*log(d0)-d0)*dV;
	      if(bCalcForces) species->addToForce(i,log(d0)*dV);
	    } else {
	    if(bCalcForces) species->addToForce(i,log(SMALL_VALUE)*dV);
	  }
	}
    }
  return Fideal;
}

double DFT::F_External(bool bCalcForces)
{
  double Fx = 0.0;
  
  for(auto &species : allSpecies_)
    Fx += species->externalField(bCalcForces);

  return Fx;
}


template <class T>
double DFT_FMT<T>::calculateFreeEnergyAndDerivatives(bool onlyFex)
{
  double F = 0;
  
  //dF.zeros(density.Ntot());

  try{
    F += fmt_.calculateFreeEnergyAndDerivatives(allSpecies_);
  } catch( Eta_Too_Large_Exception &e) {
    throw e;
  }
  if(!onlyFex) // add in ideal gas and external field and chem potential
    {
      F += F_IdealGas(true);       // Ideal gas 
      F += F_External(true);   // External field + chem potential
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
DFT_VDW<T>::DFT_VDW(int Nx, int Ny, int Nz)
  : DFT_FMT<T>(Nx,Ny,Nz), vdw_(0,0)
{
  //  dft_fmt_ = new DFT_FMT<T>(Nx,Ny,Nz);
  // initialize working space for calculations ...
  v_mean_field_.initialize(Nx,Ny,Nz);
}

template <class T>
void DFT_VDW<T>::addSpecies(VDW_Species* species,  double kT)  
{
  DFT::addSpecies(species);
  
  double hsd   = species->getHSD();
  double a_vdw = species->get_VDW_Constant();
  
  vdw_.set_VDW_Parameter(a_vdw);
  vdw_.set_HardSphere_Diameter(hsd);
}

template <class T>
double DFT_VDW<T>::calculateFreeEnergyAndDerivatives(bool onlyFex)
{
  //  dF.zeros(Ntot);
  double F = 0;

  // Hard sphere contributions to F and dF
  try {
    F = DFT_FMT<T>::calculateFreeEnergyAndDerivatives(onlyFex);
  } catch( Eta_Too_Large_Exception &e) {
    throw e;
  }

  // Mean field contribution to F and dF
  // Divide by Ntot because of normalization of fft
  for(auto &s: DFT::allSpecies_)
    F += ((VDW_Species *) s)->getInteractionEnergyAndForces(v_mean_field_);

  return F;
}



// These lines are needed to force instantiation of the relevant classes. 


template class DFT_FMT<WhiteBearI>;
template class DFT_FMT<WhiteBearII>;
template class DFT_FMT<RSLT>;
template class DFT_FMT<RSLT2>;

template class DFT_VDW<WhiteBearI>;
template class DFT_VDW<WhiteBearII>;
template class DFT_VDW<RSLT>;
template class DFT_VDW<RSLT2>;
