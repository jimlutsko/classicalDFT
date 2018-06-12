#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

using namespace std;

#ifdef USE_OMP
#include <omp.h>
#endif

#include "DFT.h"


double DFT_IdealGas::calculateFreeEnergyAndDerivatives(Density& density, double mu, DFT_Vec& dF, bool includeIdealGas)
{
  double dV = density.dV();

  double F = 0;
  dF.zeros(density.Ntot());

  // Hard sphere contributions to F and dF
  F = 0;

  // Ideal gas contribution to F and dF
  if(includeIdealGas)
    for(long i=0;i<density.Ntot();i++)
      {
	double d0 = density.getDensity(i);
	if(d0 > -SMALL_VALUE)
	  {
	    F += (d0*log(d0)-d0)*dV;
	    dF.addTo(i,log(d0)*dV);
	  } else {
	  dF.addTo(i,log(SMALL_VALUE)*dV);
	}
      }
  // External field
  F  += density.getExternalFieldEnergy()*dV - density.getNumberAtoms()*mu;
  dF.Increment_Shift_And_Scale(density.getExternalField(),dV,mu);
  return F;
}

double DFT_IdealGas::getDF_DRHO(Density &density, double mu, DFT_Vec& dF_dRho)
{
  double F = 0;

  try {
    F = calculateFreeEnergyAndDerivatives(density, mu,dF_dRho,true);
  } catch( Eta_Too_Large_Exception &e) {
     throw e;
  }

  dF_dRho.multBy(1.0/density.dV());

  return F;
}
