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


double DFT_IdealGas::calculateFreeEnergyAndDerivatives(bool includeIdealGas)
{
  double F = 0;
  for(auto &s: allSpecies_)
    s->zeroForce();

  if(includeIdealGas) F += DFT::F_IdealGas();
  F += F_External();

  return F;
}

