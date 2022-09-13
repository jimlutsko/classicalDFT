#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>


using namespace std;

#ifdef USE_MPI
#include <mpi.h>  
#endif

#include "Species.h"

FMT_Species_Extended::FMT_Species_Extended(Density& density, double hsd, double hsd2, double mu, int seq)
  : FMT_Species(density,hsd,mu,seq), psi_(density_.Nx(), density_.Ny(), density_.Nz()), vsurf_(3)
										    
{
  long Nx = density_.Nx();
  long Ny = density_.Ny();
  long Nz = density_.Nz();

  for(FMT_Weighted_Density &d: vsurf_)
    d.initialize(Nx, Ny, Nz);

  generateWeights(hsd2, NULL, NULL, vsurf_.data(), NULL);
  
  for(FMT_Weighted_Density &d: vsurf_)    
    d.transformWeights();  
}

