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


FMT_Species_EOS::FMT_Species_EOS(double D_EOS, EOS &eos, Density& density, double hsd, double mu, int seq)
  : FMT_Species(density,hsd,mu,seq), eos_weighted_density_(1), D_EOS_(D_EOS), eos_(eos)
										    
{
  long Nx = density_->Nx();
  long Ny = density_->Ny();
  long Nz = density_->Nz();

  eos_weighted_density_[0].initialize(Nx, Ny, Nz);
  generateWeights(D_EOS*hsd, eos_weighted_density_);
  eos_weighted_density_[0].transformWeights();  
}

void FMT_Species_EOS::calculateFundamentalMeasures(bool needsTensor)
{
  FMT_Species::calculateFundamentalMeasures(needsTensor);  
  const DFT_Vec_Complex &rho_k = density_->get_density_fourier();    
  eos_weighted_density_[0].convoluteWith(rho_k);          
}

void FMT_Species_EOS::calculateForce(bool needsTensor, void *param)
{
  FMT_Species::calculateForce(needsTensor);

  double dV = getLattice().dV();
  int Nx    = getLattice().Nx();
  int Ny    = getLattice().Ny();
  int Nz    = getLattice().Nz();
  
  DFT_FFT dPhi(Nx,Ny,Nz);
  dPhi.Four().zeros();
  
  eos_weighted_density_[0].add_weight_schur_dPhi_to_arg(dPhi.Four());  

  dPhi.do_fourier_2_real();
  dPhi.Real().MultBy(dV);
  addToForce(dPhi.cReal());   
}

