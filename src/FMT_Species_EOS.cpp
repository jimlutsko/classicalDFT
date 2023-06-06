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
#include "FMT.h"

FMT_Species_EOS::FMT_Species_EOS(double D_EOS, EOS &eos, double avdw, Density& density, double hsd, double mu, int seq)
  : FMT_Species(density,hsd,mu,seq), eos_weighted_density_(1), D_EOS_(D_EOS), eos_(eos), avdw_(avdw)
										    
{
  long Nx = density_->Nx();
  long Ny = density_->Ny();
  long Nz = density_->Nz();

  eos_weighted_density_[0].initialize(Nx, Ny, Nz);
  generateWeights(D_EOS*hsd, eos_weighted_density_);
  eos_weighted_density_[0].transformWeights();  
}

double FMT_Species_EOS::effDensity(long I)
{
  double eta = eos_weighted_density_[0].real(I);
  double hsd3 = hsd_*hsd_*hsd_*D_EOS_*D_EOS_*D_EOS_;
  double x   = 6*eta/(M_PI*hsd3);

  return x;
}

double FMT_Species_EOS::dfex(long pos, void* param)
{
  double x = effDensity(pos);
  double e = M_PI*x*hsd_*hsd_*hsd_/6;
  
  return eos_.fex(x) - ((FMT*) param)->get_fex(e)-0.5*avdw_*x*x;
}

void FMT_Species_EOS::calculateFundamentalMeasures(bool needsTensor)
{
  FMT_Species::calculateFundamentalMeasures(needsTensor);  
  const DFT_Vec_Complex &rho_k = density_->get_density_fourier();    
  eos_weighted_density_[0].convoluteWith(rho_k);          
}

void FMT_Species_EOS::set_fundamental_measure_derivatives(long pos, FundamentalMeasures &fm, void* param)
{
  FMT_Species::set_fundamental_measure_derivatives(pos,fm,param);
  eos_weighted_density_[0].Set_dPhi(pos,dfex(pos, param));
}

void FMT_Species_EOS::calculateForce(bool needsTensor, void *param)
{
  FMT_Species::calculateForce(needsTensor);

  double dV = getLattice().dV();
  int Nx    = getLattice().Nx();
  int Ny    = getLattice().Ny();
  int Nz    = getLattice().Nz();
  
  DFT_FFT dF_local(Nx,Ny,Nz);
  dF_local.Four().zeros();
  
  eos_weighted_density_[0].add_weight_schur_dPhi_to_arg(dF_local.Four());

  dF_local.do_fourier_2_real();
  dF_local.Real().MultBy(dV);
  addToForce(dF_local.cReal());   
}

