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



// This class adds a new weighted density which is the packing fraction over an extended volume defined by the effective hsd = D_EOS*hsd.
// Typically, one uses D_EOS = 2.
// It is intended to shift the free energy of the uniform system from F_dft(rho) to F_EOS(rho).
// This is done by defining dF_EOS(x) = F_dft(x)-F_EOS(x) and then calculating
// the contribution to the free energy as sum_I dF_eos(rho_eos(I))dV
// Here, rho_eos(I) is determined from the new packing fraction as eta_eos(I)*(6/M_PI)/(D_EOS*hsd)^3.
//
// Note that F_dft(x) = F_HS(x) + 0.5*avdW*x*x.
//
//
// For the forces, we need
//        d/drho_I sum_J dfex(rho_eos_J) = sum_I dfex'(rho_eos_J) w^{EOS}_JI
//

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

//  rho_eos(I) = eta_eos(I)*(6/M_PI)/(D_EOS*hsd)^3.
double FMT_Species_EOS::effDensity(long I)
{
  double eta = eos_weighted_density_[0].real(I);
  double hsd3 = hsd_*hsd_*hsd_*D_EOS_*D_EOS_*D_EOS_;
  double x   = 6*eta/(M_PI*hsd3);

  return x;
}

// dF_EOS(x) = F_EOS(x)-F_dft(x) = F_EOS(x)_excess - F_HS_excess(x) - 0.5*a*x*x 
double FMT_Species_EOS::dfex(long pos, void* param)
{
  double x = effDensity(pos);
  double eta = M_PI*x*hsd_*hsd_*hsd_/6;  
  double fdft = x*((FMT*) param)->get_fex(eta) + avdw_*x*x;

  return eos_.fex(x) - fdft;
}

void FMT_Species_EOS::calculateFundamentalMeasures(bool needsTensor)
{
  FMT_Species::calculateFundamentalMeasures(needsTensor);  
  const DFT_Vec_Complex &rho_k = density_->get_density_fourier();    
  eos_weighted_density_[0].convoluteWith(rho_k);          
}



// Here, we need d
void FMT_Species_EOS::set_fundamental_measure_derivatives(long pos, FundamentalMeasures &fm, void* param)
{
  if(pos == 10) cout << "Set EOS phi" << endl;
  FMT_Species::set_fundamental_measure_derivatives(pos,fm,param);

  double x = effDensity(pos);
  double eta = M_PI*x*hsd_*hsd_*hsd_/6;  
  double dfdft = ((FMT*) param)->get_fex(eta) + eta*((FMT*) param)->get_dfex_deta(eta) + 2*avdw_*x;
  // HERE
  // convert to df/deta_EOS
  eos_weighted_density_[0].Set_dPhi(pos,(eos_.f1ex(x) - dfdft)*6/(M_PI*hsd_*hsd_*hsd_*D_EOS_*D_EOS_*D_EOS_));
}

void FMT_Species_EOS::calculateForce(bool needsTensor, void *param)
{
  cout << "EOS force" << endl;
  
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

