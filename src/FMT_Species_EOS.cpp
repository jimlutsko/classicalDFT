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
//        d/drho_I sum_J dfex(rho_eos_J) = sum_I (d/d_eta_EOS_J) dfex(rho_eos_J) w^{EOS}_JI
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

double FMT_Species_EOS::get_bulk_dfex(double x, const void *param) const
{
  double eta = M_PI*x*hsd_*hsd_*hsd_/6;
  FMT &fmt   = *((FMT*) param);  
  double fdft = x*fmt.get_fex(eta) + avdw_*x*x;
  return eos_.fex(x) - fdft;
}

double FMT_Species_EOS::get_bulk_ddfex_deta(double x, const void *param) const
{
  double HSeta = M_PI*x*hsd_*hsd_*hsd_/6;
  FMT &fmt     = *((FMT*) param);

  double dfdft_dx = fmt.get_fex(HSeta) + HSeta*fmt.get_dfex_deta(HSeta) + 2*avdw_*x;

  double dx_deta  = 6/(M_PI*hsd_*hsd_*hsd_*D_EOS_*D_EOS_*D_EOS_);

  return dx_deta*(eos_.f1ex(x) - dfdft_dx);
}

double FMT_Species_EOS::get_bulk_d2dfex_deta2(double x, const void *param) const
{
  double HSeta = M_PI*x*hsd_*hsd_*hsd_/6;
  FMT &fmt     = *((FMT*) param);

  double dHSeta_dx = M_PI*hsd_*hsd_*hsd_/6;

  double d2fdft_dx2 = dHSeta_dx*(2*fmt.get_dfex_deta(HSeta) + HSeta*fmt.get_d2fex_deta2(HSeta)) + 2*avdw_;

  double dx_deta = 6/(M_PI*hsd_*hsd_*hsd_*D_EOS_*D_EOS_*D_EOS_);
  
  return dx_deta*dx_deta*(eos_.f2ex(x) - d2fdft_dx2);
}

// dF_EOS(x) = F_EOS(x)-F_dft(x) = F_EOS(x)_excess - F_HS_excess(x) - 0.5*a*x*x
double FMT_Species_EOS::dfex(long pos, void* param)
{
  double x   = effDensity(pos);
  return get_bulk_dfex(x,param);
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
  FMT_Species::set_fundamental_measure_derivatives(pos,fm,param);

  double x   = effDensity(pos);
  //  double eta = M_PI*x*hsd_*hsd_*hsd_/6;
  //  FMT &fmt   = *((FMT*) param);

  //  double dfdft = (fmt.get_fex(eta) + eta*fmt.get_dfex_deta(eta)) + 2*avdw_*x;

  // convert to df/deta_EOS
  //  double fac = 6/(M_PI*hsd_*hsd_*hsd_*D_EOS_*D_EOS_*D_EOS_);
  //  eos_weighted_density_[0].Set_dPhi(pos,(eos_.f1ex(x) - dfdft)*fac);
  eos_weighted_density_[0].Set_dPhi(pos,get_bulk_ddfex_deta(x, param));
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

// For this we need sum_J sum_K dV (d2dfex(K)/deta(K) deta(K)) w_eta(K-I) w_eta(K-J)v(J)
// Define psi(K) = dV conv(w_eta,v;K) = sum_J w_eta(K-J)V(J) dV
// Then we need sum_K (d2dfex(K)/deta(K) deta(K)) w_eta(K-I) psi(K)
// or conv((d2dfex(K)/deta(K) deta(K))psi(K), w; I)
// or we define Lambda(K) = (d2PHI(K)/deta(K) deta(K))psi(K)
// and then calculate conv(Lam, w; I)
//
// d2F is the force FOR THIS SPECIES
void FMT_Species_EOS::add_second_derivative(const DFT_FFT &v, DFT_Vec &d2F, const void *param)
{
  const Density &density1 = getDensity();

  long Ntot    = density1.Ntot();
  int Nx       = density1.Nx();
  int Ny       = density1.Ny();
  int Nz       = density1.Nz();
  double dV    = density1.dV();
  double hsd   = getHSD();
  
  // Some working space
  DFT_FFT result(Nx,Ny,Nz);

  // Get Psi: psi(K) = conv(w_eta,v;K) = sum_J w_eta(K-J)V(J)
  DFT_FFT Psi(Nx,Ny,Nz);
  
  // N.B. the fmt weights have all necessary normalization factors built in, including dV
  bool bConjugate = false;
  convolute_eta_weight_with(v, result, bConjugate);
  Psi.Real().IncrementBy(result.Real());

  // Get Lambda: Lambda(K) = (d2dfex(K)/deta(K) deta(K))psi(K)  
  DFT_FFT Lambda(Nx,Ny,Nz);
  
  long pos;
#ifdef USE_OMP
#pragma omp parallel for  private(pos) schedule(static)
#endif
  for(pos = 0; pos < Ntot; pos++)
    Lambda.Real().set(pos, get_bulk_d2dfex_deta2(effDensity(pos), param)*Psi.Real().get(pos));

  // Do the last convolution to get d2F
  Lambda.do_real_2_fourier();

  //the fmt weights have all necessary normalization factors built in, including dV
  bConjugate = true;
  convolute_eta_weight_with(Lambda, result, bConjugate);
  d2F.IncrementBy_Scaled_Vector(result.Real(),dV);
}
