#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

#include <mgl2/mgl.h>
#include <mgl2/fltk.h>

#ifdef USE_OMP
#include <omp.h>
#endif

using namespace std;


#include "Minimizer.h"

void DDFT_IF_Open::initialize()
{
  DDFT::initialize();

  int Nx = density_.Nx();
  int Ny = density_.Ny();
  int Nz = density_.Nz();

  sin_Ntot_ = (Nx-1)*(Ny-1)*(Nz-1);
  sin_Norm_ = 8*Nx*Ny*Nz;

  sin_in_ = new double[sin_Ntot_];
  sin_out_ = new double[sin_Ntot_];

  sin_plan_ = fftw_plan_r2r_3d(Nx-1, Ny-1, Nz-1,
					sin_in_, sin_out_,
					 FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00,
					 FFTW_MEASURE);
  
  F_ = dft_.calculateFreeEnergyAndDerivatives(density_, 0.0, dF_,true);   
  DFT_Vec dummy;
  F_ += dft_.F_IdealGas(density_, dummy);
  F_ += dft_.F_External(density_,0.0,dummy);

  cout << "Initial value of F = " << F_ << endl;


  pack_for_sin_transform(density_.getData(),background_);
  fftw_execute(sin_plan_);
  memcpy(sin_in_, sin_out_, sin_Ntot_*sizeof(double));
  fftw_execute(sin_plan_);
  for(unsigned k=0;k<sin_Ntot_;k++) sin_out_[k] /= sin_Norm_;
  double *d = new double[density_.Ntot()];
  unpack_after_transform(d,background_);

}

double DDFT_IF_Open::step_string(double &dt, Density &original_density, bool verbose)
{
  int Nx = original_density.Nx();
  int Ny = original_density.Ny();
  int Nz = original_density.Nz();

  long Ntot = original_density.Ntot();
  
  dt_ = dt;

  F_ = dft_.calculateFreeEnergyAndDerivatives(original_density,0.0, dF_,false);

  if(verbose) cout << "Initial F = " << F_ << endl;

  const DFT_Vec &d0 = original_density.getDensity();   
  DFT_Vec RHS0(Ntot);
  calcNonlinearTerm(d0, dF_,RHS0);
  double *RHS0_sin_transform = new double[sin_Ntot_];

  pack_for_sin_transform(RHS0.memptr(),0.0);
  fftw_execute(sin_plan_);
  memcpy(RHS0_sin_transform,sin_out_,sin_Ntot_*sizeof(double));
    
  DFT_Vec d1(d0);
  DFT_Vec RHS1; RHS1.zeros(Ntot);  
  double *RHS1_sin_transform = new double[sin_Ntot_];

  double deviation = 1;

  bool reStart;
  bool decreased_time_step = false;
  
  do {
    reStart = false;
    double old_error = 0;
    
    for(int i=0;i<20 && deviation > tolerence_fixed_point_ && !reStart;i++)
      {
	density_.set(d1);
	F_ = dft_.calculateFreeEnergyAndDerivatives(density_,0.0, dF_,false);
	calcNonlinearTerm(d1, dF_,RHS1);

	pack_for_sin_transform(RHS1.memptr(),0.0);
	fftw_execute(sin_plan_);
	memcpy(RHS1_sin_transform,sin_out_,sin_Ntot_*sizeof(double));
	
	density_.set(d0);       
	density_.doFFT();
	
	deviation = fftDiffusion(d1,RHS0_sin_transform,RHS1_sin_transform);
	if(verbose) cout << "\tdeviation = " << deviation << " dt = " << dt_ << endl;

	// decrease time step and restart if density goes negative or if error is larger than previous step
	if(d1.min() < 0 || (i > 0 && old_error < deviation))
	  {
	    reStart = true; dt_ /= (fabs(d1.min()) > 0.1 ? 10 : 2) ; d1.set(d0); decreased_time_step = true;
	  }
	old_error = deviation;	       	
      }
    if(!reStart && deviation > tolerence_fixed_point_)
      {reStart = true; dt_ /= 10; d1.set(d0); decreased_time_step = true;}
  } while(reStart);
  
  
  original_density.set(d1);
  original_density.doFFT();
  calls_++;

  dt = dt_;

  delete RHS0_sin_transform;
  delete RHS1_sin_transform;

  return F_;
}



double DDFT_IF_Open::step()
{
  double olddt = dt_;
  double dt = dt_;
  F_ = step_string(dt, density_, true);
  dt_ = olddt;
  
    // Adaptive time-step: try to increase time step if the present one works 5 times 
  if(dt < dt_*1.01) successes_ = 0;
  else successes_++;
  if(successes_ >= 5 && dt_ < dtMax_) { dt_ = min(2*dt, dtMax_); successes_ = 0;}

  
  /*
  bool bLessThanZero;
  cout << "===================================================================================================================" << endl;

  F_ = dft_.calculateFreeEnergyAndDerivatives(density_,0.0, dF_,false);

  cout << "Initial F = " << F_ << endl;

  DFT_Vec d0(density_.getDensity());   
  DFT_FFT RHS0(density_.Nx(), density_.Ny(), density_.Nz());
  calcNonlinearTerm(d0, dF_,RHS0.Real());
  RHS0.do_real_2_fourier();
    
  DFT_Vec d1(d0);
  DFT_FFT RHS1(density_.Nx(), density_.Ny(), density_.Nz());

  double deviation = 1;

  double dt = dt_;

  bool reStart;
  bool decreased_time_step = false;
  
  do {
    reStart = false;
    double old_error = 0;
    
    for(int i=0;i<100 && deviation > tolerence_fixed_point_ && !reStart;i++)
      {
	density_.set(d1);
	F_ = dft_.calculateFreeEnergyAndDerivatives(density_,0.0, dF_,false);
	calcNonlinearTerm(d1, dF_,RHS1.Real());
	RHS1.do_real_2_fourier();

	density_.set(d0);       
	density_.doFFT();
	
	deviation = fftDiffusion(density_,d1,RHS0,RHS1);
	cout << "\tdeviation = " << deviation << " dt = " << dt_ << endl;

	// decrease time step and restart if density goes negative or if error is larger than previous step
	if(d1.min() < 0 || (i > 0 && old_error < deviation)) {reStart = true; dt_ /= 10; d1.set(d0); decreased_time_step = true;}

	old_error = deviation;	       	
      }
    if(!reStart && deviation > tolerence_fixed_point_)
      {reStart = true; dt_ /= 10; d1.set(d0); decreased_time_step = true;}
  } while(reStart);

  // Adaptive time-step: try to increase time step if the present one works 5 times 
  if(decreased_time_step) successes_ = 0;
  else successes_++;
  if(successes_ >= 5 && dt_ < dtMax_) { dt_ = min(2*dt, dtMax_); successes_ = 0;}

  density_.set(d1);

  density_.doFFT();
  calls_++;

  F_ = dft_.calculateFreeEnergyAndDerivatives(density_,0.0, dF_,false);  

  cout << "F = " << F_ << endl;
    
  double d00 = dF_.min()/density_.dV();
  double d11 = dF_.max()/density_.dV();
    
  if(grace_)
    Display(F_,d00,d11,density_.getNumberAtoms());
  */
  return F_;
}


/**
 This function takes the input density and calculates a new density, d1, by propagating the density a time step dt. The nonlinear terms, RHS0 and RHS1, are treated explicitly.
*/
double DDFT_IF_Open::fftDiffusion(DFT_Vec &d1, const double *RHS0_sin_transform, const double *RHS1_sin_transform)
{
  double dx = density_.getDX();
  double dy = density_.getDY();
  double dz = density_.getDZ();

  double Dx = 1.0/(dx*dx);
  double Dy = 1.0/(dy*dy);
  double Dz = 1.0/(dz*dz);

  double dV = dx*dy*dz;

  int Nx = density_.Nx();
  int Ny = density_.Ny();
  int Nz = density_.Nz();

  double deviation = 0;
  double maxdeviation = 0;
  
  pack_for_sin_transform(density_.getData(),background_);
  fftw_execute(sin_plan_);
  
  // fft of density is now in sin_out_;
  // we construct the rhs in sin_in_ since we will fft it to get the density in real space
  
  unsigned pos = 0;
  for(int ix=0;ix<Nx-1;ix++)
    for(int iy=0;iy<Ny-1;iy++)
      for(int iz=0;iz<Nz-1;iz++)
	{
	  double kx = (M_PI*(ix+1))/Nx;
	  double ky = (M_PI*(iy+1))/Ny;
	  double kz = (M_PI*(iz+1))/Nz;

	  double facx = 2*Dx*(cos(kx)-1);
	  double facy = 2*Dy*(cos(ky)-1);
	  double facz = 2*Dz*(cos(kz)-1);

	  double Lambda = facx+facy+facz;

	  double x = sin_out_[pos];
	  
	  double fac = exp(Lambda*dt_);
	  x *= fac;
	  x += ((fac-1)/Lambda)*RHS0_sin_transform[pos];
	  x += ((fac-1-dt_*Lambda)/(Lambda*Lambda*dt_))*(RHS1_sin_transform[pos]-RHS0_sin_transform[pos]);

	  sin_in_[pos] = x;
	  pos++;
	}

  fftw_execute(sin_plan_);

  pos = 0;
  for(int ix=1;ix<Nx;ix++)
    for(int iy=1;iy<Ny;iy++)
      for(int iz=1;iz<Nz;iz++)
	{
	  
	  double d = (sin_out_[pos++] /= sin_Norm_)+background_;

	  double d1val = d1.get(density_.pos(ix,iy,iz));
	  
	  deviation += (d1val-d)*(d1val-d);
	  double u = fabs(d1val-d);
	  if(u > maxdeviation) maxdeviation = u;	  
	}
  unpack_after_transform(d1.memptr(), background_);
  return maxdeviation; 
}

void DDFT_IF_Open::pack_for_sin_transform(const double *x, double val)
{
  int Nx = density_.Nx();
  int Ny = density_.Ny();
  int Nz = density_.Nz();

  unsigned pos = 0;
  for(int ix=1;ix<Nx;ix++)
    for(int iy=1;iy<Ny;iy++)
      for(int iz=1;iz<Nz;iz++)
	  sin_in_[pos++] = x[density_.pos(ix,iy,iz)] - val;
}

void DDFT_IF_Open::unpack_after_transform(double *x, double val)
{
  int Nx = density_.Nx();
  int Ny = density_.Ny();
  int Nz = density_.Nz();

  unsigned pos = 0;
  for(int ix=0;ix<Nx;ix++)
    for(int iy=0;iy<Ny;iy++)
      for(int iz=0;iz<Nz;iz++)
	{
	  if(ix == 0 || iy == 0 || iz == 0)
	    x[density_.pos(ix,iy,iz)] = val;
	  else
	    {
	      double v = sin_out_[pos++] + val;
	      x[density_.pos(ix,iy,iz)] = v;
	    }
	}
}

void DDFT_IF_Open::calcNonlinearTerm(const DFT_Vec &d2, const DFT_Vec &dF, DFT_Vec &RHS1)
{
  double dx = density_.getDX();
  double dy = density_.getDY();
  double dz = density_.getDZ();

  double Dx = 1.0/(dx*dx);
  double Dy = 1.0/(dy*dy);
  double Dz = 1.0/(dz*dz);

  double dV = dx*dy*dz;

  int Nx = density_.Nx();
  int Ny = density_.Ny();
  int Nz = density_.Nz();
  
  int chunk = Ny/10;
  int iy;

#pragma omp parallel for			\
  shared(chunk,dx,dy,dz,Dx,Dy,Dz,dV)		\
  private(iy)					\
  schedule(static,chunk)			  
  for(iy = 0;iy<Ny;iy++)
    for(int iz = 0;iz<Nz;iz++)
      for(int ix=0;ix<Nx;ix++)
	{
	  long i0 = density_.get_PBC_Pos(ix,iy,iz);

	  long ipx = density_.get_PBC_Pos(ix+1,iy,iz);
	  long imx = density_.get_PBC_Pos(ix-1,iy,iz);

	  long ipy = density_.get_PBC_Pos(ix,iy+1,iz);
	  long imy = density_.get_PBC_Pos(ix,iy-1,iz);

	  long ipz = density_.get_PBC_Pos(ix,iy,iz+1);
	  long imz = density_.get_PBC_Pos(ix,iy,iz-1);
	    
	  double r0 = d2.get(i0);

	  double rpx = d2.get(ipx);
	  double rmx = d2.get(imx);

	  double rpy = d2.get(ipy);
	  double rmy = d2.get(imy);
	    
	  double rpz = d2.get(ipz);
	  double rmz = d2.get(imz);

	  double f0 = dF.get(i0);

	  double fpx = dF.get(ipx);
	  double fmx = dF.get(imx);

	  double fpy = dF.get(ipy);
	  double fmy = dF.get(imy);

	  double fpz = dF.get(ipz);
	  double fmz = dF.get(imz);

	  double RHS_F = Dx*((rpx+r0)*(fpx-f0)-(r0+rmx)*(f0-fmx))
	    +Dy*((rpy+r0)*(fpy-f0)-(r0+rmy)*(f0-fmy))
	    +Dz*((rpz+r0)*(fpz-f0)-(r0+rmz)*(f0-fmz));

	  // factor 1/2 because of density average used in prev line
	  // factor dV because f0, etc carries a dV
	  RHS_F *= 1.0/(2*dV);

	  RHS_F -= Dx*(rpx+rmx-2*r0)+Dy*(rpy+rmy-2*r0)+Dz*(rpz+rmz-2*r0);

	  // N.B. getFieldDeriv returns the derivative at iz-0.5.

	  // Attention: make sure the right density is used in this evaluation ....
	  //	  int Nz = density_.Nz();
	  //	  double dvpz = density_.getFieldDeriv(0,0,(iz+1 > Nz ? 0 : iz+1),2);
	  //	  double dvmz = density_.getFieldDeriv(0,0,iz,2);
	  //	  RHS_F += dz*Dz*((rpz+r0)*dvpz-(r0+rmz)*dvmz)/2;

	  RHS1.set(i0,RHS_F);
	}
}








