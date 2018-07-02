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

void DDFT_IF::initialize()
{
  DDFT::initialize();
  
  F_ = dft_.calculateFreeEnergyAndDerivatives(density_, 0.0, dF_,true);   
  DFT_Vec dummy;
  F_ += dft_.F_IdealGas(density_, dummy);
  F_ += dft_.F_External(density_,0.0,dummy);

  //  cout << "Initial value of F = " << F_ << endl;
}

double DDFT_IF::step_string(double &dt, Density &original_density, double self_consistency_threshold, bool verbose)
{
  int Nx = original_density.Nx();
  int Ny = original_density.Ny();
  int Nz = original_density.Nz();
  
  dt_ = dt;

  F_ = dft_.calculateFreeEnergyAndDerivatives(original_density,0.0, dF_,false);

  if(verbose) cout << "Initial F = " << F_ << endl;

  const DFT_Vec &d0 = original_density.getDensity();   
  DFT_FFT RHS0(Nx,Ny,Nz);
  calcNonlinearTerm(d0, dF_,RHS0.Real());
  RHS0.do_real_2_fourier();
    
  DFT_Vec d1(d0);
  DFT_FFT RHS1(Nx,Ny,Nz);

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
	calcNonlinearTerm(d1, dF_,RHS1.Real());
	RHS1.do_real_2_fourier();

	density_.set(d0);       
	density_.doFFT();
	
	deviation = fftDiffusion(d1,RHS0,RHS1);
	if(verbose) cout << "\tdeviation = " << deviation << " dt = " << dt_ << endl;

	// decrease time step and restart if density goes negative or if error is larger than previous step
	if(d1.min() < 0 || (i > 0 && old_error < deviation)) {
	  cout << "d0.min = " << d0.min() << " d1.min = " << d1.min() << endl;
	  reStart = true; dt_ /= 10; d1.set(d0); decreased_time_step = true;}

	old_error = deviation;	       	
      }
    if(!reStart && deviation > tolerence_fixed_point_)
      {reStart = true; dt_ /= 10; d1.set(d0); decreased_time_step = true;}
  } while(reStart);
  /*
  // Adaptive time-step: try to increase time step if the present one works 5 times 
  if(decreased_time_step) successes_ = 0;
  else successes_++;
  if(successes_ >= 5 && dt_ < dtMax_) { dt_ = min(2*dt, dtMax_); successes_ = 0;}
  */
  
  original_density.set(d1);

  if(fixedBorder_)
    restore_values_on_border(original_density, d0);
  
  original_density.doFFT();
  calls_++;

  //  F_ = dft_.calculateFreeEnergyAndDerivatives(original_density,0.0, dF_,false);  
    
  dt = dt_;

  return F_;
}



double DDFT_IF::step()
{
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
	
	deviation = fftDiffusion(d1,RHS0,RHS1);
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

  if(fixedBorder_)
    restore_values_on_border(density_,d0);


  density_.doFFT();
  calls_++;

  F_ = dft_.calculateFreeEnergyAndDerivatives(density_,0.0, dF_,false);  

  cout << "F = " << F_ << endl;
    
  double d00 = dF_.min()/density_.dV();
  double d11 = dF_.max()/density_.dV();
    
  if(grace_)
    Display(F_,d00,d11,density_.getNumberAtoms());

  return F_;
}

void DDFT_IF::restore_values_on_border(Density &density, const DFT_Vec &d0)
{
  
  int Nx = density.Nx();
  int Ny = density.Ny();
  int Nz = density.Nz();


  int jx=0,jy=0,jz=0;
  double max_change = 0;
  double d_before = 0;
  //  cout << "Before restoring values on border: " << density.getDensity(Nx-1,(Ny-1)/2,(Nz-1)/2);    





  
  for(int ix=0;ix<Nx;ix++)
    for(int iy=0;iy<Ny;iy++)
      {
	long i0 = density.pos(ix,iy,0);

	double dd = fabs(density.getDensity(i0) - d0.get(i0));
	if(dd > max_change){ max_change = dd; jx = ix; jy = iy; jz = 0; d_before = density.getDensity(i0);}
	
	density.set_Density_Elem(i0,d0.get(i0));
	
      }

  for(int ix=0;ix<Nx;ix++)
    for(int iz=0;iz<Nz;iz++)
      {
	long i0 = density.pos(ix,0,iz);
	
	double dd = fabs(density.getDensity(i0) - d0.get(i0));
	if(dd > max_change){ max_change = dd; jx = ix; jy = 0; jz = iz; d_before = density.getDensity(i0);}
	
	density.set_Density_Elem(i0,d0.get(i0));
      }
  for(int iy=0;iy<Ny;iy++)
    for(int iz=0;iz<Nz;iz++)
      {
	long i0 = density.pos(0,iy,iz);
	
	double dd = fabs(density.getDensity(i0) - d0.get(i0));
	if(dd > max_change){ max_change = dd; jx = 0; jy = iy; jz = iz; d_before = density.getDensity(i0);}
	
	density.set_Density_Elem(i0,d0.get(i0));
      }

  cout << "With modified_ = " << modified_ << " max change on boundary was " << max_change << " at position " << jx << "," << jy << "," << jz << " where density was " << d0.get(density.pos(jx,jy,jz)) << " and became " <<  d_before << endl;
}

/**
 This function takes the input density and calculates a new density, d1, by propagating the density a time step dt. The nonlinear terms, RHS0 and RHS1, are treated explicitly.
*/
double DDFT_IF::fftDiffusion(DFT_Vec &d1, const DFT_FFT &RHS0, const DFT_FFT &RHS1)
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
  
  DFT_FFT work(Nx,Ny,Nz);

  DFT_Vec_Complex &cwork = work.Four();

  for(int ix = 0;ix<Nx;ix++)
    for(int iy=0;iy<Ny;iy++)
      for(int iz=0;iz<Nz/2+1;iz++)
	{
	  double kx = (2*M_PI*ix)/Nx;
	  double ky = (2*M_PI*iy)/Ny;
	  double kz = (2*M_PI*iz)/Nz;

	  double facx = 2*Dx*(cos(kx)-1);
	  double facy = 2*Dy*(cos(ky)-1);
	  double facz = 2*Dz*(cos(kz)-1);

	  double Lambda = facx+facy+facz;

	  unsigned pos = iz + (1+Nz/2)*(iy +  Ny*ix);

	  complex<double> x = density_.DensityK(pos); 
	  
	  if(pos > 0)
	    {
	      double fac = exp(Lambda*dt_);
	      x *= fac;
	      x += ((fac-1)/Lambda)*RHS0.cFour().get(pos);
	      x += ((fac-1-dt_*Lambda)/(Lambda*Lambda*dt_))*(RHS1.cFour().get(pos)-RHS0.cFour().get(pos));
	    }
	  cwork.set(pos,x);
	}
  
  work.do_fourier_2_real();

  long Ntot = density_.Ntot();   
  for(unsigned i=0;i<d1.size();i++)
    {
      deviation += (d1.get(i)-work.Real().get(i)/Ntot)*(d1.get(i)-work.Real().get(i)/Ntot);
      double u = fabs(d1.get(i)-work.Real().get(i)/Ntot);
      if(u > maxdeviation) maxdeviation = u;
    }
  d1.set(work.Real());
  d1.multBy(1.0/Ntot);
  return maxdeviation; //sqrt(deviation/Ntot);
}


void DDFT_IF::calcNonlinearTerm(const DFT_Vec &d2, const DFT_Vec &dF, DFT_Vec &RHS1)
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

	  // In the modified version, on the boundaries the RHS is thrown away and only the (negative) diffusive part kept
	  if(modified_)
	    if(ix == 0 || iy == 0 || iz == 0)
	      RHS_F = 0.0;

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








