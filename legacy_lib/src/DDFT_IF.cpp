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


// Let D be the divergence operator (D = (d/dx, d/dy, d/dz)).

// The equation we are solving is
// (d/dt) rho_{r,t} = D^2 rho_{r,t} + D cdot rho_{r,t} D (delta F_ex/ delta rho_{rt}) 
// or, more succintly,
// (d/dt) rho_{r,t} = D^2 rho_{r,t} + R_{r}[rho_t]
// where we understand that R depends on rho at all spatial points.
// If we FFT we get
// (d/dt) rho_{k,t} = L_{k} rho_{k,t} + R_{k}[rho_t]
// and if R_=0 this can be integrated exactly. So we write
// rho_[k,t+dt] = exp(L_{k} dt) rho_{k,t} + exp(L_{k} (t+dt)) \int_{t}^{t+dt} exp(-L_{k} s) R_{k}[rho_s] ds
// and now we approximate R_{k}[rho_s] = R_{k}[rho_t] + ((s-t)/dt)*R_{k}[rho_{t+s}]
// giving
// rho_[k,t+dt] = exp(L_{k} dt) rho_{k,t} + exp(L_{k} dt) \int_{0}^{dt} exp(-L_{k} s) R_{k}[rho_s] ds
// rho_[k,t+dt] = exp(L_{k} dt) rho_{k,t} - ((exp(L_{k} dt)-1)/L_{k}) R_{k}[rho_t] +((exp(L_{k} dt) - 1 - L dt))/L_{k}^2) R_{k}[rho_{t+dt}]
// and these equations are solved iteratively. 

// This implementation presently only works for a single species!!!

void DDFT_IF::initialize()
{
  DDFT::initialize();
  
  F_ = dft_->calculateFreeEnergyAndDerivatives(false); //density_, 0.0, dF_,true);   

  cout << "Initial value of F = " << F_ << endl;
}

// Here we solve the nonlinear equation
// rho_[k,t+dt] = exp(L_{k} dt) rho_{k,t} - ((exp(L_{k} dt)-1)/L_{k}) R_{k}[rho_t] +((exp(L_{k} dt) - 1 - L dt))/L_{k}^2) R_{k}[rho_{t+dt}]
// by iteration.
// Define:
// d0_{k} = rho_{t,k}
// d1_{k} = rho_{t+dt,k}
// R0_{k} = R_{k}[d0]
// R1_{1} = R_{k}[d1]
// Then call to fftDiffusion(d1, R0,R1) to apply the diffusive prefactors and to update d1 according to 
//  d1 ->  exp(L_{k} dt) d0_{k} + R0 + R1
// Check for convergence and react accordingly

// This presently only works for a single species!!!

double DDFT_IF::step()
{
  bool bLessThanZero;
  cout << "===================================================================================================================" << endl;

  F_ = 0;
  try {
    F_ = dft_->calculateFreeEnergyAndDerivatives(false);   
  } catch( Eta_Too_Large_Exception &e) {
    throw e;
  }
  cout << "Initial F = " << F_ << endl;

  int Jspecies = 0;
  Species *species = dft_->getSpecies(Jspecies);
  
  const Lattice &lattice = species->getLattice();
  const Density &density = species->getDensity();
  
  DFT_Vec d0(density.getDensity());   
  DFT_FFT RHS0(lattice.Nx(), lattice.Ny(), lattice.Nz());
  calcNonlinearTerm(lattice, d0, species->getDF(), RHS0.Real());  
  RHS0.do_real_2_fourier();
    
  DFT_Vec d1(d0);
  DFT_FFT RHS1(lattice.Nx(), lattice.Ny(), lattice.Nz());

  double deviation = 1;

  double dt = dt_;

  bool reStart;
  bool decreased_time_step = false;
  
  do {
    reStart = false;
    double old_error = 0;
    
    for(int i=0;i<100 && deviation > tolerence_fixed_point_ && !reStart;i++)
      {
	species->set_density(d1);
	F_ = dft_->calculateFreeEnergyAndDerivatives(false); //density_,0.0, dF_,false);
	calcNonlinearTerm(lattice, d1, species->getDF(),RHS1.Real());
	RHS1.do_real_2_fourier();

	species->set_density(d0);       
	species->fft_density();
	
	deviation = fftDiffusion(density,d1,RHS0,RHS1);
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

  if(fixedBorder_)
    restore_values_on_border(lattice, d0, d1);
  
  species->set_density(d1);
  species->fft_density();
  calls_++;

  F_ = dft_->calculateFreeEnergyAndDerivatives(false); //density_,0.0, dF_,false);  

  cout << "F = " << F_ << endl;
    
  //  if(grace_)
  //    {
  //      double d00 = dF_.min()/density_.dV();
  //      double d11 = dF_.max()/density_.dV();   
  //      Display(F_,d00,d11,density_.getNumberAtoms());
  //    }
  
  return F_;
}

void DDFT_IF::restore_values_on_border(const Lattice &lattice, const DFT_Vec &d0, DFT_Vec &density)
{
  int Nx = lattice.Nx();
  int Ny = lattice.Ny();
  int Nz = lattice.Nz();
      
  for(int ix=0;ix<Nx;ix++)
    for(int iy=0;iy<Ny;iy++)
      {
	long i0 = lattice.pos(ix,iy,0);
	density.set(i0,d0.get(i0));
      }

  for(int ix=0;ix<Nx;ix++)
    for(int iz=0;iz<Nz;iz++)
      {
	long i0 = lattice.pos(ix,0,iz);
	density.set(i0,d0.get(i0));
      }
  for(int iy=0;iy<Ny;iy++)
    for(int iz=0;iz<Nz;iz++)
      {
	long i0 = lattice.pos(0,iy,iz);
	density.set(i0,d0.get(i0));
      }
}

/**
 This function takes the input density and calculates a new density, d1, by propagating the density a time step dt. The nonlinear terms, RHS0 and RHS1, are treated explicitly.
*/
double DDFT_IF::fftDiffusion(const Density& density, DFT_Vec &new_density, const DFT_FFT &RHS0, const DFT_FFT &RHS1)
{
  double dx = density.getDX();
  double dy = density.getDY();
  double dz = density.getDZ();

  double Dx = 1.0/(dx*dx);
  double Dy = 1.0/(dy*dy);
  double Dz = 1.0/(dz*dz);

  double dV = dx*dy*dz;

  int Nx = density.Nx();
  int Ny = density.Ny();
  int Nz = density.Nz();

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

	  complex<double> x = density.get_fourier_value(pos); 
	  
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

  long Ntot = density.Ntot();   
  for(unsigned i=0;i<new_density.size();i++)
    {
      deviation += (new_density.get(i)-work.Real().get(i)/Ntot)*(new_density.get(i)-work.Real().get(i)/Ntot);
      double u = fabs(new_density.get(i)-work.Real().get(i)/Ntot);
      if(u > maxdeviation) maxdeviation = u;
    }
  new_density.set(work.Real());
  new_density.MultBy(1.0/Ntot);
  return maxdeviation; //sqrt(deviation/Ntot);
}

// Let D be the divergence operator (D = (d/dx, d/dy, d/dz)). Then the equation we are solving is
// (d/dt) rho_{r,t} = D^2 rho_{r,t} + D cdot rho_{r,t} D (delta F_ex/ delta rho_{rt}) 
// Here, we calculate the last term using finite differences

void DDFT_IF::calcNonlinearTerm(const Lattice &lattice, const DFT_Vec &density, const DFT_Vec &dF, DFT_Vec &RHS1)
{
  double dx = lattice.getDX();
  double dy = lattice.getDY();
  double dz = lattice.getDZ();

  double Dx = 1.0/(dx*dx);
  double Dy = 1.0/(dy*dy);
  double Dz = 1.0/(dz*dz);

  double dV = dx*dy*dz;

  int Nx = lattice.Nx();
  int Ny = lattice.Ny();
  int Nz = lattice.Nz();
  
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
	  long i0 = lattice.get_PBC_Pos(ix,iy,iz);

	  long ipx = lattice.get_PBC_Pos(ix+1,iy,iz);
	  long imx = lattice.get_PBC_Pos(ix-1,iy,iz);

	  long ipy = lattice.get_PBC_Pos(ix,iy+1,iz);
	  long imy = lattice.get_PBC_Pos(ix,iy-1,iz);

	  long ipz = lattice.get_PBC_Pos(ix,iy,iz+1);
	  long imz = lattice.get_PBC_Pos(ix,iy,iz-1);
	    
	  double r0  = density.get(i0);

	  double rpx = density.get(ipx);
	  double rmx = density.get(imx);

	  double rpy = density.get(ipy);
	  double rmy = density.get(imy);
	    
	  double rpz = density.get(ipz);
	  double rmz = density.get(imz);

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

	  // eliminate the diffusive (e.g. ideal gas) contribution which is handled exactly in this algorithm
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








/*
double DDFT_IF::step_string(double &dt, Density &original_density, double self_consistency_threshold)
{
  int Nx = original_density.Nx();
  int Ny = original_density.Ny();
  int Nz = original_density.Nz();
  
  dt_ = dt;

  F_ = dft_.calculateFreeEnergyAndDerivatives(original_density,0.0, dF_,false);

  cout << "Initial F = " << F_ << endl;

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
	cout << "\tdeviation = " << deviation << " dt = " << dt_ << endl;

	// decrease time step and restart if density goes negative or if error is larger than previous step
	if(d1.min() < 0 || (i > 0 && old_error < deviation)) {reStart = true; dt_ /= 10; d1.set(d0); decreased_time_step = true;}

	old_error = deviation;	       	
      }
    if(!reStart && deviation > tolerence_fixed_point_)
      {reStart = true; dt_ /= 10; d1.set(d0); decreased_time_step = true;}
  } while(reStart);
  
  // Adaptive time-step: try to increase time step if the present one works 5 times 
  //if(decreased_time_step) successes_ = 0;
//  else successes_++;
//  if(successes_ >= 5 && dt_ < dtMax_) { dt_ = min(2*dt, dtMax_); successes_ = 0;}
  
  
  original_density.set(d1);

  if(fixedBorder_)
    restore_values_on_border(original_density, d0);

  original_density.doFFT();
  calls_++;

  //  F_ = dft_.calculateFreeEnergyAndDerivatives(original_density,0.0, dF_,false);  
    
  dt = dt_;

  return F_;
}
*/

