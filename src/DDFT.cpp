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

#define SMALL_DENSITY 1e-8

#include "Minimizer.h"

void DDFT::initialize()
{
  Minimizer::initialize();

  successes_ = 0;

  
  F_ = dft_.calculateFreeEnergyAndDerivatives(density_, 0.0, dF_,true);   
  DFT_Vec dummy;
  F_ += dft_.F_IdealGas(density_, dummy);
  F_ += dft_.F_External(density_,0.0,dummy);

  cout << "Initial value of F = " << F_ << endl;
}


void DDFT::sub_step_x(DFT_Vec &ynew, const Density& original_density, bool bFixedBoundaries)
{
  double dx = original_density.getDX();
  double dy = original_density.getDY();
  double dz = original_density.getDZ();

  double Dx = dt_/(dx*dx);
  double Dy = dt_/(dy*dy);
  double Dz = dt_/(dz*dz);

  double dV = dx*dy*dz;

  int Nx = original_density.Nx();
  
  int chunk = original_density.Ny()/10;
  int iy;

#pragma omp parallel for			\
  shared(chunk,dx,dy,dz,Dx,Dy,Dz,dV)		\
  private(iy)					\
  schedule(static,chunk)			  
  for(iy = 0;iy<original_density.Ny();iy++)
    for(int iz = 0;iz<original_density.Nz();iz++)
      {
	// set up source term for tridiagonal equations
	DFT_Vec RHS(Nx);

	for(int ix=0;ix<Nx;ix++)
	  {
	    long i0 = original_density.get_PBC_Pos(ix,iy,iz);

	    long ipx = original_density.get_PBC_Pos(ix+1,iy,iz);
	    long imx = original_density.get_PBC_Pos(ix-1,iy,iz);

	    long ipy = original_density.get_PBC_Pos(ix,iy+1,iz);
	    long imy = original_density.get_PBC_Pos(ix,iy-1,iz);

	    long ipz = original_density.get_PBC_Pos(ix,iy,iz+1);
	    long imz = original_density.get_PBC_Pos(ix,iy,iz-1);
	    
	    double r0 = original_density.getDensity(i0);

	    double rpx = original_density.getDensity(ipx);
	    double rmx = original_density.getDensity(imx);

	    double rpy = original_density.getDensity(ipy);
	    double rmy = original_density.getDensity(imy);
	    
	    double rpz = original_density.getDensity(ipz);
	    double rmz = original_density.getDensity(imz);

	    double f0 = dF_.get(i0);

	    double fpx = dF_.get(ipx);
	    double fmx = dF_.get(imx);

	    double fpy = dF_.get(ipy);
	    double fmy = dF_.get(imy);

	    double fpz = dF_.get(ipz);
	    double fmz = dF_.get(imz);

	    double RHS_F = Dx*((rpx+r0)*(fpx-f0)-(r0+rmx)*(f0-fmx))
	      +Dy*((rpy+r0)*(fpy-f0)-(r0+rmy)*(f0-fmy))
	      +Dz*((rpz+r0)*(fpz-f0)-(r0+rmz)*(f0-fmz));

	    // factor 1/2 because of density average used in prev line
	    // factor dV because f0, etc carries a dV
	    RHS_F *= 1.0/(2*dV);

	    // N.B. getFieldDeriv returns the derivative at iz-0.5.
	    int Nz = original_density.Nz();
	    double dvpz = original_density.getFieldDeriv(0,0,(iz+1 > Nz ? 0 : iz+1),2);
	    double dvmz = original_density.getFieldDeriv(0,0,iz,2);
	    RHS_F += dz*Dz*((rpz+r0)*dvpz-(r0+rmz)*dvmz)/2;

	    RHS_F = 0.0;

	    // This is Douglas-Rachford: a fully implicit scheme
	    //   : the LHS is 1-Dx d2/dx2
	    //RHS.set(ix, r0+ Dy*(rpy+rmy-2*r0) + Dz*(rpz+rmz-2*r0) + RHS_F);

	    // This is Douglas-Gunn: the ADI version of Crank-Nicholson
	    // : the LHS is 1 - (1/2) DX d2/dx2
	    RHS.set(ix, r0 + 0.5*Dx*(rpx+rmx-2*r0) + Dy*(rpy+rmy-2*r0) + Dz*(rpz+rmz-2*r0) + RHS_F);

	    if(bFixedBoundaries)
	      if(ix == 0 || ix == Nx-1) RHS.set(ix,r0);
		
	  }
	if(!bFixedBoundaries)
	  solv_periodic_tridiag(RHS,0.5*Dx);
	else {
	  DFT_Vec b(Nx);

	  b.set(0,1.0); 				
	  b.set(Nx-1,1.0); 
	  
	  for(int kk=1;kk<Nx-1;kk++)
	    b.set(kk,1+2*(0.5*Dx));

	  solv_tridiag(b,RHS,0.5*Dx,bFixedBoundaries);
	}
	for(int ix=0;ix<Nx;ix++)
	  {
	    long i0 = original_density.get_PBC_Pos(ix,iy,iz);		    
	    if(bFixedBoundaries && (iy == 0 || iy == original_density.Ny()-1 || iz == 0 || iz == original_density.Nz()-1))
	      ynew.set(i0,original_density.getDensity(i0));
	    else ynew.set(i0,RHS.get(ix));
	  }	
      }
}

// ynew comes in with the result of the first step
void DDFT::sub_step_y(DFT_Vec &ynew, const Density &original_density, bool bFixedBoundaries)
{
  double dy = original_density.getDY();
  double Dy = dt_/(dy*dy);

  int Ny = original_density.Ny();
  
  int chunk = Ny/10;
  int ix;

#pragma omp parallel for			\
  shared(chunk,Dy)		\
  private(ix)					\
  schedule(static,chunk)			  
  for(ix = 0;ix<original_density.Nx();ix++)
    for(int iz = 0;iz<original_density.Nz();iz++)
      {
	// set up source term for tridiagonal equations
	DFT_Vec RHS(Ny);

	for(int iy=0;iy<Ny;iy++)
	  {
	    long i0 = original_density.get_PBC_Pos(ix,iy,iz);

	    long ipy = original_density.get_PBC_Pos(ix,iy+1,iz);
	    long imy = original_density.get_PBC_Pos(ix,iy-1,iz);

	    double r0 = original_density.getDensity(i0);

	    double rpy = original_density.getDensity(ipy);
	    double rmy = original_density.getDensity(imy);

	    //	    RHS.set(iy, ynew.get(i0) - Dy*(rpy+rmy-2*r0));
	    RHS.set(iy, ynew.get(i0) - 0.5*Dy*(rpy+rmy-2*r0));
	    if(bFixedBoundaries)
	      if(iy == 0 || iy == Ny-1) RHS.set(iy,r0);	    
	  }
	if(!bFixedBoundaries)
	  solv_periodic_tridiag(RHS,0.5*Dy);
	else {
	  DFT_Vec b(original_density.Ny());
	  b.set(0,1.0);
	  b.set(original_density.Ny()-1,1.0);
	  for(int kk=1;kk<original_density.Ny()-1;kk++)
	    b.set(kk,1+2*(0.5*Dy));
	  solv_tridiag(b,RHS,0.5*Dy, bFixedBoundaries);
	}
	
	for(int iy=0;iy<original_density.Ny();iy++)
	  {
	    long i0 = original_density.get_PBC_Pos(ix,iy,iz);		    
	    //	    ynew.set(i0,RHS.get(iy));
	    if(bFixedBoundaries && (ix == 0 || ix == original_density.Nx()-1 || iz == 0 || iz == original_density.Nz()-1))
	      ynew.set(i0,original_density.getDensity(i0));
	    else ynew.set(i0,RHS.get(iy));	    
	  }
      }
  return;
}

static double deviation = 0;


bool DDFT::sub_step_z(DFT_Vec &ynew, const Density &original_density, bool bUseForces, bool bFixedBoundaries) //const DFT_Vec &oldF)
{
  double dx = original_density.getDX();
  double dy = original_density.getDY();
  double dz = original_density.getDZ();

  double Dx = dt_/(dx*dx);
  double Dy = dt_/(dy*dy);
  double Dz = dt_/(dz*dz);

  double dV = dx*dy*dz;

  int Nz = original_density.Nz();
  
  bool bLessThanZero = false;
    
  int chunk = original_density.Nx()/10;
  int ix;

#pragma omp parallel for			\
  shared(chunk,dz,Dz)				\
  private(ix)					\
  schedule(static,chunk)			  
  for(ix = 0;ix<original_density.Nx();ix++)
    for(int iy = 0;iy<original_density.Ny();iy++)
      {
	// set up source term for tridiagonal equations
	DFT_Vec RHS(Nz);

	for(int iz=0;iz<Nz;iz++)
	  {
	    long i0 = original_density.get_PBC_Pos(ix,iy,iz);

	    long ipz = original_density.get_PBC_Pos(ix,iy,iz+1);
	    long imz = original_density.get_PBC_Pos(ix,iy,iz-1);

	    double r0 = original_density.getDensity(i0);	      		
	    double rpz = original_density.getDensity(ipz);
	    double rmz = original_density.getDensity(imz);

	    double RHS_F = 0.0;
	    
	    if(bUseForces)
	      {	   
		// We need the full gradient
		long ipx = original_density.get_PBC_Pos(ix+1,iy,iz);
		long imx = original_density.get_PBC_Pos(ix-1,iy,iz);
		
		long ipy = original_density.get_PBC_Pos(ix,iy+1,iz);
		long imy = original_density.get_PBC_Pos(ix,iy-1,iz);
	  				
		// these are the new densities & forces
		double s0 = ynew.get(i0);

		double spx = ynew.get(ipx);
		double smx = ynew.get(imx);

		double spy = ynew.get(ipy);
		double smy = ynew.get(imy);
	    
		double spz = ynew.get(ipz);
		double smz = ynew.get(imz);

		double g0 = dF_.get(i0);

		double gpx = dF_.get(ipx);
		double gmx = dF_.get(imx);

		double gpy = dF_.get(ipy);
		double gmy = dF_.get(imy);

		double gpz = dF_.get(ipz);
		double gmz = dF_.get(imz);

		RHS_F = Dx*((spx+s0)*(gpx-g0)-(s0+smx)*(g0-gmx))
		  +Dy*((spy+s0)*(gpy-g0)-(s0+smy)*(g0-gmy))
		  +Dz*((spz+s0)*(gpz-g0)-(s0+smz)*(g0-gmz));

		// now the old densities and forces

		double rpx = original_density.getDensity(ipx);
		double rmx = original_density.getDensity(imx);

		double rpy = original_density.getDensity(ipy);
		double rmy = original_density.getDensity(imy);

		double f0 = oldF_.get(i0);

		double fpx = oldF_.get(ipx);
		double fmx = oldF_.get(imx);

		double fpy = oldF_.get(ipy);
		double fmy = oldF_.get(imy);

		double fpz = oldF_.get(ipz);
		double fmz = oldF_.get(imz);

		RHS_F -= Dx*((rpx+r0)*(fpx-f0)-(r0+rmx)*(f0-fmx))
		  +Dy*((rpy+r0)*(fpy-f0)-(r0+rmy)*(f0-fmy))
		  +Dz*((rpz+r0)*(fpz-f0)-(r0+rmz)*(f0-fmz));

		// factor 1/2 because of density average used in prev line
		// factor dV because f0, etc carries a dV
		RHS_F *= 1.0/(2*dV);
	      }
	    RHS_F = 0.0;

	    //	    RHS.set(iz, ynew.get(i0) - Dz*(rpz+rmz-2*r0));
	    RHS.set(iz, ynew.get(i0) - 0.5*Dz*(rpz+rmz-2*r0)+ 0.5*dt_*RHS_F);

	    if(bFixedBoundaries)
	      if(iz == 0 || iz == Nz-1) RHS.set(iz, r0);	    	    
	  }

	if(!bFixedBoundaries)
	  solv_periodic_tridiag(RHS,0.5*Dz);
	else {
	  DFT_Vec b(original_density.Nz());
	  b.set(0,1.0);
	  b.set(original_density.Nz()-1,1.0);
	  for(int kk=1;kk<original_density.Nz()-1;kk++)
	    b.set(kk,1+2*(0.5*Dz));
	  solv_tridiag(b,RHS,0.5*Dz,bFixedBoundaries);
	}

	for(int iz=0;iz<original_density.Nz() && !bLessThanZero;iz++)
	  {
	    long i0 = original_density.get_PBC_Pos(ix,iy,iz);	
	    double val = 0; // = RHS.get(iz);
	    if(bFixedBoundaries && (ix == 0 || ix == original_density.Nx()-1 || iy == 0 || iy == original_density.Ny()-1))
	      val = original_density.getDensity(i0);
	    else val = RHS.get(iz);

	    if(val < SMALL_VALUE) val = SMALL_VALUE;
	    
	    deviation += fabs(val-ynew.get(i0));
	    ynew.set(i0,val);
	    
	    //	    if(val < 0) {bLessThanZero = true; cout << "density ==> " << val << " ix = " << ix << " iy = " << iy << " iz = " << iz << endl;}
	  }
      }
  return bLessThanZero;
}
/*
  double DDFT::step_string(double &dt, Density &original_density, double self_consistency_threshold)
  {
    dt_ = dt;

    bool bSelfConsistent = (self_consistency_threshold > 0);

    double F = dft_.calculateFreeEnergyAndDerivatives(original_density,0.0, dF_,true);

    DFT_Vec reset(original_density.getDensity());    

    if(bSelfConsistent) oldF_.set(dF_);
    
    bool bLessThanZero = false;
    do {
      DFT_Vec new_density; 
      new_density.zeros(original_density.Ntot());
    
      sub_step_x(new_density,original_density);
      sub_step_y(new_density,original_density);

      // If there is no force calculation, then (dF_ - original_Force) = 0 and does not contribute ....
      if(!bSelfConsistent)          
	bLessThanZero = sub_step_z(new_density,original_density);
      else {
	DFT_Vec pre_update_density(new_density); // only used to judge the updating of the density

	for(int iteration = 0; iteration < 100 && deviation > self_consistency_threshold; iteration++)
	  {
	    // calc new force - make sure to restore the old density
	    DFT_Vec temp(original_density.getDensity());
	    original_density.set(new_density);
	    dft_.calculateFreeEnergyAndDerivatives(original_density,0.0, dF_,true);
	    original_density.set(temp);
	  
	    deviation = 0;
	    bLessThanZero = sub_step_z(new_density,original_density, bSelfConsistent);
	    if(bLessThanZero) break;

	    deviation /= (new_density.size()*dt_);
	  
	    cout << "\titeration = " << iteration << " diff = " << deviation << endl;
	  
	  }
	if(!bLessThanZero) // if we think we have a good solution ...
	  if(deviation > self_consistency_threshold)
	    throw std::runtime_error("Could not meet self-consistency condition ... ");
      
      }
      if(bLessThanZero) // one of the updates failed
	{
	  dt_ = (dt /= 2);
	  original_density.set(reset); // this isn't needed, is it?
	  if(bSelfConsistent) dF_.set(oldF_);	// only needed when doing self-consistent calc
	} else original_density.set(new_density);    
    } while(bLessThanZero);
  
    F = dft_.calculateFreeEnergyAndDerivatives(original_density,0.0, dF_,true);
    DFT_Vec dummy;
    F += dft_.F_IdealGas(original_density,dummy);
    F += dft_.F_External(original_density,0.0,dummy);

    return F;
  }
*/

double DDFT::step_string(double &dt, Density &original_density, double self_consistency_threshold)
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
  /*
  // Adaptive time-step: try to increase time step if the present one works 5 times 
  if(decreased_time_step) successes_ = 0;
  else successes_++;
  if(successes_ >= 5 && dt_ < dtMax_) { dt_ = min(2*dt, dtMax_); successes_ = 0;}
  */
  
  original_density.set(d1);
  original_density.doFFT();
  calls_++;

  //  F_ = dft_.calculateFreeEnergyAndDerivatives(original_density,0.0, dF_,false);  
    
  dt = dt_;

  return F_;
}



  void DDFT::reverseForce(DFT_Vec *tangent) //;reverseForce(DFT_Vec *tangent)
  {
    double prod = dF_.dotWith(*tangent);
    dF_.Increment_And_Scale(*tangent,-2*prod);
  }


  // Solve a tridiagonal system of the form
  //       b[0] -D    0    0  0  0  0 ...  0    0    = RHS[0]
  //       -D   b[1] -D    0  0  0  0 ...  0    0    = RHS[1]
  //        0   -D   b[2] -D  0  0  0 ...  0    0    = RHS[2]
  //                  ...
  //        0   0     0    0  0  0  0 ... -D  b[N-1] = RHS[N-1]
  //
  //     The upper band is c[i] so that initially c[i] = -D for 0 <= i <= N-2
  //     The lower band is a[i]
  //     In the first pass, we eliminate the lower band
  //              1    c[i-1] 0   ... = RHS[i-1]
  //            a[i]   b[i]  c[i] ... = RHS[i]
  //
  //  If bFixedBoundaries = true, c[0] ==> 0 and a[N-1] ==> 0
void DDFT::solv_tridiag(const DFT_Vec &b, DFT_Vec &RHS, double D, bool bFixedBoundaries)
  {
    int N = RHS.size();
  
    DFT_Vec a(N); // lower
    DFT_Vec c(N); // upper

    for(int i=0;i<N;i++) {a.set(i,-D); c.set(i,-D);}

    if(bFixedBoundaries)
      {c.set(0,0.0); a.set(N-1,0.0);}
    
    // Scale the first line so that b[0] ==> 1
    c.scaleBy(0,b.get(0)); RHS.scaleBy(0,b.get(0));
    for(int i=1;i<N;i++)
      {
	double v = b.get(i)-a.get(i)*c.get(i-1);
	c.scaleBy(i,v);
	RHS.set(i,(RHS.get(i) - a.get(i)*RHS.get(i-1))/v);
      }
    for(int i=N-2;i>=0;i--)
      RHS.set(i, RHS.get(i) - c.get(i)*RHS.get(i+1));
  }

  // Solve a tridiagonal system of the form
  //       b[0]  -D    0    0  0  0  0 ...  0   -D    = RHS[0]
  //       -D    b[1] -D    0  0  0  0 ...  0    0    = RHS[1]
  //        0    -D   b[2] -D  0  0  0 ...  0    0    = RHS[2]
  //                  ...
  //       -D    0     0    0  0  0  0 ... -D  b[N-1] = RHS[N-1]
  //
  //     The upper band is c[i] so that initially c[i] = -D for 0 <= i <= N-2
  //     The lower band is a[i]
  //     In the first pass, we eliminate the lower band
  //              1    c[i-1] 0   ...  0  beta[i-1] = RHS[i-1]
  //            a[i]   b[i]  c[i] ...  0     0      = RHS[i]
  //                      ...
  //             alpha  0    0    ... -D     b[N-1] = RHS[N-1]       
  //
  
void DDFT::solv_periodic_tridiag_2(DFT_Vec &b, DFT_Vec &RHS, double D)
{
  int N = RHS.size();

  double a = -D; // lower band - these values are all the same and never change
  double alpha = -D; // M[N-1][j]
  DFT_Vec c(N); // upper
  DFT_Vec beta(N); // right-most column

  for(int i=0;i<N;i++) c.set(i,-D);
  beta.zeros(N);
  beta.set(0,-D);
  
  // Scale the first line so that b[0] ==> 1
  c.scaleBy(0,b.get(0)); beta.scaleBy(0,b.get(0)); RHS.scaleBy(0,b.get(0));
  for(int i=1;i<N;i++)
    {
      double v = b.get(i)-a*c.get(i-1);

      RHS.set(i,(RHS.get(i) - a*RHS.get(i-1))/v);

      if(i < N-2)
	{
	  c.scaleBy(i,v);
	  beta.set(i,-a*beta.get(i-1)/v);
	} else if(i == N-2) {
	c.set(  N-2, (  c.get(N-2) - a*beta.get(N-3))/v);
	beta.set(N-2, 0);
      }

      // also get rid of alpha at position (N-1,i-1) and create a new one at (N-1,i)
      if(i < N-1)
	{
	  RHS.set(N-1, RHS.get(N-1) -alpha* RHS.get(i-1));
	  b.set(  N-1,   b.get(N-1) -alpha*beta.get(i-1));
	  alpha =                   -alpha*   c.get(i-1);
	}
      if(i == N-2) a += alpha;
    }
  for(int i=N-2;i>=0;i--)
    RHS.set(i, RHS.get(i) - c.get(i)*RHS.get(i+1) - beta.get(i)*RHS.get(N-1));
  cout << "b[N-1] = " << b.get(N-1)-a*c.get(N-2) << " RHS[N-1] = " << RHS.get(N-1) << endl;
}

  

  void DDFT::solv_periodic_tridiag(DFT_Vec &RHS, double D)
  {
    int N = RHS.size();
  
    // modify entries
    double alpha = -D;
    double beta  = -D;
    double gamma = -1-2*D;

    DFT_Vec b(N); // diag

    for(int i=0;i<N;i++) b.set(i,1+2*D);

    b.addTo(0,-gamma);
    b.addTo(N-1,-alpha*beta/gamma);

    solv_tridiag(b,RHS,D);

    DFT_Vec u(N); u.zeros(N);

    u.set(0,gamma);
    u.set(N-1,alpha);

    solv_tridiag(b,u,D);

    double fact = (RHS.get(0)+beta*RHS.get(N-1)/gamma)/(1+u.get(0)+beta*u.get(N-1)/gamma);

    for(int i=0;i<N;i++) RHS.addTo(i, -fact*u.get(i)); 
  }


  int cc = 0;

double DDFT::step()
{
  bool bLessThanZero;
  cout << "===================================================================================================================" << endl;
  /*
    do {
    DFT_Vec new_density(density_.Ntot());

    sub_step_x(new_density,density_,bFixedBoundaries_);
    sub_step_y(new_density,density_,bFixedBoundaries_);
    bLessThanZero = sub_step_z(new_density,density_, false,bFixedBoundaries_);

    if(bLessThanZero) {dt_ /= 2; cout << "dt_ = " << dt_ << endl;} // density went negative: reduce time step
    else density_.set(new_density);    
    } while(bLessThanZero);
  */

  //    dt_ = 1e-4;
    

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
    if(deviation > tolerence_fixed_point_)
      throw std::runtime_error("Failed to find fixed point in DDFT::step()");
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

  return F_;
}

/**
 This function takes the input density and calculates a new density, d1, by propagating the density a time step dt. The nonlinear terms, RHS0 and RHS1, are treated explicitly.
*/
double DDFT::fftDiffusion(DFT_Vec &d1, const DFT_FFT &RHS0, const DFT_FFT &RHS1)
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


void DDFT::calcNonlinearTerm(const DFT_Vec &d2, const DFT_Vec &dF, DFT_Vec &RHS1)
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



double DDFT::F_string(Density &original_density, double *fmax)
{
  DFT_Vec dummy;
  double F =  dft_.calculateFreeEnergyAndDerivatives(original_density,0.0, dummy,false);

  if(fmax) *fmax = dummy.inf_norm()/original_density.dV();;
  
  return F;
}



void DDFT::test_solv_tridiag()
{
  int N = 8;

  double rmin = -10;
  double rmax =  10;
  
  double D = 0.8;
   
  DFT_Vec y(N);
  DFT_Vec rhs(N);
  DFT_Vec b(N);

  for(int i=0;i<N;i++)
    rhs.set(i,rmin + (rmax-rmin)*rand()*1.0/RAND_MAX);

  y.set(rhs);
  
  for(int i=0;i<N;i++)
      b.set(i,rmin + (rmax-rmin)*rand()*1.0/RAND_MAX);
  
  double M[N][N];
  for(int i=0;i<N;i++)
    for(int j=0;j<N;j++)
      M[i][j] = 0.0;

  for(int i=0;i<N;i++)
    {
      M[i][i] = b.get(i);
      if(i < N-1) M[i][i+1] = -D;
      if(i > 0)   M[i][i-1] = -D;
    }

  cout << "M:" << endl;
  for(int i=0;i<N;i++)
    {
      for(int j=0;j<N;j++)
	cout << M[i][j] << "\t";
      cout << endl;
    }
  cout << endl;

  cout << "Test solution with free boundaries:" << endl;
  
  bool bFixedBoundaries = false;
  solv_tridiag(b, rhs, D, bFixedBoundaries);
  
  for(int i=0;i<N;i++)
    {
      double x = 0;
      for(int j=0;j<N;j++)
	x += M[i][j]*rhs.get(j);
      cout << i << " : y = " << rhs.get(i) << " M*y = " << x << " RHS = " << y.get(i) << " diff = " << x-y.get(i) << endl;
    }
  cout << endl;

  cout << "Test periodic boundries:" << endl;

  M[0][N-1] = M[N-1][0] = -D;
  rhs.set(y);
  
  solv_periodic_tridiag(rhs, D);
  
  for(int i=0;i<N;i++)
    {
      double x = 0;
      for(int j=0;j<N;j++)
	x += (i == j ? 1+2*D : M[i][j])*rhs.get(j);
      cout << i << " : y = " << rhs.get(i) << " M*y = " << x << " RHS = " << y.get(i) << " diff = " << x-y.get(i) << endl;
    }
  cout << "Test periodic boundries alg2:" << endl;

  M[0][N-1] = M[N-1][0] = -D;
  rhs.set(y);

  for(int i=0;i<N;i++)
    {
      b.set(i,1+2*D); 
      M[i][i] = 1+2*D;
    }
  solv_periodic_tridiag_2(b,rhs, D);
  
  for(int i=0;i<N;i++)
    {
      double x = 0;
      for(int j=0;j<N;j++)
	x += M[i][j]*rhs.get(j);
      cout << i << " : y = " << rhs.get(i) << " M*y = " << x << " RHS = " << y.get(i) << " diff = " << x-y.get(i) << endl;
    }

  cout << "Test solution with fixed boundaries:" << endl;

  M[0][N-1] = M[N-1][0] = 0;
  
  M[0][0] = 1; M[0][1] = 0.0;
  M[N-1][N-1] = 1; M[N-1][N-2] = 0.0;
    
  for(int i=0;i<N;i++)
    b.set(i,M[i][i]);
  
  rhs.set(y);
  
  bFixedBoundaries = true;
  solv_tridiag(b, rhs, D, bFixedBoundaries);
  
  for(int i=0;i<N;i++)
    {
      double x = 0;
      for(int j=0;j<N;j++)
	x += M[i][j]*rhs.get(j);
      cout << i << " : y = " << rhs.get(i) << " M*y = " << x << " RHS = " << y.get(i) << " diff = " << x-y.get(i) << endl;
    }

  
}

void DDFT::Display(double F, double dFmin, double dFmax, double N)
{
  static int cc = 0;
  
  grace_->deleteDataSet(0);
  grace_->deleteDataSet(1);
  grace_->deleteDataSet(2);
  for(int i=0;i<density_.Nx();i++)
    grace_->addPoint(i,density_.getDensity(i,density_.Ny()/2, density_.Nz()/2),0);
  
  for(int i=0;i<density_.Ny();i++)
    grace_->addPoint(i,density_.getDensity(density_.Nx()/2,i, density_.Nz()/2),1);
  
  for(int i=0;i<density_.Nz();i++)
    grace_->addPoint(i,density_.getDensity(density_.Nx()/2,density_.Ny()/2,i),2);
  
  if(cc == 0)
    for(int i=0;i<density_.Nz();i++)
      grace_->addPoint(i,density_.getDensity(density_.Nx()/2,density_.Ny()/2, i),3);
  
  cc = 1;
  
  grace_->redraw();
  
  stringstream ss;
  ss << "Step = " << step_counter_ << " F = " << F << " dFMin = " << dFmin << " dFmax = " << dFmax << " N = " << N;
  cout << "Setting title to: " << ss.str() << endl;
  grace_->setTitle(ss.str().c_str());
  
  grace_->redraw(1,0);
  string s("string_graph.agr");
  grace_->store(s);
}


