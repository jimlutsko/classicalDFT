#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

#include <mgl2/mgl.h>

#ifdef USE_OMP
#include <omp.h>
#endif



using namespace std;

// required MPI include file  
//#include <mpi.h>

#include "Minimizer.h"

void DDFT_IF_Open::initialize()
{
  DDFT_IF::initialize();

  int Jspecies = 0;
  Species *species = dft_->getSpecies(Jspecies);  
  
  sin_Ntot_ = (Nx_-1)*(Ny_-1)*(Nz_-1);
  sin_Norm_ = 8*Nx_*Ny_*Nz_;

  sin_in_  = new double[sin_Ntot_];
  sin_out_ = new double[sin_Ntot_];
  
  sin_plan_ = fftw_plan_r2r_3d(Nx_-1, Ny_-1, Nz_-1,
			       sin_in_, sin_out_,
			       FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00,
			       FFTW_MEASURE);			       

  pack_for_sin_transform(species->getDensityData(),background_);
  fftw_execute(sin_plan_);
  memcpy(sin_in_, sin_out_, sin_Ntot_*sizeof(double));
  fftw_execute(sin_plan_);
  for(unsigned k=0;k<sin_Ntot_;k++) sin_out_[k] /= sin_Norm_;

  double Dx = 1.0/(dx_*dx_);
  double Dy = 1.0/(dy_*dy_);
  double Dz = 1.0/(dz_*dz_);
  
  for(int ix=0;ix<Nx_-1;ix++)
    {
      double kx = (M_PI*(ix+1))/Nx_;
      double facx = 2*Dx*(cos(kx)-1);

      Lamx_.push_back(facx);
    }
      
  for(int iy=0;iy<Ny_-1;iy++)
    {
      double ky = (M_PI*(iy+1))/Ny_;
      double facy = 2*Dy*(cos(ky)-1);

      Lamy_.push_back(facy);
    }

  for(int iz=0;iz<Nz_-1;iz++)
    {
      double kz = (M_PI*(iz+1))/Nz_;
      double facz = 2*Dz*(cos(kz)-1);

      Lamz_.push_back(facz);
    }
}

double DDFT_IF_Open::step()
{
  DDFT_IF::step();


  int Jspecies = 0;
  Species *species = dft_->getSpecies(Jspecies);
  const Density &density = species->getDensity();
  
  cout << "Initial F = " << F_ << " x = " << density.get(0,0,0) << endl;

  DFT_Vec d0(density.Ntot()); d0.set(density.getDensity());  
  DFT_Vec RHS0(density.Ntot());
  calcNonlinearTerm(d0, species->getDF(), RHS0);
  double *RHS0_sin_transform = new double[sin_Ntot_];

  pack_for_sin_transform(RHS0.memptr(),0.0);
  fftw_execute(sin_plan_);
  memcpy(RHS0_sin_transform,sin_out_,sin_Ntot_*sizeof(double));
    
  DFT_Vec d1(density.Ntot()); d1.set(d0);
  DFT_Vec RHS1; RHS1.zeros(density.Ntot());  
  double *RHS1_sin_transform = new double[sin_Ntot_];
  memcpy(RHS1_sin_transform,RHS0_sin_transform,sin_Ntot_*sizeof(double));

  double deviation = 1;

  double dt = dt_;
  
  bool reStart;
  bool decreased_time_step = false;
  
  do {
    reStart = false;
    double old_error = 0;

    for(int i=0;i<20 && deviation > tolerence_fixed_point_ && !reStart;i++)
      {
	species->set_density(d1);
	F_ = dft_->calculateFreeEnergyAndDerivatives(false); //density_,0.0, dF_,false);
	calcNonlinearTerm(d1, species->getDF(),RHS1);
	pack_for_sin_transform(RHS1.memptr(),0.0);
	fftw_execute(sin_plan_);
	memcpy(RHS1_sin_transform,sin_out_,sin_Ntot_*sizeof(double));	  

	species->set_density(d0); // Unnecessary?
	species->fft_density(); // Unnecessary?

	deviation = fftDiffusion(d1,RHS0_sin_transform,RHS1_sin_transform);
	cout << "\tdeviation = " << deviation << " dt = " << dt_ << endl;

	// decrease time step and restart if density goes negative or if error is larger than previous step
	if(d1.min() < 0 || (i > 0 && old_error < deviation))
	  {reStart = true; dt_ /= 10; d1.set(d0); decreased_time_step = true;}

	old_error = deviation;	
      }
    if(!reStart && deviation > tolerence_fixed_point_)
      {reStart = true; dt_ /= 10; d1.set(d0); decreased_time_step = true;}
  } while(reStart);

  if(decreased_time_step) successes_ = 0;
  else successes_++;
  if(successes_ >= 5 && dt_ < dtMax_) { dt_ = min(2*dt, dtMax_); successes_ = 0;}
  
  species->set_density(d1); 
  species->fft_density();
  calls_++;

  delete RHS0_sin_transform;
  delete RHS1_sin_transform;

  return F_;
}

/**
 This function takes the input density and calculates a new density, d1, by propagating the density a time step dt. The nonlinear terms, RHS0 and RHS1, are treated explicitly.
*/
double DDFT_IF_Open::fftDiffusion(DFT_Vec &d1, const double *RHS0_sin_transform, const double *RHS1_sin_transform)
{
  int Jspecies = 0;
  Species *species = dft_->getSpecies(Jspecies);
  const Lattice &lattice = species->getLattice();
  
  pack_for_sin_transform(species->getDensityData(),background_);
  fftw_execute(sin_plan_);
  
  // fft of density is now in sin_out_;
  // we construct the rhs in sin_in_ since we will fft it to get the density in real space

  vector<double> fx;
  for(int ix=0;ix<Lamx_.size();ix++)
    fx.push_back(exp(Lamx_[ix]*dt_));

  vector<double> fy;
  for(int iy=0;iy<Lamy_.size();iy++)
    fy.push_back(exp(Lamy_[iy]*dt_));
  
  vector<double> fz;
  for(int iz=0;iz<Lamz_.size();iz++)
    fz.push_back(exp(Lamz_[iz]*dt_));
  
  unsigned pos = 0;
  for(int ix=0;ix<Lamx_.size();ix++)
    for(int iy=0;iy<Lamy_.size();iy++)
      for(int iz=0;iz<Lamz_.size();iz++)
	{
	  double Lambda = Lamx_[ix]+Lamy_[iy]+Lamz_[iz];
	  double fac    = fx[ix]*fy[iy]*fz[iz];
	  
	  double x = sin_out_[pos];
	  
	  x *= fac;
	  x += ((fac-1)/Lambda)*RHS0_sin_transform[pos];
	  x += ((fac-1-dt_*Lambda)/(Lambda*Lambda*dt_))*(RHS1_sin_transform[pos]-RHS0_sin_transform[pos]);

	  sin_in_[pos] = x;
	  pos++;
	}

  fftw_execute(sin_plan_);

  double deviation = 0;
  double maxdeviation = 0;
  
  pos = 0;
  for(int ix=1;ix<Nx_;ix++)
    for(int iy=1;iy<Ny_;iy++)
      for(int iz=1;iz<Nz_;iz++)
	{
	  double d     = (sin_out_[pos++] /= sin_Norm_)+background_;
	  double d1val = d1.get(lattice.pos(ix,iy,iz));
	  double u     = fabs(d1val-d);
	  
	  deviation += u*u;
	  if(u > maxdeviation) maxdeviation = u;	  
	}
  unpack_after_transform(d1.memptr(), background_);
  return maxdeviation; 
}

void DDFT_IF_Open::pack_for_sin_transform(const double *x, double val)
{
  int Jspecies = 0;
  const Density &density = dft_->getSpecies(Jspecies)->getDensity();
  
  unsigned pos = 0;
  for(int ix=1;ix<Nx_;ix++)
    for(int iy=1;iy<Ny_;iy++)
      for(int iz=1;iz<Nz_;iz++)
	  sin_in_[pos++] = x[density.pos(ix,iy,iz)] - val;
}

void DDFT_IF_Open::unpack_after_transform(double *x, double val)
{
  int Jspecies = 0;
  const Density &density = dft_->getSpecies(Jspecies)->getDensity();
  
  unsigned pos = 0;
  for(int ix=0;ix<Nx_;ix++)
    for(int iy=0;iy<Ny_;iy++)
      for(int iz=0;iz<Nz_;iz++)
	{
	  if(ix == 0 || iy == 0 || iz == 0)
	    x[density.pos(ix,iy,iz)] = val;
	  else
	    {
	      double v = sin_out_[pos++] + val;
	      x[density.pos(ix,iy,iz)] = v;
	    }
	}
}

void DDFT_IF_Open::calcNonlinearTerm(const DFT_Vec &d2, const DFT_Vec &dF, DFT_Vec &RHS1)
{
  int Jspecies = 0;
  const Density &density = dft_->getSpecies(Jspecies)->getDensity();
  
  double Dx = 1.0/(dx_*dx_);
  double Dy = 1.0/(dy_*dy_);
  double Dz = 1.0/(dz_*dz_);

  double dV = dx_*dy_*dz_;

  int iy;

#pragma omp parallel for  shared(dV)  private(iy) schedule(static)			  
  for(iy = 0;iy<Ny_;iy++)
    for(int iz = 0;iz<Nz_;iz++)
      for(int ix=0;ix<Nx_;ix++)
	{
	  long i0 = density.get_PBC_Pos(ix,iy,iz);

	  long ipx = density.get_PBC_Pos(ix+1,iy,iz);
	  long imx = density.get_PBC_Pos(ix-1,iy,iz);

	  long ipy = density.get_PBC_Pos(ix,iy+1,iz);
	  long imy = density.get_PBC_Pos(ix,iy-1,iz);

	  long ipz = density.get_PBC_Pos(ix,iy,iz+1);
	  long imz = density.get_PBC_Pos(ix,iy,iz-1);
	    
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


/*
double DDFT_IF_Open::step_string(double &dt, Density &original_density, unsigned &time_den, bool verbose)
{
  int Nx = original_density.Nx();
  int Ny = original_density.Ny();
  int Nz = original_density.Nz();

  long Ntot = original_density.Ntot();
  
  dt_ = dt;

  F_ = calculateFreeEnergyAndDerivatives(original_density,0.0, dF_,false);

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
  memcpy(RHS1_sin_transform,RHS0_sin_transform,sin_Ntot_*sizeof(double));

  double deviation = 1;

  bool reStart;
  
  do {
    reStart = false;
    double old_error = 0;

    for(int i=0;i<20 && deviation > tolerence_fixed_point_ && !reStart;i++)
      {
	density_.set(d0); // Unnecessary?
	density_.doFFT(); // Unnecessary?
	deviation = fftDiffusion(d1,RHS0_sin_transform,RHS1_sin_transform);
	if(verbose) cout << "\tdeviation = " << deviation << " dt = " << dt_ << endl;

	// decrease time step and restart if density goes negative or if error is larger than previous step
	if(d1.min() < 0 || (i > 0 && old_error < deviation))
	    reStart = true;
	else { // prepare for next step
	  density_.set(d1);
	  F_ = dft_.calculateFreeEnergyAndDerivatives(density_,0.0, dF_,false);
	  calcNonlinearTerm(d1, dF_,RHS1);
	  pack_for_sin_transform(RHS1.memptr(),0.0);
	  fftw_execute(sin_plan_);
	  memcpy(RHS1_sin_transform,sin_out_,sin_Ntot_*sizeof(double));	  
	}	
	old_error = deviation;
      }
    if(deviation > tolerence_fixed_point_)
      reStart = true;

    if(reStart)
      {
	dt_ /= 2;
	time_den *= 2;

	d1.set(d0);
	memcpy(RHS1_sin_transform,RHS0_sin_transform,sin_Ntot_*sizeof(double));
      }
    
  } while(reStart);
    
  original_density.set(d1); 
  original_density.doFFT();// Unnecessary?
  calls_++;

  dt = dt_;

  delete RHS0_sin_transform;
  delete RHS1_sin_transform;

  return F_;
}
*/







