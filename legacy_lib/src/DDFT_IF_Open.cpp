#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

//#include <mgl2/mgl.h>

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

  //  if(species->is_background_fixed() != true)
  //    throw std::runtime_error("DDFT_IF_OPEN requires that the species have fixed background set to true");
  
  RHS0.zeros(species->getDensity().Ntot());
  RHS1.zeros(species->getDensity().Ntot());

  sin_Ntot_ = (Nx_-1)*(Ny_-1)*(Nz_-1);
  sin_Norm_ = 8*Nx_*Ny_*Nz_;
  
  RHS0_sin_transform_ = new double[sin_Ntot_];
  RHS1_sin_transform_ = new double[sin_Ntot_];
    

  sin_in_  = new double[sin_Ntot_];
  sin_out_ = new double[sin_Ntot_];
  
  sin_plan_ = fftw_plan_r2r_3d(Nx_-1, Ny_-1, Nz_-1,
			       sin_in_, sin_out_,
			       FFTW_RODFT00, FFTW_RODFT00, FFTW_RODFT00,
			       FFTW_MEASURE);			       

  pack_for_sin_transform(species->get_density_data(),background_);
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


void DDFT_IF_Open::calcNonlinearTerm(const DFT_Vec &density, Species *species, bool use_R0)
{
  if(use_R0)
    {
      DDFT_IF::calcNonlinearTerm_intern(density, species->getDF(),RHS0);
      pack_for_sin_transform(RHS0.memptr(),0.0);
      fftw_execute(sin_plan_);
      memcpy(RHS0_sin_transform_,sin_out_,sin_Ntot_*sizeof(double));      
    } else {
    DDFT_IF::calcNonlinearTerm_intern(density, species->getDF(),RHS1); // we are counting on the species to set the force on the border points to zero
	pack_for_sin_transform(RHS1.memptr(),0.0);
	fftw_execute(sin_plan_);
	memcpy(RHS1_sin_transform_,sin_out_,sin_Ntot_*sizeof(double));	  
  }
}


/*
double DDFT_IF_Open::step()
{
  DDFT_IF::step();

  int Jspecies = 0;
  Species *species = dft_->getSpecies(Jspecies);
  const Density &density = species->getDensity();
  
  DFT_Vec d0(density.Ntot()); d0.set(density.getDensity());  
  DFT_Vec d1(density.Ntot()); d1.set(d0);

  const bool USE_R0 = true;
  
  calcNonlinearTerm(d0,species,USE_R0);
   
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
	F_ = dft_->calculateFreeEnergyAndDerivatives(false);
	calcNonlinearTerm(d1,species,!USE_R0);	

	species->set_density(d0); // Unnecessary?
	species->fft_density(); // Unnecessary?

	deviation = fftDiffusion(d1);
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

  return F_;
}
*/

void DDFT_IF_Open::finish_nonlinear_calc(DFT_Vec& d0, DFT_Vec& d1)
{
  // nothing to do in this case ...
}



/**
 This function takes the input density and calculates a new density, d1, by propagating the density a time step dt. The nonlinear terms, RHS0 and RHS1, are treated explicitly.
*/
double DDFT_IF_Open::fftDiffusion(DFT_Vec &d1)
{
  int Jspecies = 0;
  Species *species = dft_->getSpecies(Jspecies);
  const Lattice &lattice = species->getLattice();
  
  pack_for_sin_transform(species->get_density_data(),background_);
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
	  x += ((fac-1)/Lambda)*RHS0_sin_transform_[pos];
	  x += ((fac-1-dt_*Lambda)/(Lambda*Lambda*dt_))*(RHS1_sin_transform_[pos]-RHS0_sin_transform_[pos]);

	  sin_in_[pos] = x;
	  pos++;
	}

  fftw_execute(sin_plan_);

  double deviation = 0;
  double maxdeviation = 0;

  // Notice that sin_out is indexed using sign_out[1] as element 0 ... this differs from the closed case
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







