#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <time.h>

//#include <mgl2/mgl.h>
//#include <mgl2/fltk.h>

#include <complex>

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


DDFT_IF_Periodic::DDFT_IF_Periodic(DFT *dft, bool showGraphics)
  : DDFT_IF(dft,showGraphics)
{
  RHS0.initialize(Nx_,Ny_,Nz_);
  RHS1.initialize(Nx_,Ny_,Nz_);
  
  double Dx = 1.0/(dx_*dx_);
  double Dy = 1.0/(dy_*dy_);
  double Dz = 1.0/(dz_*dz_);

  for(int ix=0;ix<=Nx_-1;ix++)
    {
      double kx = (2*M_PI*ix)/Nx_;
      double facx = 2*Dx*(cos(kx)-1);
      
      Lamx_.push_back(facx);
    }
      
  for(int iy=0;iy<=Ny_-1;iy++)
    {
      double ky = (2*M_PI*iy)/Ny_;
      double facy = 2*Dy*(cos(ky)-1);

      Lamy_.push_back(facy);
    }

  for(int iz=0;iz<Nz_/2+1;iz++) // this is not a mistake
    {
      double kz = (2*M_PI*iz)/Nz_;
      double facz = 2*Dz*(cos(kz)-1);
      
      Lamz_.push_back(facz);
    }    
  
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

void DDFT_IF_Periodic::calcNonlinearTerm(const DFT_Vec &density, Species *species, bool use_R0)
{
  if(use_R0)
    {
      DDFT_IF::calcNonlinearTerm_intern(density, species->getDF(), RHS0.Real());    
      RHS0.do_real_2_fourier();
    } else {
    DDFT_IF::calcNonlinearTerm_intern(density, species->getDF(),RHS1.Real());  
    RHS1.do_real_2_fourier();
  }
}

/**
 This function takes the input density and calculates a new density, d1, by propagating the density a time step dt. The nonlinear terms, RHS0 and RHS1, are treated explicitly.
*/
double DDFT_IF_Periodic::fftDiffusion(DFT_Vec &new_density)
{
  int Jspecies = 0;
  Species *species = dft_->getSpecies(Jspecies);
  const Density &density = species->getDensity();

  species->doFFT();
  
  // save some evaluations of the exponent
  vector<double> fx;
  for(int ix=0;ix<Lamx_.size();ix++)
    fx.push_back(exp(Lamx_[ix]*dt_));

  vector<double> fy;
  for(int iy=0;iy<Lamy_.size();iy++)
    fy.push_back(exp(Lamy_[iy]*dt_));
  
  vector<double> fz;
  for(int iz=0;iz<Lamz_.size();iz++)
    fz.push_back(exp(Lamz_[iz]*dt_));  

  DFT_FFT work(Nx_,Ny_,Nz_);
  DFT_Vec_Complex &cwork = work.Four();
  
  unsigned pos = 0;
  for(int ix = 0;ix<Lamx_.size();ix++)
    for(int iy=0;iy<Lamy_.size();iy++)
      for(int iz=0;iz<Lamz_.size();iz++)
	{
	  complex<double> x = density.get_fourier_value(pos); 
	  
	  if(pos > 0) // pos == 0 corresponds to K=0 - and mass is conserved so nothing to do ...
	    {
	      double Lambda = Lamx_[ix]+Lamy_[iy]+Lamz_[iz];
	      double exp_dt = fx[ix]*fy[iy]*fz[iz];
	      x *= exp_dt;
	      x += ((exp_dt-1)/Lambda)*RHS0.cFour().get(pos);
	      x += ((exp_dt-1-dt_*Lambda)/(Lambda*Lambda*dt_))*(RHS1.cFour().get(pos)-RHS0.cFour().get(pos));
	    }
	  cwork.set(pos,x);
	  pos++;	  
	}
  
  work.do_fourier_2_real();
  work.Real().MultBy(1.0/density.Ntot());

  double deviation = 0;
  double maxdeviation = 0;
    
  for(unsigned i=0;i<new_density.size();i++)
    {
      double d = new_density.get(i);
      double w = work.Real().get(i);
      double u = fabs(d-w);
      
      deviation += u*u;    
      if(u > maxdeviation) maxdeviation = u;
    }
  new_density.set(work.Real());
  return maxdeviation; 
}


/*
double DDFT_IF_Periodic::step_string(double &dt, Density &original_density, double self_consistency_threshold)
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

