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

  int Jspecies = 0;
  Species *species = dft_->getSpecies(Jspecies);  
  const Lattice &lattice = species->getLattice();
  
  Nx_ = lattice.Nx();
  Ny_ = lattice.Ny();
  Nz_ = lattice.Nz();

  dx_ = lattice.getDX();
  dy_ = lattice.getDY();
  dz_ = lattice.getDZ();

  double Dx = 1.0/(dx_*dx_);
  double Dy = 1.0/(dy_*dy_);
  double Dz = 1.0/(dz_*dz_);

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
  cout << "===================================================================================================================" << endl;


  if(dft_->getNumberOfSpecies() > 1)
    throw std::runtime_error("DDFT only implemented for single species systems ... aborting");
  
  F_ = 0;
  try {
    F_ = dft_->calculateFreeEnergyAndDerivatives(false);   
  } catch( Eta_Too_Large_Exception &e) {
    throw e;
  }
  cout << "Initial F = " << F_ << endl;

  return 0;
}


void DDFT_IF::calcNonlinearTerm(const DFT_Vec &d2, DFT_Vec &dF, DFT_Vec &RHS1)
{
  int Jspecies = 0;
  const Density &density = dft_->getSpecies(Jspecies)->getDensity();
  
  double Dx = 1.0/(dx_*dx_);
  double Dy = 1.0/(dy_*dy_);
  double Dz = 1.0/(dz_*dz_);

  double D[] = {Dx,Dy,Dz};

  double dV = dx_*dy_*dz_;

  if(fixedBorder_ && !isOPEN_) update_forces_fixed_background(density,d2,dF,D);
  cout << "fixedBorder_ = " << fixedBorder_ << " isOPEN_ = " << isOPEN_ << endl;

  A_dot_x(dF, RHS1, density, D, true);

}

double  static get(DFT_Vec &x, int ix, int iy, int iz, bool is_full, const Density& density)
{
  double r = 0;
  if(is_full)
    r = x.get(density.get_PBC_Pos(ix,iy,iz));
  else if(ix == 0 || iy == 0 || iz == 0)
    r = x.get(density.boundary_pos(ix,iy,iz));
  return r;
}

void DDFT_IF::A_dot_x(DFT_Vec& x, DFT_Vec& Ax, const Density &density, double D[], bool do_subtract_ideal) 
{
  const bool Ax_is_full = (Ax.size() == density.Ntot());
  const bool x_is_full  = ( x.size() == density.Ntot());

  const double dV = dx_*dy_*dz_;
  
  long pos;

  double sum_forces = 0.0;
  double sum_forces2 = 0.0;
  double sum_sq = 0.0;
    
  
#pragma omp parallel for  private(pos) schedule(static) reduction(+:sum_forces) reduction(+:sum_forces2) reduction(+:sum_sq)
  for(pos = 0;pos<Ax.size();pos++)
    {
      int ix, iy,iz; // working space

      if(Ax_is_full) density.cartesian(pos,ix,iy,iz);
      else density.boundary_cartesian(pos,ix,iy,iz);

      double x0  = get(x,ix,iy,iz,x_is_full, density);

      double xpx = get(x,ix+2,iy,iz,x_is_full, density);
      double xmx = get(x,ix-2,iy,iz,x_is_full, density);

      double xpy = get(x,ix,iy+2,iz,x_is_full, density);
      double xmy = get(x,ix,iy-2,iz,x_is_full, density);

      double xpz = get(x,ix,iy,iz+2,x_is_full, density);
      double xmz = get(x,ix,iy,iz-2,x_is_full, density);

      double r0  = density.get(ix,iy,iz);

      double rpx = density.get(ix+1,iy,iz);
      double rmx = density.get(ix-1,iy,iz);

      double rpy = density.get(ix,iy+1,iz);
      double rmy = density.get(ix,iy-1,iz);

      double rpz = density.get(ix,iy,iz+1);
      double rmz = density.get(ix,iy,iz-1);

      double RHS = D[0]*(rpx*xpx+rmx*xmx-(rpx+rmx)*x0)
	+D[1]*(rpy*xpy+rmy*xmy-(rpy+rmy)*x0)
	+D[2]*(rpz*xpz+rmz*xmz-(rpz+rmz)*x0);

      RHS /= (4*dV);

      if(do_subtract_ideal)
	{
	  if(ix != 0 && iy != 0 && iz != 0) sum_forces += RHS;
	  else sum_forces2 += RHS;

	  sum_sq += RHS*RHS;
	  
	  RHS -= D[0]*(rpx+rmx-2*r0)+D[1]*(rpy+rmy-2*r0)+D[2]*(rpz+rmz-2*r0) ;
	}
      Ax.set(pos,RHS);

    }
  if(do_subtract_ideal) cout << "sum_forces = " << sum_forces << " sum_forces2 = " << sum_forces2 << " rms = " << sqrt(sum_sq) << endl;
}

// This function corrects the forces on the border for the case of fixed background
// using conjuage gradients
void DDFT_IF::update_forces_fixed_background(const Density &density,const DFT_Vec &d2, DFT_Vec &dF, const double DD[3])
{
   long Nboundary = Nx_*Ny_+Nx_*Nz_+Ny_*Nz_-Nx_-Ny_-Nz_+1;
   double D[3] = {DD[0],DD[1],DD[2]};
  
  DFT_Vec r; r.zeros(Nboundary);   // residuals
  DFT_Vec p; p.zeros(Nboundary);   // CG intermediate quantity
  DFT_Vec x; x.zeros(Nboundary);   // Current approximation

  for(long pboundary = 0;pboundary<Nboundary;pboundary++)
    {
      int cartesian[3]; // working space      
      density.boundary_cartesian(pboundary,cartesian);
      x.set(pboundary, dF.get(density.get_PBC_Pos(cartesian)));
    }

  const bool boundary_only = true;

  A_dot_x(dF,r,density,D); 
  r.MultBy(-1);
  p.set(r);

  double error = sqrt(r.euclidean_norm()/r.size())/(dx_*dy_*dz_);
  double error_old = error;
  //  cout << "error = " << error << endl;

  while(error > 1e-13)
    {
      DFT_Vec Ap; Ap.zeros(Nboundary);

      A_dot_x(p,Ap,density,D);

      double r2  = r.dotWith(r);
      double pAp = p.dotWith(Ap);
      
      double alf = r2/pAp;
      r.IncrementBy_Scaled_Vector(Ap,-alf);
      x.IncrementBy_Scaled_Vector(p,alf);
      
      double beta = r.dotWith(r)/r2;
      p.MultBy(beta);
      p.IncrementBy(r);

      error = sqrt(r.euclidean_norm()/r.size())/(dx_*dy_*dz_);
      //      cout << "error = " << error << " pAp = " << pAp << " dError/error = : " << fabs(error-error_old)/error_old << endl;
      if(fabs(error-error_old) < 0.01*error_old) throw std::runtime_error("Convergence problem in DDFT_IF_CLOSED::update_forces_fixed_background");
      error_old = error;
    }

    for(long pboundary = 0;pboundary<Nboundary;pboundary++)
    {
      int cartesian[3]; // working space      
      density.boundary_cartesian(pboundary,cartesian);
      dF.set(density.get_PBC_Pos(cartesian), x.get(pboundary));
      if(pboundary == 64)
	cout << "HERE: " << cartesian[0] << " " << cartesian[1] << " " << cartesian[2] << " " << density.get_PBC_Pos(cartesian) << " " << density.get_PBC_Pos(cartesian) << " " << x.get(pboundary) << endl;
      
    }

}






void DDFT_IF::calcNonlinearTerm_old(const DFT_Vec &d2, const DFT_Vec &dF, DFT_Vec &RHS1)
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

	  RHS1.set(i0,RHS_F);
	}
}


