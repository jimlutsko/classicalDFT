#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

//#include <mgl2/mgl.h>
//#include <mgl2/fltk.h>

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

DDFT::DDFT(DFT *dft, bool showGraphics)
  : Minimizer(dft), show_(showGraphics) , tolerence_fixed_point_(1e-4), successes_(0)
{
  double dx = dft_->get_lattice().getDX();
  dt_ = 10*0.1*dx*dx;
  dt_ = 0.0001*dx*dx;
  dtMax_ = 1; //1*dx*dx;

  F_ = dft_->calculateFreeEnergyAndDerivatives(false); //density_, 0.0, dF_,true);   

  int Jspecies = 0;
  Species *species = dft_->getSpecies(Jspecies);  
  const Lattice &lattice = species->getLattice();
  
  Nx_ = lattice.Nx();
  Ny_ = lattice.Ny();
  Nz_ = lattice.Nz();

  dx_ = lattice.getDX();
  dy_ = lattice.getDY();
  dz_ = lattice.getDZ();

  RHS0_.initialize(Nx_,Ny_,Nz_);
  RHS1_.initialize(Nx_,Ny_,Nz_);
  
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

double DDFT::step()
{
  if(dft_->getNumberOfSpecies() > 1)
    throw std::runtime_error("DDFT only implemented for single species systems ... aborting");
  
  F_ = 0;
  try {
    F_ = dft_->calculateFreeEnergyAndDerivatives(false);   
    calls_++;
  } catch( Eta_Too_Large_Exception &e) {
    throw e;
  }
  //  cout << "Initial F = " << F_ << endl;

  
  int Jspecies = 0;
  Species *species = dft_->getSpecies(Jspecies);
  
  const Lattice &lattice = species->getLattice();
  const Density &density = species->getDensity();

  // copy the current density
  DFT_Vec d0(density.Ntot()); d0.set(density.get_density_real());
  DFT_Vec d1(density.Ntot()); d1.set(d0);

  double deviation = 1;
  double dt = dt_;
  bool   reStart = false;
  bool   decreased_time_step = false;

  calcNonlinearTerm(d0, species, RHS0_);
  
  do {
    reStart = false;
    double old_error = 0;
    
    for(int i=0;i<100 && deviation > tolerence_fixed_point_ && !reStart;i++)
      {
	species->set_density(d1);
	F_ = dft_->calculateFreeEnergyAndDerivatives(false);
	calcNonlinearTerm(d1, species, RHS1_);  

	double nn = species->getDensity().getNumberAtoms();
	
	species->set_density(d0);       
	species->fft_density();
	
	deviation = fftDiffusion(d1);

	// decrease time step and restart if density goes negative or if error is larger than previous step
	if(d1.min() < 0 || (i > 0 && old_error < deviation)) {reStart = true; dt_ /= 10; d1.set(d0); decreased_time_step = true;}

	old_error = deviation;	       	
      }
    if(!reStart && deviation > tolerence_fixed_point_)
      {reStart = true; dt_ /= 10; d1.set(d0); decreased_time_step = true;}
  } while(reStart);

  time_ += dt_;
  
  // Adaptive time-step: try to increase time step if the present one works 5 times 
  if(decreased_time_step) successes_ = 0;
  else successes_++;
  if(successes_ >= 5 && dt_ < dtMax_) { dt_ = min(2*dt, dtMax_); successes_ = 0;}

  species->set_density(d1);
  species->fft_density();

  F_ = dft_->calculateFreeEnergyAndDerivatives(false); 

  //  cout << setprecision(12) << "F = " << F_ << " Natoms = " << species->getDensity().getNumberAtoms() << endl;
  
  return F_;  
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

void DDFT::calcNonlinearTerm(const DFT_Vec &density, Species *species, DFT_FFT& RHS) const
{
  g_dot_x(species->getDF(), RHS.Real());    
  subtract_ideal_gas(density,RHS.Real());

  RHS.do_real_2_fourier();  
}

// This function takes the input density and calculates a new density, d1, by propagating the density a time step dt.
// The nonlinear terms, RHS0 and RHS1, are treated explicitly.
double DDFT::fftDiffusion(DFT_Vec &new_density)
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
	      x += ((exp_dt-1)/Lambda)*RHS0_.cFour().get(pos);
	      x += ((exp_dt-1-dt_*Lambda)/(Lambda*Lambda*dt_))*(RHS1_.cFour().get(pos)-RHS0_.cFour().get(pos));
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

// Note that this must always be consistent with the definition of Lam in the constructor.
void DDFT::subtract_ideal_gas(const DFT_Vec &density, DFT_Vec& RHS) const
{
  const double D[]       = {1.0/(dx_*dx_), 1.0/(dy_*dy_), 1.0/(dz_*dz_)};
  
  unsigned pos;
  //#pragma omp parallel for  private(pos) schedule(static)
  for(pos = 0;pos<RHS.size();pos++)
    {      
      double dpx,dmx,dpy,dmy,dpz,dmz; // density
      double d0 = get_neighbors(density, 0, pos, 1, dpx,dmx,dpy,dmy,dpz,dmz);
      RHS.IncrementBy(pos, -D[0]*(dpx+dmx-2*d0)-D[1]*(dpy+dmy-2*d0)-D[2]*(dpz+dmz-2*d0));
    }
}

double DDFT::get_neighbors(const DFT_Vec &x, int species, long pos, int stride,
			      double &xpx, double &xmx, double &xpy, double &xmy, double &xpz, double &xmz) const
{
  const Density &density = dft_->getDensity(species);
  
  int ix, iy,iz;
  density.cartesian(pos,ix,iy,iz);

  xpx = x.get(density.get_PBC_Pos(ix+stride,iy,iz));
  xmx = x.get(density.get_PBC_Pos(ix-stride,iy,iz));

  xpy = x.get(density.get_PBC_Pos(ix,iy+stride,iz));
  xmy = x.get(density.get_PBC_Pos(ix,iy-stride,iz));

  xpz = x.get(density.get_PBC_Pos(ix,iy,iz+stride));
  xmz = x.get(density.get_PBC_Pos(ix,iy,iz-stride));

  return x.get(density.get_PBC_Pos(ix,iy,iz));
}


// This calculates (del_IJ rho_J del_JK) x_K i.e. matrix g is discretized (del rho del) operator.
// It could be more efficient: as it stands, for each new entry at pos, we retrieve 13 values from the two matrices:
// if instead there was a triple loop over the three directions, then we would only need two new values for each position.
// For now, I prefer the more compact version that is cleaner to read. 
//
// N.B.: For FIXED boundaries, the input vector x must have zero entries on the boundary.
void DDFT::g_dot_x(const DFT_Vec& x, DFT_Vec& gx) const
{
  if(dft_->getNumberOfSpecies() > 1) throw std::runtime_error("DDFT::g_dot_x is not implemented for more than one species");
  int species = 0;

  const int stride       = 1; // 2 for centered differences  
  const Density &density = dft_->getDensity(species);
  const double D[]       = {0.5/(stride*dx_*dx_), 0.5/(stride*dy_*dy_), 0.5/(stride*dz_*dz_)};

  long pos;
  //#pragma omp parallel for  private(pos) schedule(static)
  for(pos = 0;pos<gx.size();pos++)
    {      
      if(is_fixed_boundary() && density.is_boundary_point(pos))
	{
	  gx.set(pos,0.0);
	  if(fabs(x.get(pos) > 0.0))
	    throw std::runtime_error("For fixed boundaries, input vector to DDFT::g_dot_x must have zero boundary entries");
	} else {      
	double xpx,xmx,xpy,xmy,xpz,xmz;
	double x0 = get_neighbors(x,species,pos,stride,xpx,xmx,xpy,xmy,xpz,xmz);
	
	double dpx,dmx,dpy,dmy,dpz,dmz; // density
	double d0 = density.get_neighbors(pos,dpx,dmx,dpy,dmy,dpz,dmz);
	
	// centered differences: replace dpx+d0 ->dpx, etc, i.e. just set d0 = 0 ... 
	gx.set(pos,D[0]*((dpx+d0)*(xpx-x0)-(d0+dmx)*(x0-xmx))
	       + D[1]*((dpy+d0)*(xpy-x0)-(d0+dmy)*(x0-xmy))
	       + D[2]*((dpz+d0)*(xpz-x0)-(d0+dmz)*(x0-xmz)));
      }
    }      
}

void DDFT::matrix_dot_v(const vector<DFT_FFT> &v, vector<DFT_Vec> &result, void *param) const
{
  if(dft_->getSpecies(0)->is_fixed_boundary() != is_fixed_boundary()) throw std::runtime_error("DDFT::matrix_dot_v: DDFT and DFT object out of synch on fixed boundaries");
  
  dft_->matrix_dot_v(v, result);
  DFT_Vec intermediate_result(result[0]);
  result[0].zeros();
  g_dot_x(intermediate_result, result[0]);
}






//void DDFT::reverseForce(DFT_Vec *tangent) 
//{
//  double prod = dF_.dotWith(*tangent);
//  dF_.Increment_And_Scale(*tangent,-2*prod);
//}

void DDFT::Display(double F, double dFmin, double dFmax, double N)
{
  /*
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
  */
}

  /*
double DDFT::F_string(Density &original_density, double *fmax)
{
  double F = 0;

  DFT_Vec dummy;
  double F =  dft_.calculateFreeEnergyAndDerivatives(original_density,0.0, dummy,false);

  if(fmax) *fmax = dummy.inf_norm()/original_density.dV();;
  return F;
}
  */


