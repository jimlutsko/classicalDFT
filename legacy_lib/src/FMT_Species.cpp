#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>

#include <gsl/gsl_integration.h>

using namespace std;

#ifdef USE_MPI
#include <mpi.h>  
#endif

#include "Species.h"
#include "myColor.h"

#include "FMT_Helper.h"

double FMT_Helper::SMALL_NUM = 1e-15;

int Species::SequenceNumber_ = 0;


// These functions just alias the density so that it is always non-zero and in fact bigger than SMALL_VALUE.
void Species::set_density_from_alias(const DFT_Vec &x)
{
  long pos;
  const double dmin = SMALL_VALUE;
  
#pragma omp parallel for  private(pos)  schedule(static)
  for(pos=0;pos<x.size();pos++)
    density_.set(pos,dmin+x.get(pos)*x.get(pos));
}
  
void Species::get_density_alias(DFT_Vec &x) const
{
  long pos;
  const double dmin = SMALL_VALUE;
  
#pragma omp parallel for  private(pos)  schedule(static)				    
  for(pos=0;pos<x.size();pos++)
    x.set(pos, sqrt(std::max(0.0, density_.get(pos)-1e-20)));    
}

void Species::convert_to_alias_deriv(DFT_Vec &x, DFT_Vec &dF_dRho) const
{
  dF_dRho.Schur(x,dF_dRho);
  dF_dRho.MultBy(2.0);
}


// These impose density limits based on the fact that in the most extreme case,
// if the density at some point is dmax and all other points it is dmin, then
// eta at that point is (4*pi/3)*hsd*hsd*hsd*dmin + (dmax-dmin)*dV.
// So the largest possible value of dmax, based on eta < 1, is
//    dmax = dmin + (eta-etamin)/dV
// Note that the limits used here are really only valid for FMT_Species
// but for now, those are really the only ones used.

void FMT_Species::set_density_from_alias(const DFT_Vec &x)
{
  //Species::set_density_from_alias(x);
  
  long pos;
  const double dmin = SMALL_VALUE;
  const double etamin = dmin*(4*M_PI*getHSD()*getHSD()*getHSD()/3);
  const double c = (1.0-etamin)/density_.dV();
  //const double c = (0.99-etamin)/density_.dV();
  
#pragma omp parallel for  private(pos)  schedule(static)				    
  for(pos=0;pos<x.size();pos++)
    {
      double y = x.get(pos);
      //double z = dmin +c*(1-exp(-y*y));
      double z = dmin +c*y*y/(1+y*y);
      density_.set(pos,z);
    }
  
}

  
void FMT_Species::get_density_alias(DFT_Vec &x) const
{
  //Species::get_density_alias(x);
   
  long pos;
  double dv = density_.dV();
  
  const double dmin = SMALL_VALUE;
  const double etamin = dmin*(4*M_PI*getHSD()*getHSD()*getHSD()/3);
  const double c = (1.0-etamin)/density_.dV();
  //const double c = (0.99-etamin)/density_.dV();
  
#pragma omp parallel for  private(pos)  schedule(static)				    
  for(pos=0;pos<x.size();pos++)
    {
      double z = (density_.get(pos) - dmin)/c;
      //double y = sqrt(fabs(log(1.0/(1-z))));
      double y = sqrt(fabs(z/(1-z)));
      x.set(pos, y);          
    }
  
}

void FMT_Species::convert_to_alias_deriv(DFT_Vec &x, DFT_Vec &dF_dRho) const
{
  //Species::convert_to_alias_deriv(x,dF_dRho);  
  
  long pos;
  const double dmin = SMALL_VALUE;
  const double etamin = dmin*(4*M_PI*getHSD()*getHSD()*getHSD()/3);
  const double c = (1.0-etamin)/density_.dV();
  //  const double c = (0.99-etamin)/density_.dV();  

#pragma omp parallel for  private(pos)  schedule(static)				      
  for(pos=0;pos<x.size();pos++)
    {
      double y = x.get(pos);
      double df = dF_dRho.get(pos);
      //dF_dRho.set(pos, df*(2*c*y*exp(-y*y)));
      dF_dRho.set(pos, df*(2*c*y/((1+y*y)*(1+y*y))));
    }
  
}

FMT_Species::FMT_Species(Density& density, double hsd, double mu, int seq): Species(density,mu,seq), hsd_(hsd), fmt_weighted_densities(11)
{
  long Nx = density_.Nx();
  long Ny = density_.Ny();
  long Nz = density_.Nz();

  for(FMT_Weighted_Density &d: fmt_weighted_densities)
    d.initialize(Nx, Ny, Nz);

  //  generateWeights(hsd, fmt_weighted_densities);
  FMT_Weighted_Density *data = fmt_weighted_densities.data();
  generateWeights(hsd, data + EI(), data + SI(), data+VI(0), data+TI(0,0));
  
  for(FMT_Weighted_Density &d: fmt_weighted_densities)
    d.transformWeights();
}

void FMT_Species::generateWeights(double hsd, vector<FMT_Weighted_Density> &fmt_weights)
{
  generateWeights(hsd, fmt_weights.data()+EI(), fmt_weights.data()+SI(), fmt_weights.data()+VI(0), fmt_weights.data()+TI(0,0));

}
void FMT_Species::generateWeights(double hsd, FMT_Weighted_Density* eta, FMT_Weighted_Density* s, FMT_Weighted_Density* v, FMT_Weighted_Density *T)
{
  cout << "OK: Generating weights" << endl;
  
  double dx = density_.getDX();
  double dy = density_.getDY();
  double dz = density_.getDZ();

  double dV = dx*dy*dz;
  
  double dS = dx*dx;

  double hsr = hsd/2; // the hard-sphere radius
  
  // This saves having to a code a useless special case
  if(hsd < min(min(dx,dy),dz))
    throw std::runtime_error("hsd is less than the lattice spacing ... aborting");
  
  int Sx_max = 2+int(hsr/dx);
  int Sy_max = 2+int(hsr/dy);
  int Sz_max = 2+int(hsr/dz);

  long pmax = (Sx_max+1)*(Sy_max+1)*(Sz_max+1);
  
  int I[2] = {-1,1};

  cout << endl;
  cout << myColor::GREEN;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << "/////  Generating weights using analytic formulae" << endl;
  cout << myColor::RESET << endl;

  long counter = 0;

  //REGARDING dx,dy,dz: In the following, I allow for different spacings in the different directions. In order to preserve - exactly - the previous results, I essentially
  // scale everything by dx. A more rationale implementation would involve changes as noted below in comments marked with !!!. 

  //!!! do not scale by dx  
  double R = hsr/dx; 
  dy /= dx; 
  dz /= dx;
  dx = 1;
  
  for(int Sx = 0; Sx <= Sx_max; Sx++)
    for(int Sy = 0; Sy <= Sy_max; Sy++)
      for(int Sz = 0; Sz <= Sz_max; Sz++)
	{
	  counter++;
	  if(counter%1000 == 0) {if(counter > 0) cout << '\r'; cout << "\t" << int(double(counter)*100.0/pmax) << "% finished: " << counter << " out of " << pmax; cout.flush();}

	  double R2_min = dx*dx*(Sx-(Sx == 0 ? 0 : 1))*(Sx-(Sx == 0 ? 0 : 1))+dy*dy*(Sy-(Sy == 0 ? 0 : 1))*(Sy-(Sy == 0 ? 0 : 1))+dz*dz*(Sz-(Sz == 0 ? 0 : 1))*(Sz-(Sz == 0 ? 0 : 1));
	  double R2_max = dx*dx*(Sx+1)*(Sx+1)+dy*dy*(Sy+1)*(Sy+1)+dz*dz*(Sz+1)*(Sz+1);
	  
	  double w_eta = 0.0;
	  double w_s   = 0.0;
	  double w_v[3] = {0.0,0.0,0.0};
	  double w_T[3][3] = {{0.0,0.0,0.0},
			      {0.0,0.0,0.0},
			      {0.0,0.0,0.0}};

	  // The weight for point Sx,Sy,Sz has contributions from all adjoining cells
	  // The furthest corner is Sx+1,Sy+1,Sz+1 and if this less than hsr*hsr, the volume weights are 1, and the surface weights are zero
	  // else, the nearest corner is Sx-1,Sy-1,Sz-1 (unless Sx, Sy or Sz = 0) and if this is less than hsr*hsr, then the boundary is between these limits and we must compute
	  // else, all hsd boundary is less than the nearest corner and all weights are zero.

	  //Note: Special cases of I-sums, e.g. when one or more components of S are zero, are handled in the called functions.
	  
	  if(R*R > R2_max) {w_eta = dV;}
	  else if(R*R > R2_min)
	    for(int ax:I)
	      {
		int ix = ax;
		double Tx = dx*(Sx+ix);
		int px = (Tx < 0 ? -1 : 1);
		if(Tx < 0) {ix = 1; Tx = dx;}		   
		for(int ay:I)
		  {
		    int iy = ay;
		    double Ty = dy*(Sy+iy);
		    int py = (Ty < 0 ? -1 : 1);
		    if(Ty < 0) {iy = 1; Ty = dy;}	       		    
		    for(int az:I)
		      {
			int iz = az;
			double Tz = dz*(Sz+iz);
			int pz = (Tz < 0 ? -1 : 1);
			if(Tz < 0) {iz = 1; Tz = dz;}
  
			int vx[2] = {0,ix};
			int vy[2] = {0,iy};
			int vz[2] = {0,iz};

			double j = 0.0;

			for(int jx = 0; jx < 2; jx++)
			  for(int jy = 0; jy < 2; jy++)
			    for(int jz = 0; jz < 2; jz++)
			      {
				double Vx = Tx - dx*vx[jx];
				double Vy = Ty - dy*vy[jy];
				double Vz = Tz - dz*vz[jz];
			    
				int sgn = 1;
				if(jx == 1) sgn *= -1;
				if(jy == 1) sgn *= -1;
				if(jz == 1) sgn *= -1;
			    
				// test up to machine precision
				if((R*R - (Vx*Vx+Vy*Vy+Vz*Vz)) < std::nextafter(0.d,1.d)) continue;

				//!!! Replace dV and dS in the following by (1/dV) in all cases
				
				if(eta) w_eta += dV*sgn*(FMT_Helper::G_eta(R,Vx,Vy,Vz,Tx,Ty,Tz) - FMT_Helper::G_eta(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));

				if(s) w_s   += dS*sgn*(FMT_Helper::G_s(R,Vx,Vy,Vz,Tx,Ty,Tz)   - FMT_Helper::G_s(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));

				if(v)
				  {
				    w_v[0] += dS*px*(1.0/R)*sgn*(FMT_Helper::G_vx(R,Vx,Vy,Vz,Tx,Ty,Tz) - FMT_Helper::G_vx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_v[1] += dS*py*(1.0/R)*sgn*(FMT_Helper::G_vx(R,Vy,Vx,Vz,Ty,Tx,Tz) - FMT_Helper::G_vx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz));
				    w_v[2] += dS*pz*(1.0/R)*sgn*(FMT_Helper::G_vx(R,Vz,Vy,Vx,Tz,Ty,Tx) - FMT_Helper::G_vx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx));
				  }
				if(T)
				  {
				    w_T[0][0] += dS*(1.0/(R*R))*sgn*(FMT_Helper::G_txx(R,Vx,Vy,Vz,Tx,Ty,Tz) - FMT_Helper::G_txx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_T[1][1] += dS*(1.0/(R*R))*sgn*(FMT_Helper::G_txx(R,Vy,Vx,Vz,Ty,Tx,Tz) - FMT_Helper::G_txx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz));
				    w_T[2][2] += dS*(1.0/(R*R))*sgn*(FMT_Helper::G_txx(R,Vz,Vy,Vx,Tz,Ty,Tx) - FMT_Helper::G_txx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx));			    
				    
				    w_T[0][1] += dS*px*py*(1.0/(R*R))*sgn*(FMT_Helper::G_txy(R,Vx,Vy,Vz,Tx,Ty,Tz) - FMT_Helper::G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_T[0][2] += dS*px*pz*(1.0/(R*R))*sgn*(FMT_Helper::G_txy(R,Vx,Vz,Vy,Tx,Tz,Ty) - FMT_Helper::G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vz,Vy,Tx,Tz,Ty));
				    w_T[1][2] += dS*py*pz*(1.0/(R*R))*sgn*(FMT_Helper::G_txy(R,Vy,Vz,Vx,Ty,Tz,Tx) - FMT_Helper::G_txy(R,sqrt(R*R-Vz*Vz-Vx*Vx),Vz,Vx,Ty,Tz,Tx));
				  }
			      }
		      }
		  }
	      }	
	  // Add in for all octants of the sphere: take account of parity of vector and tensor quantities
	  for(int ix = 0; ix < (Sx == 0 ? 1 : 2); ix++)
	    for(int iy = 0; iy < (Sy == 0 ? 1 : 2); iy++)
	      for(int iz = 0; iz < (Sz == 0 ? 1 : 2); iz++)
		{		  
		  long pos = density_.get_PBC_Pos((1-2*ix)*Sx,(1-2*iy)*Sy,(1-2*iz)*Sz);	  
		  if(eta) eta->addToWeight(pos,w_eta);
		  if(eta && isnan(eta->getWeight(pos)))
		    {
		      cout << ix << " " << iy << " " << iz << " " << Sx << " " << Sy << " " << Sz << endl;
		      throw std::runtime_error("Found NAN");
		    }
		  if(s) s->addToWeight(pos,w_s);

		  if(v)
		    for(int iv = 0;iv < 3;iv++)
		      v[iv].addToWeight(pos,(iv == 0 ? (1-2*ix) : (iv == 1 ? (1-2*iy) : (1-2*iz)))*w_v[iv]);

		  if(T)
		    for(int iv = 0;iv < 3;iv++)
		      for(int it=iv;it<3;it++)
			T[TI(iv,it) - TI(0,0)].addToWeight(pos,(iv == 0 ? (1-2*ix) : (iv == 1 ? (1-2*iy) : (1-2*iz)))*(it == 0 ? (1-2*ix) : (it == 1 ? (1-2*iy) : (1-2*iz)))*w_T[iv][it]);
		}		  
	}

  cout << endl;
  cout << myColor::GREEN;
  cout << "/////  Finished.  " << endl;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << myColor::RESET << endl;
}


////////////////////////////
// AO model
FMT_AO_Species:: FMT_AO_Species(Density& density, double hsd, double Rp, double reservoir_density, double mu, int seq)
  : Rp_(Rp), reservoir_density_(reservoir_density), FMT_Species(density,hsd,mu,seq), fmt_weighted_densitiesAO_(11)
{
  // get the polymer species weights
  long Nx = density_.Nx();
  long Ny = density_.Ny();
  long Nz = density_.Nz();

  for(FMT_Weighted_Density &d: fmt_weighted_densitiesAO_)
    d.initialize(Nx, Ny, Nz);  
  generateWeights(2*Rp_, fmt_weighted_densitiesAO_);
  for(FMT_Weighted_Density &d: fmt_weighted_densitiesAO_)
    d.transformWeights();

  PSI_.initialize(Nx,Ny,Nz);
}


// The job of this function is to take the information concerning dF/dn_{a}(pos) (stored in DPHI) and to construct the dPhi_dn for partial measures that are stored in fmt_weighted_densities.
// This is done via the call to FMT_Species::set_fundamental_measure_derivatives. 
void FMT_AO_Species::set_fundamental_measure_derivatives(FundamentalMeasures &DPHI, long pos, bool needsTensor)
{
  FMT_Species::set_fundamental_measure_derivatives(DPHI,pos,needsTensor);

  double hsdp = 2*Rp_;
  
  double dPhi_dEta = DPHI.eta;
  double dPhi_dS = (DPHI.s0/(hsdp*hsdp)) + DPHI.s1*((1/hsdp)) + DPHI.s2;
  double dPhi_dV[3] = {DPHI.v2[0] + DPHI.v1[0]/hsdp,
			  DPHI.v2[1] + DPHI.v1[1]/hsdp,
			  DPHI.v2[2] + DPHI.v1[2]/hsdp};
  
  fmt_weighted_densitiesAO_[EI()].Set_dPhi(pos,dPhi_dEta);
  fmt_weighted_densitiesAO_[SI()].Set_dPhi(pos,dPhi_dS);

  //note parity factor for the vector
  for(int j=0;j<3;j++)
    {
      fmt_weighted_densitiesAO_[VI(j)].Set_dPhi(pos, -dPhi_dV[j]);	
      if(needsTensor)
	for(int k=j;k<3;k++)
	  fmt_weighted_densitiesAO_[TI(j,k)].Set_dPhi(pos,(j == k ? 1 : 2)*DPHI.T[j][k]); // taking account that we only use half the entries
    }
}

// This is where we finish the construction of Psi(r) and
// evaluate the contribution to the free energy
double FMT_AO_Species::free_energy_post_process(bool needsTensor)  
{  
  PSI_.Four().zeros();  
  
  // The call to add_to_dPhi s does the fft of dPhi/dn_{a} for each fm and adds to array PSI
  int number_of_weights = (needsTensor ? fmt_weighted_densities.size() : 5);

  int i1;  
  //#pragma omp parallel for   private(i1) schedule(static)	 
  for(i1=0;i1<number_of_weights;i1++)
    fmt_weighted_densitiesAO_[i1].add_weight_schur_dPhi_to_arg(PSI_.Four()); 

  PSI_.do_fourier_2_real();

  // Here, we compute the contribution to the free energy and begin to reuse PSI in preparation for computing the forces
  double F = 0;

  long Ntot = PSI_.cReal().size();
  long i;
  
#pragma omp parallel for \
  shared( PSI_ )  \
  private(i)		 \
  schedule(static)	 \
  reduction(+:F)
  for(long i=0;i<Ntot;i++)
    {
      double val = exp(-PSI_.cReal().get(i));
      PSI_.Real().set(i,val);
      F += val;
    }
  F *= -reservoir_density_;

  // This is the "standard" shift of the free energy. The Ntot eventually gets multiplied by dV to become V. 
  F += reservoir_density_*density_.Ntot();
  
  // This is to prepare for the force calculation which comes later. Strictly speaking, it is not a part of the free energy calculation. 
  // Upsilon requires convoluting the (now exponentiated) PSI with the various weights
  // We do this here because the next step is carried out in FMT::calculateFreeEnergyAndDerivatives and so PSI_ is not accessible.
  PSI_.do_real_2_fourier(); // do FFT

  for(auto &x: fmt_weighted_densitiesAO_)
    x.convoluteWith(PSI_.cFour()); // point-wise multiplication of FFT's and call to fourier_2_real (norm factor was included in definition of weights)
  
  return F;
}
