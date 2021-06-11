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

  generateWeights(hsd, fmt_weighted_densities);

  for(FMT_Weighted_Density &d: fmt_weighted_densities)
    d.transformWeights();

  //  if(hsd > 0.99)
  //    for(int i=2;i<fmt_weighted_densities.size();i++)
  //      fmt_weighted_densities[i].zero();

  
}



static void getI(double X, double A, double I[])
{
  double a = asin(X/A);
  double b = sqrt(fabs(A*A-X*X));
  
  I[0] = 0.5*X*b+0.5*A*A*a;
  I[1] = -(1.0/3)*(A*A-X*X)*b;
  I[2] = 0.125*X*(2*X*X-A*A)*b+0.125*A*A*A*A*a;
  I[3] = -(1.0/15)*(A*A-X*X)*(3*X*X+2*A*A)*b;
  I[4] = (1.0/48)*(8*X*X*X*X-2*A*A*X*X-3*A*A*A*A)*X*b+0.0625*A*A*A*A*A*A*a;
}

static const double SMALL_NUM = 1e-15;
  
static void getJ(double X, double V, double R, double J[])
{
  double aV  = asin(max(-1.0,min(1.0,V/max(SMALL_NUM,sqrt(R*R-X*X)))));
  double aX  = asin(max(-1.0,min(1.0,X/sqrt(R*R-V*V))));
  double aVX = asin(max(-1.0,min(1.0,X*V/max(SMALL_NUM,sqrt((R*R-V*V)*(R*R-X*X))))));
  double b   = sqrt(fabs(R*R-V*V-X*X));
  
  J[0] = X*aV+V*aX-R*aVX;
  J[1] = -0.5*V*b+0.5*(X*X-R*R)*aV;
  J[2] = -(1.0/6)*V*(V*V-3*R*R)*aX+(1.0/3)*X*X*X*aV-(1.0/3)*R*R*R*aVX-(1.0/6)*V*X*b;
  J[3] = 0.25*(X*X*X*X-R*R*R*R)*aV+(M_PI/8)*R*R*R*R-(V/12.0)*(5*R*R+X*X-2*V*V)*b;
  J[4] = 0.025*X*V*(3*V*V-2*X*X-7*R*R)*b+0.025*V*(3*V*V*V*V-10*V*V*R*R+15*R*R*R*R)*aX+0.2*X*X*X*X*X*aV-0.2*R*R*R*R*R*aVX;

  if(isnan(J[1]) || isnan(J[2]) || isnan(J[3]) || isnan(J[4]))
    {
      cout << X << " " << V << " " << R << " | " << b << " " << R*R-V*V << " " << R*R-X*X << " " << R*R - V*V -X*X << endl;
      throw std::runtime_error("done");
    }
  
}


double G_eta(double R, double X, double Vy, double Vz, double Tx, double Ty, double Tz)
{
  double Iy[5];
  double Iz[5];
  double Jy[5];
  double Jz[5];

  // A temporary workaround to avoid NAN.
  R += 1e-15;
  
  getI(X, sqrt(R*R-Vy*Vy), Iy);
  getI(X, sqrt(R*R-Vz*Vz), Iz);
  getJ(X,Vy,R,Jy);
  getJ(X,Vz,R,Jz);

  double A = Ty*Tz*Vy*Vz+0.25*Vy*Vy*Vz*Vz+0.125*R*R*R*R
    -(1.0/6)*(Ty*Vy*Vy*Vy+Tz*Vz*Vz*Vz) -0.5*Tz*Vz*Vy*Vy-0.5*Ty*Vy*Vz*Vz
    +0.125*Vy*Vy*Vy*Vy+0.125*Vz*Vz*Vz*Vz+0.5*R*R*(Ty*Vy+Tz*Vz-0.5*Vy*Vy-0.5*Vz*Vz+0.5*M_PI*Ty*Tz);
  double B = 0.5*Ty*Vy+0.5*Tz*Vz-0.25*Vy*Vy-0.25*Vz*Vz+(M_PI/4)*Ty*Tz+0.25*R*R;
  double C = 0.5*Ty*Vy+(1.0/3)*(R*R-Vy*Vy);
  double D = 0.5*Tz*Vz+(1.0/3)*(R*R-Vz*Vz);

  double g = Tx*A*X-0.5*A*X*X-(1.0/3)*B*Tx*X*X*X+0.25*X*X*X*X*B
    +0.025*Tx*X*X*X*X*X-(1.0/48)*X*X*X*X*X*X;

  if(isnan(g))
    throw std::runtime_error("Here 1");
  
  g -= 0.5*Ty*Tz*(Tx*R*R*Jy[0]-R*R*Jy[1]-Tx*Jy[2]+Jy[3]);
  if(isnan(g))
    throw std::runtime_error("Here 2");  
  g -= 0.5*Ty*Tz*(Tx*R*R*Jz[0]-R*R*Jz[1]-Tx*Jz[2]+Jz[3]);
  if(isnan(g))
    {
      cout << Jz[0] << " " << Jz[1] << " " << Jz[2] << " " << Jz[3] << endl;
      throw std::runtime_error("Here 3");
    }
  g -= (Tx*Tz*C*Iy[0]-Tz*C*Iy[1]-(1.0/3)*Tx*Tz*Iy[2]+(1.0/3)*Tz*Iy[3]);
  if(isnan(g))
    throw std::runtime_error("Here 4");  
  g -= (Tx*Ty*D*Iz[0]-Ty*D*Iz[1]-(1.0/3)*Tx*Ty*Iz[2]+(1.0/3)*Ty*Iz[3]);


  
  return g;
}

double G_s(double R, double X, double Vy, double Vz, double Tx, double Ty, double Tz)
{
  double Iy[5];
  double Iz[5];
  double Jy[5];
  double Jz[5];

  // A temporary workaround to avoid NAN.
  R += 1e-15;
  
  getI(X, sqrt(R*R-Vy*Vy), Iy);
  getI(X, sqrt(R*R-Vz*Vz), Iz);
  getJ(X,Vy,R,Jy);
  getJ(X,Vz,R,Jz);

  double A = Ty*Vy+Tz*Vz-0.5*Vy*Vy-0.5*Vz*Vz+0.5*R*R+(M_PI/2)*Ty*Tz;
  
  double g = 0.125*R*X*X*X*X-(1.0/6)*R*Tx*X*X*X-0.5*R*A*X*X+R*A*Tx*X;

  g -= R*Tz*(Tx*Iy[0]-Iy[1]);
  g -= R*Ty*(Tx*Iz[0]-Iz[1]);

  g -= R*Ty*Tz*(Tx*Jy[0]-Jy[1]);
  g -= R*Ty*Tz*(Tx*Jz[0]-Jz[1]);

  return g;
}

double G_vx(double R, double X, double Vy, double Vz, double Tx, double Ty, double Tz)
{
  double Iy[5];
  double Iz[5];
  double Jy[5];
  double Jz[5];

  // A temporary workaround to avoid NAN.
  R += 1e-15;  

  getI(X, sqrt(R*R-Vy*Vy), Iy);
  getI(X, sqrt(R*R-Vz*Vz), Iz);
  getJ(X,Vy,R,Jy);
  getJ(X,Vz,R,Jz);

  double A = Ty*Vy+Tz*Vz-0.5*Vy*Vy-0.5*Vz*Vz+0.5*R*R+(M_PI/2)*Ty*Tz;
  
  double g = 0.1*R*X*X*X*X*X-0.125*R*Tx*X*X*X*X-(1.0/3)*R*A*X*X*X+0.5*R*A*Tx*X*X;

  g -= R*Tz*(Tx*Iy[1]-Iy[2]);
  g -= R*Ty*(Tx*Iz[1]-Iz[2]);

  g -= R*Ty*Tz*(Tx*Jy[1]-Jy[2]);
  g -= R*Ty*Tz*(Tx*Jz[1]-Jz[2]);

  
  return g;
}

double G_txx(double R, double X, double Vy, double Vz, double Tx, double Ty, double Tz)
{
  double Iy[5];
  double Iz[5];
  double Jy[5];
  double Jz[5];

  // A temporary workaround to avoid NAN.
  R += 1e-15;
  
  getI(X, sqrt(R*R-Vy*Vy), Iy);
  getI(X, sqrt(R*R-Vz*Vz), Iz);
  getJ(X,Vy,R,Jy);
  getJ(X,Vz,R,Jz);

  double A = Ty*Vy+Tz*Vz-0.5*Vy*Vy-0.5*Vz*Vz+0.5*R*R+(M_PI/2)*Ty*Tz;
  
  double g = (1.0/12)*R*X*X*X*X*X*X-0.1*R*Tx*X*X*X*X*X-0.25*R*A*X*X*X*X+(1.0/3)*R*A*Tx*X*X*X;

  g -= R*Tz*(Tx*Iy[2]-Iy[3]);
  g -= R*Ty*(Tx*Iz[2]-Iz[3]);

  g -= R*Ty*Tz*(Tx*Jy[2]-Jy[3]);
  g -= R*Ty*Tz*(Tx*Jz[2]-Jz[3]);

  return g;
}

double G_txy(double R, double X, double Vy, double Vz, double Tx, double Ty, double Tz)
{
  double Iy[5];
  double Iz[5];
  double Jy[5];
  double Jz[5];

  // A temporary workaround to avoid NAN.
  R += 1e-15;
  
  getI(X, sqrt(R*R-Vy*Vy), Iy);
  getI(X, sqrt(R*R-Vz*Vz), Iz);
  getJ(X,Vy,R,Jy);
  getJ(X,Vz,R,Jz);

  double A = (1.0/6)*(3*Ty*R*R+2*Vy*Vy*Vy-3*Ty*Vy*Vy-3*Ty*Vz*Vz+6*Ty*Tz*Vz+1.5*M_PI*R*R*Tz);
  double B = 0.25*M_PI*Tz+0.5*Ty;
  double C = -(1.0/6)*(2*R*R-2*Vz*Vz+3*Tz*Vz);
  
  double g = -0.5*R*A*Tx*X*X+(1.0/3)*R*A*X*X*X+0.25*R*B*Tx*X*X*X*X-0.2*R*B*X*X*X*X*X;

  g += 0.5*R*Tz*(Tx*R*R*Jy[1]-R*R*Jy[2]-Tx*Jy[3]+Jy[4]);
  g += 0.5*R*Tz*(Tx*R*R*Jz[1]-R*R*Jz[2]-Tx*Jz[3]+Jz[4]);

  g += 0.5*R*Tz*(2*Ty-Vy)*(Tx*Iy[1]-Iy[2]);
  g += R*(-Tx*C*Iz[1]+C*Iz[2]-(1.0/3)*Tx*Iz[3]+(1.0/3)*Iz[4]);

  return g;
}

void FMT_Species::generateWeights(double hsd, vector<FMT_Weighted_Density> &fmt_weights)
{
  cout << "OK: Generating weights" << endl;
  
  double dx = density_.getDX();
  double dy = density_.getDY();
  double dz = density_.getDZ();

  double dV = dx*dy*dz;

  double R = hsd/2; // the hard-sphere radius
  
  // This saves having to a code a useless special case
  if(hsd < min(min(dx,dy),dz))
    throw std::runtime_error("hsd is less than the lattice spacing ... aborting");
  
  int Sx_max = 2+int(R/dx);
  int Sy_max = 2+int(R/dy);
  int Sz_max = 2+int(R/dz);

  long pmax = (Sx_max+1)*(Sy_max+1)*(Sz_max+1);
  
  int I[2] = {-1,1};

  cout << endl;
  cout << myColor::GREEN;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << "/////  Generating weights using analytic formulae" << endl;
  cout << myColor::RESET << endl;

  long counter = 0;

  int numWeights = fmt_weights.size();
  
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
				
				w_eta += (1/dV)*sgn*(G_eta(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_eta(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				if(numWeights > 1)
				  {
				    w_s   += (1/dV)*sgn*(G_s(R,Vx,Vy,Vz,Tx,Ty,Tz)   - G_s(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));

				    w_v[0] += (1/dV)*px*(1.0/R)*sgn*(G_vx(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_vx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_v[1] += (1/dV)*py*(1.0/R)*sgn*(G_vx(R,Vy,Vx,Vz,Ty,Tx,Tz) - G_vx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz));
				    w_v[2] += (1/dV)*pz*(1.0/R)*sgn*(G_vx(R,Vz,Vy,Vx,Tz,Ty,Tx) - G_vx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx));

				    w_T[0][0] += (1/dV)*(1.0/(R*R))*sgn*(G_txx(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_txx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_T[1][1] += (1/dV)*(1.0/(R*R))*sgn*(G_txx(R,Vy,Vx,Vz,Ty,Tx,Tz) - G_txx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz));
				    w_T[2][2] += (1/dV)*(1.0/(R*R))*sgn*(G_txx(R,Vz,Vy,Vx,Tz,Ty,Tx) - G_txx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx));			    

				    w_T[0][1] += (1/dV)*px*py*(1.0/(R*R))*sgn*(G_txy(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_T[0][2] += (1/dV)*px*pz*(1.0/(R*R))*sgn*(G_txy(R,Vx,Vz,Vy,Tx,Tz,Ty) - G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vz,Vy,Tx,Tz,Ty));
				    w_T[1][2] += (1/dV)*py*pz*(1.0/(R*R))*sgn*(G_txy(R,Vy,Vz,Vx,Ty,Tz,Tx) - G_txy(R,sqrt(R*R-Vz*Vz-Vx*Vx),Vz,Vx,Ty,Tz,Tx));
				  }
			      }
		      }
		  }
	      }	
	  
	  /*
	  // Check weights and G functions (must set Sx,Sy,Sz manually)
	  // This has been used to debug the dx!=dy!=dz scenario
	  if ((Sx==10 && Sy==0  && Sz==0 ) ||
	      (Sx==0  && Sy==0  && Sz==10))
	  {
		  for (int iii=0; iii<4; iii++)  cout << endl;
		  cout << "Sx = " << Sx << " Sy = " << Sy << " Sz = " << Sz << endl;
		  cout << endl;
		  double Tx = dx*(Sx+1);
		  double Ty = dy*(Sy+1);
		  double Tz = dz*(Sz+1);
		  if (Sx==10) Tx = dx*(Sx-1);
		  if (Sz==10) Tz = dz*(Sz-1);
		  double Vx = Tx;
		  double Vy = Ty;
		  double Vz = Tz;
		  cout << "GG_eta = " << G_eta(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_eta(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz) << endl;
		  cout << "GG_s   = " <<   G_s(R,Vx,Vy,Vz,Tx,Ty,Tz) -   G_s(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz) << endl;
		  cout << "GG_vx  = " <<  G_vx(R,Vx,Vy,Vz,Tx,Ty,Tz) -  G_vx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz) << endl;
		  cout << "GG_vz  = " <<  G_vx(R,Vz,Vy,Vx,Tz,Ty,Tx) -  G_vx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx) << endl;
		  cout << "GG_txx = " << G_txx(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_txx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz) << endl;
		  cout << "GG_tzz = " << G_txx(R,Vz,Vy,Vx,Tz,Ty,Tx) - G_txx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx) << endl;
		  cout << "GG_txy = " << G_txy(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz) << endl;
		  cout << "GG_tzy = " << G_txy(R,Vz,Vy,Vx,Tz,Ty,Tx) - G_txy(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx) << endl;
		  cout << endl;
		  cout << "w_eta  = " << w_eta << endl;
		  cout << "w_s    = " << w_s << endl;
		  cout << "w_v[0] = " << w_v[0] << endl;
		  cout << "w_v[1] = " << w_v[1] << endl;
		  cout << "w_v[2] = " << w_v[2] << endl;
		  cout << "w_T[0][0] = " << w_T[0][0] << endl;
		  cout << "w_T[1][1] = " << w_T[1][1] << endl;
		  cout << "w_T[2][2] = " << w_T[2][2] << endl;
		  cout << "w_T[0][1] = " << w_T[0][1] << endl;
		  cout << "w_T[0][2] = " << w_T[0][2] << endl;
		  cout << "w_T[1][2] = " << w_T[1][2] << endl;
		  for (int iii=0; iii<4; iii++)  cout << endl;
	  }
	  */
	  
	  // Add in for all octants of the sphere: take account of parity of vector and tensor quantities
	  for(int ix = 0; ix < (Sx == 0 ? 1 : 2); ix++)
	    for(int iy = 0; iy < (Sy == 0 ? 1 : 2); iy++)
	      for(int iz = 0; iz < (Sz == 0 ? 1 : 2); iz++)
		{		  
		  long pos = density_.get_PBC_Pos((1-2*ix)*Sx,(1-2*iy)*Sy,(1-2*iz)*Sz);	  
		  fmt_weights[EI()].addToWeight(pos,w_eta);
		  if(isnan(fmt_weights[EI()].getWeight(pos)))
		    {
		      cout << ix << " " << iy << " " << iz << " " << Sx << " " << Sy << " " << Sz << endl;
		      throw std::runtime_error("Found NAN");
		    }
		  if(numWeights > 1)
		    {
		      fmt_weights[SI()].addToWeight(pos,w_s);
		      for(int iv = 0;iv < 3;iv++)
			{
			  fmt_weights[VI(iv)].addToWeight(pos,(iv == 0 ? (1-2*ix) : (iv == 1 ? (1-2*iy) : (1-2*iz)))*w_v[iv]);
			  for(int it=iv;it<3;it++)
			    fmt_weights[TI(iv,it)].addToWeight(pos,(iv == 0 ? (1-2*ix) : (iv == 1 ? (1-2*iy) : (1-2*iz)))*(it == 0 ? (1-2*ix) : (it == 1 ? (1-2*iy) : (1-2*iz)))*w_T[iv][it]);
			}
		    }
		}		  
	}

  cout << endl;
  cout << myColor::GREEN;
  cout << "/////  Finished.  " << endl;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << myColor::RESET << endl;
}


FMT_Species_EOS::FMT_Species_EOS(double D_EOS, Density& density, double hsd, double mu, int seq)
  : FMT_Species(density,hsd,mu,seq), eos_weighted_density_(1), D_EOS_(D_EOS)
										    
{
  long Nx = density_.Nx();
  long Ny = density_.Ny();
  long Nz = density_.Nz();

  eos_weighted_density_[0].initialize(Nx, Ny, Nz);
  generateWeights(D_EOS*hsd, eos_weighted_density_);
  eos_weighted_density_[0].transformWeights();  
}
/*
// This is an exact copy of FMT_Species::generateWeights. It is necessary because we want to do
// the same but with a different hsd. Obviously, these could (and probably SHOULD) be combined
// with some additional logic to only do the work required. For now, this is a quick solution
// with the only drawback being that essentially the same code must be maintained in two different
// places (e.g. if any errors are found). Since the code is fairly compact, I do not see this as a major problem.

void FMT_Species_EOS::generate_additional_Weight()
{
  cout << "Generating additional eos weight" << endl;
  
  double dx = density_.getDX();
  double dy = density_.getDY();
  double dz = density_.getDZ();

  double dV = dx*dy*dz;
  double dS = dx*dy;

  double hsr = hsd_; // twide the hard-sphere radius
  
  // This saves having to a code a useless special case
  if(2*hsd_ < dx)
    throw std::runtime_error("2*hsd is less than the lattice spacing ... aborting");
  
  int Sx_max = 2+int(hsr/dx);
  int Sy_max = 2+int(hsr/dy);
  int Sz_max = 2+int(hsr/dz);

  long pmax = (Sx_max+1)*(Sy_max+1)*(Sz_max+1);
  
  int I[2] = {-1,1};

  cout << endl;
  cout << myColor::GREEN;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << "/////  Generating eos weight using analytic formulae" << endl;
  cout << myColor::RESET << endl;

  long counter = 0;
  
  for(int Sx = 0; Sx <= Sx_max; Sx++)
    for(int Sy = 0; Sy <= Sy_max; Sy++)
      for(int Sz = 0; Sz <= Sz_max; Sz++)
	{
	  counter++;
	  if(counter%1000 == 0) {if(counter > 0) cout << '\r'; cout << "\t" << int(double(counter)*100.0/pmax) << "% finished: " << counter << " out of " << pmax; cout.flush();}

	  double R2_min = (Sx-(Sx == 0 ? 0 : 1))*(Sx-(Sx == 0 ? 0 : 1))+(Sy-(Sy == 0 ? 0 : 1))*(Sy-(Sy == 0 ? 0 : 1))+(Sz-(Sz == 0 ? 0 : 1))*(Sz-(Sz == 0 ? 0 : 1));
	  double R2_max = (Sx+1)*(Sx+1)+(Sy+1)*(Sy+1)+(Sz+1)*(Sz+1);
	  
	  double w_eta = 0.0;

	  // The weight for point Sx,Sy,Sz has contributions from all adjoining cells
	  // The furthest corner is Sx+1,Sy+1,Sz+1 and if this less than hsr*hsr, the volume weights are 1, and the surface weights are zero
	  // else, the nearest corner is Sx-1,Sy-1,Sz-1 (unless Sx, Sy or Sz = 0) and if this is less than hsr*hsr, then the boundary is between these limits and we must compute
	  // else, all hsd boundary is less than the nearest corner and all weights are zero.

	  double R = hsr/dx;

	  //Note: Special cases of I-sums, e.g. when one or more components of S are zero, are handled in the called functions.
	  
	  if(R*R > R2_max) {w_eta = dV;}
	  else if(R*R > R2_min)
	    for(int ax:I)
	      {
		int ix = ax;
		double Tx = Sx+ix;
		int px = (Tx < 0 ? -1 : 1);
		if(Tx < 0) Tx = ix = 1;		   
		for(int ay:I)
		  {
		    int iy = ay;
		    double Ty = Sy+iy;
		    int py = (Ty < 0 ? -1 : 1);
		    if(Ty < 0) {Ty = 1; iy = 1;}	       		    
		    for(int az:I)
		      {
			int iz = az;
			double Tz = Sz+iz;
			int pz = (Tz < 0 ? -1 : 1);
			if(Tz < 0) Tz = iz = 1;
  
			int vx[2] = {0,ix};
			int vy[2] = {0,iy};
			int vz[2] = {0,iz};

			double j = 0.0;

			for(int jx = 0; jx < 2; jx++)
			  for(int jy = 0; jy < 2; jy++)
			    for(int jz = 0; jz < 2; jz++)
			      {
				double Vx = Tx - vx[jx];
				double Vy = Ty - vy[jy];
				double Vz = Tz - vz[jz];
			    
				int sgn = 1;
				if(jx == 1) sgn *= -1;
				if(jy == 1) sgn *= -1;
				if(jz == 1) sgn *= -1;
			    
				// test up to machine precision
				if((R*R - (Vx*Vx+Vy*Vy+Vz*Vz)) < std::nextafter(0.d,1.d)) continue;

				w_eta += dV*sgn*(G_eta(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_eta(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
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
		  eos_weighted_density_.addToWeight(pos,w_eta);
		  if(isnan(eos_weighted_density_.getWeight(pos)))
		    {
		      cout << ix << " " << iy << " " << iz << " " << Sx << " " << Sy << " " << Sz << endl;
		      throw std::runtime_error("Found NAN");
		    }
		}		  
	}

  cout << endl;
  cout << myColor::GREEN;
  cout << "/////  Finished.  " << endl;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << myColor::RESET << endl;
}
*/


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
