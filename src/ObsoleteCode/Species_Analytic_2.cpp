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


FMT_Species_Analytic_2::FMT_Species_Analytic_2(Density& density, double hsd, double mu, int seq):
  FMT_Species(density,hsd,mu,seq)
{
  generateWeights();

  for(FMT_Weighted_Density &d: d_)
    d.transformWeights();
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

void FMT_Species_Analytic_2::generateWeights()
{
  cout << "Generating weights" << endl;
  
  double dx = density_.getDX();
  double dy = density_.getDY();
  double dz = density_.getDZ();

  double dV = dx*dy*dz;
  double dS = dx*dy;

  double hsr = hsd_/2; // the hard-sphere radius
  
  // This saves having to a code a useless special case
  if(hsd_ < dx)
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
  
  for(int Sx = 0; Sx <= Sx_max; Sx++)
    for(int Sy = 0; Sy <= Sy_max; Sy++)
      for(int Sz = 0; Sz <= Sz_max; Sz++)
	{
	  counter++;
	  if(counter%1000 == 0) {if(counter > 0) cout << '\r'; cout << "\t" << int(double(counter)*100.0/pmax) << "% finished: " << counter << " out of " << pmax; cout.flush();}

	  double R2_min = (Sx-(Sx == 0 ? 0 : 1))*(Sx-(Sx == 0 ? 0 : 1))+(Sy-(Sy == 0 ? 0 : 1))*(Sy-(Sy == 0 ? 0 : 1))+(Sz-(Sz == 0 ? 0 : 1))*(Sz-(Sz == 0 ? 0 : 1));
	  double R2_max = (Sx+1)*(Sx+1)+(Sy+1)*(Sy+1)+(Sz+1)*(Sz+1);
	  
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
				w_s   += dS*sgn*(G_s(R,Vx,Vy,Vz,Tx,Ty,Tz)   - G_s(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));

				w_v[0] += px*(dS/R)*sgn*(G_vx(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_vx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				w_v[1] += py*(dS/R)*sgn*(G_vx(R,Vy,Vx,Vz,Ty,Tx,Tz) - G_vx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz));
				w_v[2] += pz*(dS/R)*sgn*(G_vx(R,Vz,Vy,Vx,Tz,Ty,Tx) - G_vx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx));

				w_T[0][0] += (dS/(R*R))*sgn*(G_txx(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_txx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				w_T[1][1] += (dS/(R*R))*sgn*(G_txx(R,Vy,Vx,Vz,Ty,Tx,Tz) - G_txx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz));
				w_T[2][2] += (dS/(R*R))*sgn*(G_txx(R,Vz,Vy,Vx,Tz,Ty,Tx) - G_txx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx));			    

				w_T[0][1] += px*py*(dS/(R*R))*sgn*(G_txy(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				w_T[0][2] += px*pz*(dS/(R*R))*sgn*(G_txy(R,Vx,Vz,Vy,Tx,Tz,Ty) - G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vz,Vy,Tx,Tz,Ty));
				w_T[1][2] += py*pz*(dS/(R*R))*sgn*(G_txy(R,Vy,Vz,Vx,Ty,Tz,Tx) - G_txy(R,sqrt(R*R-Vz*Vz-Vx*Vx),Vz,Vx,Ty,Tz,Tx));			    			    			    
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
		  d_[EI()].addToWeight(pos,w_eta);
		  if(isnan(d_[EI()].getWeight(pos)))
		    {
		      cout << ix << " " << iy << " " << iz << " " << Sx << " " << Sy << " " << Sz << endl;
		      throw std::runtime_error("FOund NAN");
		    }
		  d_[SI()].addToWeight(pos,w_s);
		  for(int iv = 0;iv < 3;iv++)
		    {
		      d_[VI(iv)].addToWeight(pos,(iv == 0 ? (1-2*ix) : (iv == 1 ? (1-2*iy) : (1-2*iz)))*w_v[iv]);
		      for(int it=iv;it<3;it++)
			d_[TI(iv,it)].addToWeight(pos,(iv == 0 ? (1-2*ix) : (iv == 1 ? (1-2*iy) : (1-2*iz)))*(it == 0 ? (1-2*ix) : (it == 1 ? (1-2*iy) : (1-2*iz)))*w_T[iv][it]);
		    }
		}		  
	}

  cout << endl;
  cout << myColor::GREEN;
  cout << "/////  Finished.  " << endl;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << myColor::RESET << endl;
}
