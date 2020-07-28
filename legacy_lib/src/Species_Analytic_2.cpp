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
  double b = sqrt(A*A-X*X);
  
  I[0] = 0.5*X*b+0.5*A*A*a;
  I[1] = -(1.0/3)*(A*X-X*X)*b;
  I[2] = 0.125*X*(2*X*X-A*A)*b+0.125*a;
  I[3] = -(1.0/15)*(A*A-X*X)*(3X*X+2*A*A)*b;  
}

static void getJ(double X, double V, double J[])
{
  double aV  = asin(V/sqrt(R*R-X*X));
  double aX  = asin(X/sqrt(R*R-V*V));
  double aVX = asin(X*V/sqrt((R*R-V*V)*(R*R-X*X)));
  double b   = sqrt(R*R-V*V-X*X);
  
  J[0] = X*aV+V*aX-R*aVX;
  J[1] = -0.5*C*b+0.5*(X*X-R*R)*aV;
  J[2] = -(1.0/6)*V*(V*V-3*R*R)*aX+(1.0/3)*X*X*X*aV-(1.0/3)*R*R*R*aVX-(1.0/6)*V*X*b;
  J[3] = 0.25*(X*X*X*X-R*R*R*R)*aV+(M_PI/8)*R*R*R*R-(V/12.0)*(5*R*R+X*X-2*V*V)*b;
}


double G_eta(double R, double X, double Vy, double Vz, double Tx, double Ty, double Tz)
{
  double Iy[4];
  double Iz[4];
  double Jy[4];
  double Jz[4];

  getI(X, sqrt(R*R-Vy*Vy), Iy);
  getI(X, sqrt(R*R-Vz*Vz), Iz);
  getJ(X,Vy,Jy);
  getJ(X,Vz,Jz);

  double A = 2*Ty*Tz*Vy*Vz+0.5*Vy*Vy*Vz*Vz+0.25*R*R*R*R
    -(1.0/3)*(Ty*Vy*Vy*Vy+Tz*Vz*Vz*Vz) -Tz*Vz*Vy*Vy-Ty*Vy*Vz*Vz;
  double B = Ty*Vy+Tz*Vz-0.5*Vy*Vy-0.5*Vz*Vz+(M_PI/2)*Ty*Tz+0.5*R*R;
  double C = Ty*Vy+(2.0/3)*(R*R-Vy*Vy);
  double D = Tz*Vz+(2.0/3)*(R*R-Vz*Vz);

  
  doule g = Tx*A*X-0.5*A*X*X-(1.0/3)*B*TxX*X*X+0.25*X*X*X*X
    +0.05*Tx*X*X*X*X*X-(1.0/24)*X*X*X*X*X*X;

  g -= Ty*Tz*(Tx*R*R*Jy[0]-R*R*Jy[1]-Tx*Jy[2]+Jy[3]);
  g -= Ty*Tz*(Tx*R*R*Jz[0]-R*R*Jz[1]-Tx*Jz[2]+Jz[3]);

  g -= (Tx*Tz*C*Iy[0]-Tz*C*Iy[1]-(2.0/3)*Yx*Tz*Iy[2]+(2.0/3)*Tz*Iy[3]);
  g -= (Tx*Ty*C*Iz[0]-Ty*C*Iz[1]-(2.0/3)*Yx*Ty*Iz[2]+(2.0/3)*Ty*Iz[3]);

  return g;
}

double G_s(double R, double X, double Vy, double Vz, double Tx, double Ty, double Tz)
{
  double Iy[4];
  double Iz[4];
  double Jy[4];
  double Jz[4];

  getI(X, sqrt(R*R-Vy*Vy), Iy);
  getI(X, sqrt(R*R-Vz*Vz), Iz);
  getJ(X,Vy,Jy);
  getJ(X,Vz,Jz);

  double A = Ty*Vy+Tz*Vz-0.5*Vy*Vy-0.5*Vz*Vz+0.5*R*R+(M_PI/2)*Ty*Tz;
  
  doule g = 0.125*R*X*X*X*X-(1.0/6)*R*Tx*X*X*X-0.5*R*A*X*X+R*A*Tx*X;

  g -= R*Tz*(Tx*Iy[0]-Iy[1]);
  g -= R*Tt*(Tx*Iz[0]-Iz[1]);

  g -= R*Ty*Tz*(Tx*Jy[0]-Jy[1]);
  g -= R*Ty*Tz*(Tx*Jz[0]-Jz[1]);

  return g;
}

double G_vx(double R, double X, double Vy, double Vz, double Tx, double Ty, double Tz)
{
  double Iy[4];
  double Iz[4];
  double Jy[4];
  double Jz[4];

  getI(X, sqrt(R*R-Vy*Vy), Iy);
  getI(X, sqrt(R*R-Vz*Vz), Iz);
  getJ(X,Vy,Jy);
  getJ(X,Vz,Jz);

  double A = Ty*Vy+Tz*Vz-0.5*Vy*Vy-0.5*Vz*Vz+0.5*R*R+(M_PI/2)*Ty*Tz;
  
  doule g = 0.1*R*X*X*X*X*X-0.125*R*Tx*X*X*X*X-(1.0/3)*R*A*X*X*X+0.5*R*A*Tx*X*X;

  g -= R*Tz*(Tx*Iy[1]-Iy[2]);
  g -= R*Tt*(Tx*Iz[1]-Iz[2]);

  g -= R*Ty*Tz*(Tx*Jy[1]-Jy[2]);
  g -= R*Ty*Tz*(Tx*Jz[1]-Jz[2]);

  return g;
}

double G_txx(double R, double X, double Vy, double Vz, double Tx, double Ty, double Tz)
{
  double Iy[4];
  double Iz[4];
  double Jy[4];
  double Jz[4];

  getI(X, sqrt(R*R-Vy*Vy), Iy);
  getI(X, sqrt(R*R-Vz*Vz), Iz);
  getJ(X,Vy,Jy);
  getJ(X,Vz,Jz);

  double A = Ty*Vy+Tz*Vz-0.5*Vy*Vy-0.5*Vz*Vz+0.5*R*R+(M_PI/2)*Ty*Tz;
  
  doule g = (1.0/12)*R*X*X*X*X*X*X-0.1*R*Tx*X*X*X*X*X-0.25*R*A*X*X*X*X+(1.0/3)R*A*Tx*X*X*X;

  g -= R*Tz*(Tx*Iy[2]-Iy[3]);
  g -= R*Tt*(Tx*Iz[2]-Iz[3]);

  g -= R*Ty*Tz*(Tx*Jy[2]-Jy[3]);
  g -= R*Ty*Tz*(Tx*Jz[2]-Jz[3]);

  return g;
}

double G_txy(double R, double X, double Vy, double Vz, double Tx, double Ty, double Tz)
{
  double Iy[4];
  double Iz[4];
  double Jy[4];
  double Jz[4];

  getI(X, sqrt(R*R-Vy*Vy), Iy);
  getI(X, sqrt(R*R-Vz*Vz), Iz);
  getJ(X,Vy,Jy);
  getJ(X,Vz,Jz);

  double A = (1.0/6)*(3*Ty*R*R+2*Vy*Vy*Vy-3*Ty*Vy*Vy-3*Ty*Vz*Vz+6*Ty*Tz*Vz_1.5*M_PI*R*R*Tz);
  double B = 0.25*M_PI*Tz+0.5*Ty;
  double C = -(1.0/6)*(2*R*R-2*Vz*Vz+3*Tz*Vz);
  double D = Tx*C;
  
  doule g = -0.5*R*A*Tx*X*X+(1.0/3)*R*A*X*X*X+0.25*R*B*Tx*X*X*X*X-0.2*R*B*X*X*X*X*X;

  g += 0.5*R*Tz*(Tx*R*R*Jy[1]-R*R*Jy[2]-Tx*Jy[3]+Jy[4]);
  g += 0.5*R*Tz*(Tx*R*R*Jz[1]-R*R*Jz[2]-Tx*Jz[3]+Jz[4]);

  g += 0.5*R*Tz*(2*Ty-Vy)*(Tx*Iy[1]-Iy[2]);
  g += R*(-Tx*C*Iz[1]+C*Iz[2]-(1.0/3)*Tx*Iz[3]+(1.0/3)*Iz[4]);

  return g;
}


double J_eta(double R, double Tx, double Ty, double Tz, int Ix, int Iy, int Iz)
{
  // Handle case of Sx = 0, Ix = -1 ==> Tx = -1, etc. For eta, parity is always equal to 1
  
  if(Tx == -1)
    Tx = Ix = 1;

  if(Ty == -1)
    Ty = Iy = 1;

  if(Tz == -1)
    Tz = Iz = 1;
  
  int vx[2] = {0,Ix};
  int vy[2] = {0,Iy};
  int vz[2] = {0,Iz};

  double j = 0.0;

  for(int ix = 0; ix < 2; ix++)
    for(int iy = 0; iy < 2; iy++)
      for(int iz = 0; iz < 2; iz++)
	{
	  double Vx = Tx - vx[ix];
	  double Vy = Ty - vy[iy];
	  double Vz = Tz - vz[iz];

	  int sgn = 1;
	  if(ix == 1) sgn *= -1;
	  if(iy == 1) sgn *= -1;
	  if(iz == 1) sgn *= -1;
			  
	  // if(R*R < a*a+b*b+c*c) continue;
	  // test up to machine precision
	  if((R*R - (Vx*Vx+Vy*Vy+Vz*Vz)) < std::nextafter(0.d,1.d)) continue;

	  j += sgn*(G_eta(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_eta(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
	}
  
  return j;
}

double J_Vz_helper(double R, double Sx, double Sy, double Sz, double c, double b, double a)
{
  double j1 = 0.0;
  
  j1 += a*(15*a*a*a*Sz+3*a*a*a*a+40*b*b*b*Sy+80*a*a*b*Sy+60*b*b*c*Sx-120*a*b*Sy*Sz-120*b*c*Sx*Sy)/240.0;
  j1 += b*b*(60*a*a*Sz+20*a*a*a+30*b*b*Sz+30*b*b*a-40*b*Sy*Sz-40*a*b*Sy-60*c*Sx*Sz)/240.0;
  j1 += c*c*(30*b*b*Sz+30*b*b*a-60*b*Sy*Sz-60*a*b*Sy)/240.0;

  j1 += R*R*(3*(4*b*Sy+R*R)*(a+Sz)-6*a*a*Sz-2*a*a*a-24*a*b*Sy)/48.0;
  j1 += -R*R*b*(b-Sy)*(a+Sz)/4;
  
  j1 += -Sx*Sy*(3*(R*R-a*a)*(R*R-a*a)+4*a*Sz*(3*R*R-a*a))*(a*b/((R*R-a*a)*sqrt(R*R-a*a-b*b)))/24.0;
  j1 += -Sx*Sy*(2*a*a*a+3*Sz*(R*R-a*a))*(asin(b/sqrt(R*R-a*a))-M_PI/4)/6.0;
  j1 += -Sx*(3*(R*R-b*b)*(R*R-b*b)+4*b*Sy*(3*R*R-b*b))*(asin(a/sqrt(R*R-b*b))-M_PI/4)/24.0;
  j1 += -Sx*Sz*(3*(R*R-b*b)*(R*R-b*b)+4*b*Sy*(3*R*R-b*b))*(1.0/sqrt(R*R-a*a-b*b))/24.0;
  j1 += -Sy*(3*(R*R-c*c)*(R*R-c*c)+4*c*Sx*(3*R*R-c*c))*(asin(b/sqrt(R*R-c*c))-M_PI/4)/24.0;
  
  j1 += -Sx*(16*a*b*b-30*Sy*a*b+20*b*Sy*Sz+22*a*a*a-30*Sz*a*a-7*R*R*a+25*Sz*R*R)*sqrt(R*R-a*a-b*b)/120.0;
  j1 += -Sx*(-9*a*b*b-25*Sz*b*b+20*b*Sy*Sz+20*a*b*Sy)*sqrt(R*R-a*a-b*b)/120.0;
  j1 += -(4*R*R*R*R+8*c*c*b*b-25*Sy*c*c*b+20*c*b*Sy*Sx+8*c*c*c*c-10*Sx*c*c*c-16*R*R*c*c+25*Sx*R*R*c)*sqrt(R*R-c*c-b*b)/120.0;  
  j1 += Sx*(4*R*R*R*R+8*a*a*b*b-25*Sy*a*a*b+20*a*b*Sy*Sz+8*a*a*a*a-10*Sz*a*a*a-16*R*R*a*a+25*Sz*R*R*a)*(a/sqrt(R*R-a*a-b*b))/120.0;
  j1 += Sx*(4*R*R*R*R+8*a*a*b*b-25*Sz*a*b*b+20*a*b*Sy*Sz+8*b*b*b*b-10*Sy*b*b*b-16*R*R*b*b+25*Sy*R*R*b)*(a/sqrt(R*R-a*a-b*b))/120.0;
  
  return j1;
}

double J_Vz(double R, double Sx, double Sy, double Sz, int wx, int wy, int wz)
{
  int eps = 1;
  if(Sx < 0)
  {
    if(wx != -1) throw std::runtime_error("What is this x1?");
    if(Sx != -1) throw std::runtime_error("What is this x2?");
    else Sx = wx = 1;
  }

  if(Sy < 0)
  {
    if(wy != -1) throw std::runtime_error("What is this y1?");
    if(Sy != -1) throw std::runtime_error("What is this y2?");
    else Sy = wy = 1;
  }

  if(Sz < 0)
  {
    if(wz != -1) throw std::runtime_error("What is this z1?");
    if(Sz != -1) throw std::runtime_error("What is this z2?");
    else {Sz = wz = 1; eps = -1;}
  }  

  
  int va[2] = {0,wz};
  int vb[2] = {0,wy};
  int vc[2] = {0,wx};

  double j = 0.0;

  for(int ia = 0; ia < 2; ia++)
    for(int ib = 0; ib < 2; ib++)
      for(int ic = 0; ic < 2; ic++)
	{
	  double a = Sz - va[ia];
	  double b = Sy - vb[ib];
	  double c = Sx - vc[ic];

	  int sgn = 1;
	  if(ia == 1) sgn *= -1;
	  if(ib == 1) sgn *= -1;
	  if(ic == 1) sgn *= -1;
			  
	  // if(R*R < a*a+b*b+c*c) continue;
	  // test up to machine precision
	  if((R*R - (a*a+b*b+c*c)) < std::nextafter(0.d,1.d)) continue;	  

	  double j0 = 0.0;
	  
	  j0 += -a*(R*R-a*a-b*b-c*c)*(R*R-a*a-b*b-c*c)/8.0;
	  j0 += b*c*Sx*Sy*Sz+a*b*c*Sx*Sy;
	  j0 += -(M_PI/6)*Sx*Sy*R*R*R;

	  j0 += Sx*Sy*R*R*R*(asin(a*b/(sqrt(R*R-a*a)*sqrt(R*R-b*b)))
			     +asin(a*c/(sqrt(R*R-a*a)*sqrt(R*R-c*c)))
			     +asin(b*c/(sqrt(R*R-b*b)*sqrt(R*R-c*c))))/3.0;
	  j0 += Sx*Sy*Sz*R*R*R*R*((b/((R*R-a*a)*sqrt(R*R-a*a-b*b)))
				  +(c/((R*R-a*a)*sqrt(R*R-a*a-c*c))))/3.0;
	  
	  j0 += J_Vz_helper(R,Sx,Sy,Sz,c,b,a);
	  j0 += J_Vz_helper(R,Sy,Sx,Sz,b,c,a);

	  j += sgn * j0;
	}
  
  return j*eps;
}

double J_s_helper(double R, double Sx, double Sy, double Sz, double c, double b, double a)
{
  double j1 = 0.0;
  j1 += R*(6*a*Sz*R*R+12*a*b*Sy*Sz-2*a*a*a*Sz-12*a*a*b*Sy)/24.0;
  j1 += -Sx*Sy*R*(R*R-a*a+2*a*Sz)*(asin(b/sqrt(R*R-a*a))-M_PI/4)/2.0;
  j1 += -Sx*Sy*(3*(R*R-a*a)*(R*R-a*a)+4*a*Sz*(3*R*R-a*a))*(-R*b/((R*R-a*a)*sqrt(R*R-a*a-b*b)))/24.0;
  j1 += Sx*Sy*Sz*R*R*asin(a*b/(sqrt(R*R-a*a)*sqrt(R*R-b*b)))/2.0;
  j1 += Sx*Sy*Sz*R*R*R*(a*b/((R*R-a*a)*(R*R-b*b)))*((a*a+b*b-2*R*R)/sqrt(R*R-a*a-b*b))/6.0;
  j1 += -Sx*(16*R*R*R-32*R*a*a+50*Sz*R*a)*sqrt(R*R-a*a-b*b)/120.0;
  j1 += -Sx*(4*R*R*R*R+8*a*a*b*b-25*Sy*a*a*b+20*a*b*Sy*Sz+8*a*a*a*a-10*Sz*a*a*a-16*R*R*a*a+25*Sz*R*R*a)*(R/sqrt(R*R-a*a-b*b))/120.0;
  return j1;
}

double J_s(double R, double Sx, double Sy, double Sz, int wx, int wy, int wz)
{
  if(Sx < 0)
  {
    if(wx != -1) throw std::runtime_error("What is this x1?");
    if(Sx != -1) throw std::runtime_error("What is this x2?");
    else Sx = wx = 1;
  }

  if(Sy < 0)
  {
    if(wy != -1) throw std::runtime_error("What is this y1?");
    if(Sy != -1) throw std::runtime_error("What is this y2?");
    else Sy = wy = 1;
  }

  if(Sz < 0)
  {
    if(wz != -1) throw std::runtime_error("What is this z1?");
    if(Sz != -1) throw std::runtime_error("What is this z2?");
    else Sz = wz = 1;
  }  

  
  int va[2] = {0,wz};
  int vb[2] = {0,wy};
  int vc[2] = {0,wx};

  double j = 0.0;

  for(int ia = 0; ia < 2; ia++)
    for(int ib = 0; ib < 2; ib++)
      for(int ic = 0; ic < 2; ic++)
	{
	  double a = Sz - va[ia];
	  double b = Sy - vb[ib];
	  double c = Sx - vc[ic];

	  int sgn = 1;
	  if(ia == 1) sgn *= -1;
	  if(ib == 1) sgn *= -1;
	  if(ic == 1) sgn *= -1;
			  
	  // if(R*R < a*a+b*b+c*c) continue;
	  // test up to machine precision
	  if((R*R - (a*a+b*b+c*c)) < std::nextafter(0.d,1.d)) continue;

	  double j0 = 0.0;
	  
	  j0 += R*(R*R-a*a-b*b-c*c)*(R*R-a*a-b*b-c*c)/8.0;
	  j0 -= (M_PI/2)*Sx*Sy*Sz*R*R;

	  j0 += J_s_helper(R,Sx,Sy,Sz,c,b,a);
	  j0 += J_s_helper(R,Sx,Sz,Sy,c,a,b);
	  j0 += J_s_helper(R,Sy,Sx,Sz,b,c,a);
	  j0 += J_s_helper(R,Sy,Sz,Sx,b,a,c);
	  j0 += J_s_helper(R,Sz,Sx,Sy,a,c,b);
	  j0 += J_s_helper(R,Sz,Sy,Sx,a,b,c);

	  j += sgn * j0;
	}
  
  return j;
}

void FMT_Species_Analytic_2::generateWeights()
{
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

	  // The weight for point Sx,Sy,Sz has contributions from all adjoining cells
	  // The furthest corner is Sx+1,Sy+1,Sz+1 and if this less than hsr*hsr, the volume weights are 1, and the surface weights are zero
	  // else, the nearest corner is Sx-1,Sy-1,Sz-1 (unless Sx, Sy or Sz = 0) and if this is less than hsr*hsr, then the boundary is between these limits and we must compute
	  // else, all hsd boundary is less than the nearest corner and all weights are zero.

	  double R = hsr/dx;

	  //Note: Special cases of I-sums, e.g. when one or more components of S are zero, are handled in the called functions.
	  
	  if(R*R > R2_max) {w_eta = dV;}
	  else if(R*R > R2_min)
	    for(int ix:I)
	      for(int iy:I)
		for(int iz:I)
		  {
		    w_eta  += dV*J_eta(R,Sx+ix,Sy+iy,Sz+iz,ix,iy,iz);
		    w_s    += dS*  J_s(R,Sx+ix,Sy+iy,Sz+iz,ix,iy,iz);
		    w_v[2] += dS* J_Vz(R,Sx+ix,Sy+iy,Sz+iz,ix,iy,iz);
		    w_v[1] += dS* J_Vz(R,Sx+ix,Sz+iz,Sy+iy,ix,iz,iy);
		    w_v[0] += dS* J_Vz(R,Sz+iz,Sy+iy,Sx+ix,iz,iy,ix);
		  }
	  
	  // Add in for all symmetries
	  for(int ix = 0; ix < (Sx == 0 ? 1 : 2); ix++)
	    for(int iy = 0; iy < (Sy == 0 ? 1 : 2); iy++)
	      for(int iz = 0; iz < (Sz == 0 ? 1 : 2); iz++)
		{		  
		  long pos = density_.get_PBC_Pos((1-2*ix)*Sx,(1-2*iy)*Sy,(1-2*iz)*Sz);	  
		  d_[EI()].addToWeight(pos,w_eta);
		  d_[SI()].addToWeight(pos,w_s);
		  for(int iv = 0;iv < 3;iv++)
		    d_[VI(iv)].addToWeight(pos,(iv == 0 ? (1-2*ix) : (iv == 1 ? (1-2*iy) : (1-2*iz)))*w_v[iv]);		  
		}
	}
  cout << endl;
  cout << myColor::GREEN;
  cout << "/////  Finished.  " << endl;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << myColor::RESET << endl;

}

