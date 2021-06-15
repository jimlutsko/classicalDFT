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


EOS_Species::EOS_Species(Density& density, double hsd, int seq): Species(density,0,seq), hsd_(hsd), moments_(5)
{
  long Nx = density_.Nx();
  long Ny = density_.Ny();
  long Nz = density_.Nz();

  for(FMT_Weighted_Density &d: moments_)
    d.initialize(Nx, Ny, Nz);

  generateWeights(hsd, moments_);

  for(FMT_Weighted_Density &d: moments_)
    d.transformWeights();

}


void EOS_Species::report(FMT_Species &sp)
{
  double dmax = 0;
  double r2   = 0;
  long   pmax = -1;
  double eta  = 0;
  
  for(long pos = 0; pos < density_.Ntot(); pos++)
  //  long pos = density_.pos(31,31,31);
    {
      double I0  = moments_[0].getDensity(pos);
      double I1x = moments_[1].getDensity(pos);
      double I1y = moments_[2].getDensity(pos);
      double I1z = moments_[3].getDensity(pos);
      double I2  = moments_[4].getDensity(pos);

      double s  = sp.getS(pos);
      double vx = sp.getV(0,pos);
      double vy = sp.getV(1,pos);
      double vz = sp.getV(2,pos);
      
      double dd = s*(0.25*hsd_*hsd_+(I1x*I1x+I1y*I1y+I1z*I1z)/(I0*I0));
      dd -= hsd_*(I1x*vx+I1y*vy+I1z*vz)/(I0);
      dd /= (hsd_*hsd_*M_PI);

      dd *= (M_PI/6)*hsd_*hsd_*hsd_;
      
      /*
      double I12 = I1x*I1x+I1y*I1y+I1z*I1z;
      double I1v = I1x*vx +I1y*vy +I1z*vz;

      double dd = (0.25*hsd_*hsd_-(I2/I0)+2*I12/(I0*I0))*(s/I0)-hsd_*I1v/(I0*I0);

      dd *= 5/(3*hsd_);
      dd *= I0;

      */

      
      //      dd = s*s-(vx*vx+vy*vy+vz*vz);
      //      dd /= ((hsd_*hsd_*M_PI)*(hsd_*hsd_*M_PI));

      //      dd *= I0; // effective eta
      if(dd > dmax) { dmax = dd;  pmax = pos; eta = I0;}      
      
      /*
      
      double dd = I2 - (I1x*I1x+I1y*I1y+I1z*I1z)/I0;
      // dd /= I0;
      double ss = dd;
      dd *= 20/(3*hsd_*hsd_);
      
      if(dd > dmax) { dmax = dd; r2 = ss/I0; pmax = pos; eta = I0;}
      */
      
    }
  long i,j,k;
  density_.cartesian(pmax, i, j, k);
  cout << "etamax = " << dmax << " <r^2> = " <<  r2 << " eta = " << eta << " pos = " << i << " " << j << " " << k << endl;
}


void EOS_Species::generateWeights(double hsd, vector<FMT_Weighted_Density> &moment_weights)
{
  cout << "EOS_Species: Generating weights" << endl;
  
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
  double dx0 = dx;
  double dy0 = dy;
  double dz0 = dz;
  
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
	  
	  double w_I0    = 0.0;
	  double w_I2    = 0.0;
	  double w_I1[3] = {0.0,0.0,0.0};

	  // The weight for point Sx,Sy,Sz has contributions from all adjoining cells
	  // The furthest corner is Sx+1,Sy+1,Sz+1 and if this less than hsr*hsr, the volume weights are 1, and the surface weights are zero
	  // else, the nearest corner is Sx-1,Sy-1,Sz-1 (unless Sx, Sy or Sz = 0) and if this is less than hsr*hsr, then the boundary is between these limits and we must compute
	  // else, all hsd boundary is less than the nearest corner and all weights are zero.

	  //Note: Special cases of I-sums, e.g. when one or more components of S are zero, are handled in the called functions.
	  
	  if(R*R > R2_max) {w_I0 = dV; w_I1[0] = Sx*dx0*dV; w_I1[1] = Sy*dy0*dV; w_I1[2] = Sz*dz0*dV; w_I2 = 0.5*(((dx0*dx0+dy0*dy0+dz0*dz0)/3)+2*Sx*Sx*dx0*dx0+2*Sy*Sy*dy0*dy0+2*Sz*Sz*dz0*dz0)*dV;}
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
				
				w_I0    += dV*sgn*(FMT_Helper::G_eta(R,Vx,Vy,Vz,Tx,Ty,Tz) - FMT_Helper::G_eta(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));

				w_I1[0] += dV*dx0*px*sgn*(FMT_Helper::G_I1x(R,Vx,Vy,Vz,Tx,Ty,Tz) - FMT_Helper::G_I1x(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				w_I1[1] += dV*dy0*py*sgn*(FMT_Helper::G_I1x(R,Vy,Vx,Vz,Ty,Tx,Tz) - FMT_Helper::G_I1x(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz));
				w_I1[2] += dV*dz0*pz*sgn*(FMT_Helper::G_I1x(R,Vz,Vy,Vx,Tz,Ty,Tx) - FMT_Helper::G_I1x(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx));

				w_I2    += dV*sgn*dx0*dx0*px*px*(FMT_Helper::G_I2(R,Vx,Vy,Vz,Tx,Ty,Tz)   - FMT_Helper::G_I2(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				w_I2    += dV*sgn*dy0*dy0*py*py*(FMT_Helper::G_I2(R,Vy,Vx,Vz,Ty,Tx,Tz)   - FMT_Helper::G_I2(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz));
				w_I2    += dV*sgn*dz0*dz0*pz*pz*(FMT_Helper::G_I2(R,Vz,Vy,Vx,Tz,Ty,Tx)   - FMT_Helper::G_I2(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx));

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
		  moment_weights[0].addToWeight(pos,w_I0);
		  moment_weights[1].addToWeight(pos,(1-2*ix)*w_I1[0]);
		  moment_weights[2].addToWeight(pos,(1-2*iy)*w_I1[1]);
		  moment_weights[3].addToWeight(pos,(1-2*iz)*w_I1[2]);		  
		  moment_weights[4].addToWeight(pos,w_I2);
		  if(isnan(moment_weights[0].getWeight(pos)))
		    {
		      cout << ix << " " << iy << " " << iz << " " << Sx << " " << Sy << " " << Sz << endl;
		      throw std::runtime_error("Found NAN in EOS_Species::getWeight");
		    }
		}		  
	}

  cout << endl;
  cout << myColor::GREEN;
  cout << "/////  Finished.  " << endl;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << myColor::RESET << endl;
}
