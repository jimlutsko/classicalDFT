#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>


using namespace std;

#include "Gaussian_Density.h"


void GaussianDensity::get_dN(int ig, double &dN_dx, double &dN_dalf) const
{
  const Gaussian &g = gaussians_[ig];
  dN_dx = g.dprefactor_dx();
  dN_dalf = g.dprefactor_dalf();
}


void GaussianDensity::get_images(double Rx1, double Ry1, double Rz1, double R2max, vector<vector<double>> &images) const
{
  // find nearest image
  
  while(Rx1 > L_[0]/2) Rx1 -= L_[0]; while(Rx1 < -L_[0]/2) Rx1 += L_[0];
  while(Ry1 > L_[1]/2) Ry1 -= L_[1]; while(Ry1 < -L_[1]/2) Ry1 += L_[1];
  while(Rz1 > L_[2]/2) Rz1 -= L_[2]; while(Rz1 < -L_[2]/2) Rz1 += L_[2];

  for(double Rx = Rx1; Rx*Rx+Ry1*Ry1+Rz1*Rz1 <  R2max; Rx += L_[0])
    {
      for(double Ry = Ry1; Rx*Rx+Ry*Ry+Rz1*Rz1 < R2max; Ry += L_[1])
	{
	  for(double Rz = Rz1; Rx*Rx+Ry*Ry+Rz*Rz < R2max; Rz += L_[2])
	    {
	      vector<double> R{Rx,Ry,Rz};
	      images.push_back(R);
	    }
	  for(double Rz = Rz1-L_[2]; Rx*Rx+Ry*Ry+Rz*Rz < R2max; Rz -= L_[2])
	    {
	      vector<double> R{Rx,Ry,Rz};
	      images.push_back(R);
	    }	    
	}
      for(double Ry = Ry1-L_[1]; Rx*Rx+Ry*Ry+Rz1*Rz1 < R2max; Ry -= L_[1])
	{
	  for(double Rz = Rz1; Rx*Rx+Ry*Ry+Rz*Rz < R2max; Rz += L_[2])
	    {
	      vector<double> R{Rx,Ry,Rz};
	      images.push_back(R);
	    }	    
	  for(double Rz = Rz1-L_[2]; Rx*Rx+Ry*Ry+Rz*Rz < R2max; Rz -= L_[2])
	    {
	      vector<double> R{Rx,Ry,Rz};
	      images.push_back(R);
	    }	    
	}
    }
  for(double Rx = Rx1-L_[0]; Rx*Rx+Ry1*Ry1+Rz1*Rz1 <  R2max; Rx -= L_[0])
    {
      for(double Ry = Ry1; Rx*Rx+Ry*Ry+Rz1*Rz1 < R2max; Ry += L_[1])
	{
	  for(double Rz = Rz1; Rx*Rx+Ry*Ry+Rz*Rz < R2max; Rz += L_[2])
	    {
	      vector<double> R{Rx,Ry,Rz};
	      images.push_back(R);
	    }	    
	  for(double Rz = Rz1-L_[2]; Rx*Rx+Ry*Ry+Rz*Rz < R2max; Rz -= L_[2])
	    {
	      vector<double> R{Rx,Ry,Rz};
	      images.push_back(R);
	    }	    
	}
      for(double Ry = Ry1-L_[1]; Rx*Rx+Ry*Ry+Rz1*Rz1 < R2max; Ry -= L_[1])
	{
	  for(double Rz = Rz1; Rx*Rx+Ry*Ry+Rz*Rz < R2max; Rz += L_[2])
	    {
	      vector<double> R{Rx,Ry,Rz};
	      images.push_back(R);
	    }	    
	  for(double Rz = Rz1-L_[2]; Rx*Rx+Ry*Ry+Rz*Rz < R2max; Rz -= L_[2])
	    {
	      vector<double> R{Rx,Ry,Rz};
	      images.push_back(R);
	    }	    
	}
    }  
}

void GaussianDensity::get_measures(double rx1, double ry1, double rz1, double hsd, FundamentalMeasures &fm) const
{
  // maximum range for a non-zero contribution is A+R
  // so we need max (2A+2R)/L images
  
  for(auto &g: gaussians_)
    {
      double rx = rx1 - g.Rx();
      double ry = ry1 - g.Ry();
      double rz = rz1 - g.Rz();

      double R2max = (g.A()+hsd/2)*(g.A()+hsd/2);

      vector<vector<double>> images;
      get_images(rx, ry, rz, R2max, images);

      for(auto &R: images)
	g.get_measures(R[0]+g.Rx(),R[1]+g.Ry(),R[2]+g.Rz(),fm);
    }   
         
}

void GaussianDensity::get_dmeasures_for_gaussian(int igaussian, double rx, double ry, double rz, double hsd, FundamentalMeasures dfm[5]) const
{
  const Gaussian &g = gaussians_[igaussian];

  rx -= g.Rx();
  ry -= g.Ry();
  rz -= g.Rz();  

  double R2max = (g.A()+hsd/2)*(g.A()+hsd/2);

  vector<vector<double>> images;
  get_images(rx, ry, rz, R2max, images);

  for(auto &R: images)  
    {
      g.get_dmeasures_dX(R[0]+g.Rx(),R[1]+g.Ry(),R[2]+g.Rz(),dfm[0]);
      g.get_dmeasures_dAlf(R[0]+g.Rx(),R[1]+g.Ry(),R[2]+g.Rz(),dfm[1]);
      g.get_dmeasures_dR(R[0]+g.Rx(),R[1]+g.Ry(),R[2]+g.Rz(),dfm+2);
    }
}

void GaussianDensity::fill_discrete_gaussian_field()
{
  discrete_gaussian_field.zeros();
  
  for(Gaussian &g: gaussians_)
    {
      int ix = getIX(g.Rx());
      int iy = getIY(g.Ry());
      int iz = getIZ(g.Rz());

      double N_max_x = g.Rmax()/dx_;
      double N_max_y = g.Rmax()/dy_;
      double N_max_z = g.Rmax()/dz_;

      for(int jx = ix-N_max_x; jx <= ix+N_max_x;jx++)
	for(int jy = iy-N_max_y; jy <= iy+N_max_y;jy++)
	  for(int jz = iz-N_max_z; jz <= iz+N_max_z;jz++)
	    {
	      double Rx = getX(jx);
	      double Ry = getY(jy);
	      double Rz = getZ(jz);

	      double d = g.density(Rx,Ry,Rz);
	      long pos = get_PBC_Pos(jx,jy,jz);

	      discrete_gaussian_field.Real().IncrementBy(pos,d);
	    }	        
    }  
}
