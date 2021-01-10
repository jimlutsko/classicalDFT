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

void Gaussian::set_parameters(double x, double alf, double Rx, double Ry, double Rz, double hsd, double L[])
{
  alf_ = alf;
  Rx_  = Rx;
  Ry_  = Ry;
  Rz_  = Rz;

  // the vacancy concetration: aliased against a necessary but not sufficient max
  double z    = sqrt(alf_)*hsd/2;
  double xmax = 1.0/(erf(z)-M_2_SQRTPI*z*exp(-z*z));

  prefactor_       = xmax*g(x);
  dprefactor_dx_   = xmax*g1(x);
  dprefactor_dalf_ = -M_2_SQRTPI*0.25*hsd*hsd*z*exp(-z*z)*xmax*xmax*g(x);

  norm_ = pow(alf_/M_PI,1.5);

  // spatial limits concerning the influence of this Gaussian:
  // can only be determined once all parameters have been set    
  Rmax_ = 0.5*hsd-log(SMALL_VALUE/(prefactor_*norm_))/alf_; // sqrt(R2max());
  Nimage_x_ = (Rmax_ < L[0]/2 ? 0 : 1 + int(-0.5+Rmax_/L[0]));
  Nimage_y_ = (Rmax_ < L[1]/2 ? 0 : 1 + int(-0.5+Rmax_/L[1]));
  Nimage_z_ = (Rmax_ < L[2]/2 ? 0 : 1 + int(-0.5+Rmax_/L[2]));

  Nimage_x_ = max(1,Nimage_x_);
  Nimage_y_ = max(1,Nimage_y_);
  Nimage_z_ = max(1,Nimage_z_);

}

void Gaussian::get_parameters(double &x, double &alf, double &Rx, double &Ry, double &Rz, double hsd) const
{
  alf = alf_;
  Rx  = Rx_;
  Ry  = Ry_;
  Rz  = Rz_;

  double z    = sqrt(alf_)*hsd/2;
  double xmax = 1.0/(erf(z)-M_2_SQRTPI*z*exp(-z*z));

  x = ginv(prefactor_/xmax);
}


double Gaussian::density(double rx, double ry, double rz) const
{
  double dx = fabs(rx-Rx_); 
  double dy = fabs(ry-Ry_); 
  double dz = fabs(rz-Rz_); 

  double r2 = dx*dx+dy*dy+dz*dz;

  return prefactor_*norm_*exp(-alf_*r2);
}

void Gaussian::get_f1f2f3f4(double y, double hsd, double r2, double &f1, double &f2, double &f3, double &f4) const
{
  double y2 = y*y;
  double y4 = y2*y2;
  double y6= y4*y2;
      
  if(y < 0.01)
    {
      f1 = 2+(y2/3)+(y4/60);
      f2 = (2.0/3)+(y2/15)+(y4/420);
      f3 = (2.0/15)+(y2/105)+(y4/3780);
      f4 = (2.0/105)+(y2/945)+(y4/41580);
      
      double ee = exp(-alf_*(r2+0.25*hsd*hsd));
      f1 *= ee;
      f2 *= ee;
      f3 *= ee;
      f4 *= ee;
	
    } else {
    double em = exp(-y-alf_*(r2+0.25*hsd*hsd));
    double ep = exp(y-alf_*(r2+0.25*hsd*hsd));
    f1 = (ep-em)/y;
    f2 = ((1+y)*em-(1-y)*ep)/(y*y2);
    f3 = ((3-3*y+y*y)*ep-(3+3*y+y*y)*em)/(y*y4);
    f4 = ((-15+15*y-6*y*y+y*y*y)*ep+(15+15*y+6*y*y+y*y*y)*em)/(y*y6);
  }
}


void Gaussian::get_measures(double rx, double ry, double rz, double hsd, FundamentalMeasures &fm) const
{
  double dx = rx-Rx_;
  double dy = ry-Ry_;
  double dz = rz-Rz_;
  
  double r2 = dx*dx+dy*dy+dz*dz;
  double r = sqrt(r2);

  // This is to be consistent with the cutoff (i.e. with Gaussians that are not being used)
  if(r > Rmax()) return;  

  double sqalf = sqrt(alf_);

  double A = 0.25*sqalf*alf_*hsd*hsd*M_2_SQRTPI*prefactor_;

  double y = alf_*hsd*r;

  double f1 = 0;
  double f2 = 0;
  double f3 = 0;
  double f4 = 0;

  get_f1f2f3f4(y,hsd,r2,f1,f2,f3,f4);

  //  fm.eta += -(A/(alf_*hsd))*f1+prefactor_*0.5*(erf(sqalf*(r+hsd/2))-erf(sqalf*(r-hsd/2)));
  fm.eta += (-0.25*sqalf*hsd*M_2_SQRTPI*f1+0.5*(erf(sqalf*(r+hsd/2))-erf(sqalf*(r-hsd/2))))*prefactor_;

  double s = A*f1;
  fm.s0  += s/(hsd*hsd);
  fm.s1  += s/hsd;
  fm.s2  += s;

  dx *= alf_*hsd;
  dy *= alf_*hsd;
  dz *= alf_*hsd;
    
  double v = A*f2;
    
  fm.v1[0] += dx*v/hsd;
  fm.v1[1] += dy*v/hsd;
  fm.v1[2] += dz*v/hsd;

  fm.v2[0] += dx*v;
  fm.v2[1] += dy*v;
  fm.v2[2] += dz*v;

  double T1 = A*f2;
  double T2 = A*f3;
    
  fm.T[0][0] += T1 + dx*dx*T2;
  fm.T[1][1] += T1 + dy*dy*T2;
  fm.T[2][2] += T1 + dz*dz*T2;

  fm.T[0][1] += dx*dy*T2;
  fm.T[0][2] += dx*dz*T2;

  fm.T[1][0] += dy*dx*T2;
  fm.T[1][2] += dy*dz*T2;

  fm.T[2][0] += dz*dx*T2;
  fm.T[2][1] += dz*dy*T2;

  fm.calculate_derived_quantities();      
}

void Gaussian::get_dmeasures_dX(double rx, double ry, double rz, double hsd, FundamentalMeasures &dfm) const
{
  double dx = rx-Rx_;
  double dy = ry-Ry_;
  double dz = rz-Rz_;
  
  double r2 = dx*dx+dy*dy+dz*dz;
  double r = sqrt(r2);

  // This is to be consistent with the cutoff (i.e. with Gaussians that are not being used)
  if(r > Rmax()) return;

  double sqalf = sqrt(alf_);

  double A = 0.25*sqalf*alf_*hsd*hsd*M_2_SQRTPI*dprefactor_dx_;

  double y = alf_*hsd*r;

  double f1 = 0;
  double f2 = 0;
  double f3 = 0;
  double f4 = 0;

  get_f1f2f3f4(y,hsd,r2,f1,f2,f3,f4);

  dfm.eta += -(A/(alf_*hsd))*f1+dprefactor_dx_*0.5*(erf(sqalf*(r+hsd/2))-erf(sqalf*(r-hsd/2)));

  double s = A*f1;
  dfm.s0  += s/(hsd*hsd);
  dfm.s1  += s/hsd;
  dfm.s2  += s;

  dx *= alf_*hsd;
  dy *= alf_*hsd;
  dz *= alf_*hsd;
    
  double v = A*f2;
    
  dfm.v1[0] += dx*v/hsd;
  dfm.v1[1] += dy*v/hsd;
  dfm.v1[2] += dz*v/hsd;

  dfm.v2[0] += dx*v;
  dfm.v2[1] += dy*v;
  dfm.v2[2] += dz*v;

  double T1 = A*f2;
  double T2 = A*f3;
    
  dfm.T[0][0] += T1 + dx*dx*T2;
  dfm.T[1][1] += T1 + dy*dy*T2;
  dfm.T[2][2] += T1 + dz*dz*T2;

  dfm.T[0][1] += dx*dy*T2;
  dfm.T[0][2] += dx*dz*T2;

  dfm.T[1][0] += dy*dx*T2;
  dfm.T[1][2] += dy*dz*T2;

  dfm.T[2][0] += dz*dx*T2;
  dfm.T[2][1] += dz*dy*T2;

  dfm.calculate_derived_quantities();
}



void Gaussian::get_dmeasures_dR(double rx, double ry, double rz, double hsd, FundamentalMeasures dfm[3]) const
{
  double dx = rx-Rx_;
  double dy = ry-Ry_;
  double dz = rz-Rz_;
  
  double r2 = dx*dx+dy*dy+dz*dz;
  double r = sqrt(r2);

  // This is to be consistent with the cutoff (i.e. with Gaussians that are not being used)
  if(r > Rmax()) return;
  
  double sqalf = sqrt(alf_);

  double A = 0.25*sqalf*alf_*hsd*hsd*M_2_SQRTPI*prefactor_;

  double e1 = erf(sqalf*hsd/2);
  double e2 = sqalf*M_2_SQRTPI*exp(-alf_*0.25*hsd*hsd);
    
  double y = alf_*hsd*r;

  double f1 = 0;
  double f2 = 0;
  double f3 = 0;
  double f4 = 0;
  
  get_f1f2f3f4(y,hsd,r2,f1,f2,f3,f4);

  dx *= alf_*hsd;
  dy *= alf_*hsd;
  dz *= alf_*hsd;

  double yv[] = {dx,dy,dz};

  double fac = -hsd*alf_*A;
  
  for(int i=0;i<3;i++)
    {
      dfm[i].eta += A*f2*yv[i];
      
      double s = (2*A/hsd)*(f1-0.5*alf_*hsd*hsd*f2)*yv[i];

      dfm[i].s0  += s/(hsd*hsd);
      dfm[i].s1  += s/hsd;
      dfm[i].s2  += s;

      dfm[i].v1[i] -= alf_*A*f2;
      dfm[i].v2[i] -= alf_*A*f2*hsd;
      
      for(int j=0;j<3;j++)
	{
	  dfm[i].v1[j] += (2*A/hsd)*(f2-0.5*alf_*hsd*hsd*f3)*yv[i]*yv[j]/hsd;
	  dfm[i].v2[j] += (2*A/hsd)*(f2-0.5*alf_*hsd*hsd*f3)*yv[i]*yv[j];

	  dfm[i].T[j][j] += (2*A/hsd)*f2*yv[i];
	  
	  dfm[i].T[j][j] -= alf_*hsd*A*f3*yv[i];
	  dfm[j].T[i][j] -= alf_*hsd*A*f3*yv[i];
	  dfm[j].T[j][i] -= alf_*hsd*A*f3*yv[i];

	  for(int k=0;k<3;k++)
	    dfm[i].T[j][k] += (2*A/hsd)*(f3-0.5*alf_*hsd*hsd*f4)*yv[i]*yv[j]*yv[k];
	}
      dfm[i].calculate_derived_quantities();      
    }
}


void Gaussian::get_dmeasures_dAlf(double rx, double ry, double rz, double hsd, FundamentalMeasures &dfm) const
{
  double dx = rx-Rx_;
  double dy = ry-Ry_;
  double dz = rz-Rz_;

  double r2 = dx*dx+dy*dy+dz*dz;
  double r = sqrt(r2);

  // This is to be consistent with the cutoff (i.e. with Gaussians that are not being used)
  if(r > Rmax()) return;
  
  double sqalf = sqrt(alf_);

  double A = 0.25*sqalf*alf_*hsd*hsd*M_2_SQRTPI*prefactor_;

  double z1  = sqalf*hsd/2;
  double e1 = erf(z1);
  double e2 = exp(-z1*z1);
  
  double B = 1.5 - alf_*r2 - z1*z1*e1/(e1-M_2_SQRTPI*z1*e2);

  double z    = sqrt(alf_)*hsd/2;
  double xmax = 1.0/(erf(z)-M_2_SQRTPI*z*exp(-z*z));

  
  double y = alf_*hsd*r;

  double f1 = 0;
  double f2 = 0;
  double f3 = 0;
  double f4 = 0;
  
  get_f1f2f3f4(y,hsd,r2,f1,f2,f3,f4);

  dx *= alf_*hsd;
  dy *= alf_*hsd;
  dz *= alf_*hsd;

  double yv[] = {dx,dy,dz};

  dfm.eta += (A/(hsd*alf_*alf_))*((1.5-B-y*y/(alf_*hsd*hsd))*f1 - 0.5*y*y*f2);
  dfm.eta += (B-1.5+alf_*r2+0.25*alf_*hsd*hsd)*prefactor_*0.5*(erf(sqalf*(r+hsd/2))-erf(sqalf*(r-hsd/2)))/alf_;
      
  double s = (A/alf_)*(B*f1+y*y*f2);

  dfm.s0  += s/(hsd*hsd);
  dfm.s1  += s/hsd;
  dfm.s2  += s;

  double v = (A/alf_)*((1+B)*f2+y*y*f3);
  double T1 = (A/alf_)*(B*f2+y*y*f3);
  double T2 = (A/alf_)*((2+B)*f3+f2-5*f3);
  for(int i=0;i<3;i++)
    {
      dfm.v1[i] += v*yv[i]/hsd;
      dfm.v2[i] += v*yv[i]; 

      dfm.T[i][i] += T1;
      
      for(int j=0;j<3;j++)
	dfm.T[i][j] += T2*yv[i]*yv[j];      
    }

  dfm.calculate_derived_quantities();      
}


