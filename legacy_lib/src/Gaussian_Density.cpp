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

void GaussianDensity::get_measures(double rx1, double ry1, double rz1, double hsd, FundamentalMeasures &fm) const
{
  // maximum range for a non-zero contribution is A+R
  // so we need max (2A+2R)/L images
  
  for(auto &g: gaussians_)
    {
      double rx = rx1 - g.Rx();
      double ry = ry1 - g.Ry();
      double rz = rz1 - g.Rz();
      
      // find closest image
      while(rx > L_[0]/2) rx -= L_[0]; while(rx < -L_[0]/2) rx += L_[0];
      while(ry > L_[1]/2) ry -= L_[1]; while(ry < -L_[1]/2) ry += L_[1];
      while(rz > L_[2]/2) rz -= L_[2]; while(rz < -L_[2]/2) rz += L_[2];

       // sum over all contributing images
      //            for(int imx = -g.get_Nimage_x(); imx <= g.get_Nimage_x(); imx++)
      //      	for(int imy = -g.get_Nimage_y(); imy <= g.get_Nimage_y(); imy++)
      //      	  for(int imz = -g.get_Nimage_z(); imz <= g.get_Nimage_z(); imz++)	      
      //      	    g.get_measures(rx+imx*L_[0],ry+imy*L_[1],rz+imz*L_[2],fm);

      
      double R2max = (g.A()+hsd/2)*(g.A()+hsd/2);
      
      for(double Rx = rx; Rx*Rx+ry*ry+rz*rz <  R2max; Rx += L_[0])
	{
	  for(double Ry = ry; Rx*Rx+Ry*Ry+rz*rz < R2max; Ry += L_[1])
	    {
	      for(double Rz = rz; Rx*Rx+Ry*Ry+Rz*Rz < R2max; Rz += L_[2])
		g.get_measures(Rx+g.Rx(),Ry+g.Ry(),Rz+g.Rz(),fm);
	      for(double Rz = rz-L_[2]; Rx*Rx+Ry*Ry+Rz*Rz < R2max; Rz -= L_[2])
		g.get_measures(Rx+g.Rx(),Ry+g.Ry(),Rz+g.Rz(),fm);
	    }
	  for(double Ry = ry-L_[1]; Rx*Rx+Ry*Ry+rz*rz < R2max; Ry -= L_[1])
	    {
	      for(double Rz = rz; Rx*Rx+Ry*Ry+Rz*Rz < R2max; Rz += L_[2])
		g.get_measures(Rx+g.Rx(),Ry+g.Ry(),Rz+g.Rz(),fm);
	      for(double Rz = rz-L_[2]; Rx*Rx+Ry*Ry+Rz*Rz < R2max; Rz -= L_[2])
		g.get_measures(Rx+g.Rx(),Ry+g.Ry(),Rz+g.Rz(),fm);
	    }
	}
      for(double Rx = rx-L_[0]; Rx*Rx+ry*ry+rz*rz <  R2max; Rx -= L_[0])
	{
	  for(double Ry = ry; Rx*Rx+Ry*Ry+rz*rz < R2max; Ry += L_[1])
	    {
	      for(double Rz = rz; Rx*Rx+Ry*Ry+Rz*Rz < R2max; Rz += L_[2])
		g.get_measures(Rx+g.Rx(),Ry+g.Ry(),Rz+g.Rz(),fm);
	      for(double Rz = rz-L_[2]; Rx*Rx+Ry*Ry+Rz*Rz < R2max; Rz -= L_[2])
		g.get_measures(Rx+g.Rx(),Ry+g.Ry(),Rz+g.Rz(),fm);
	    }
	  for(double Ry = ry-L_[1]; Rx*Rx+Ry*Ry+rz*rz < R2max; Ry -= L_[1])
	    {
	      for(double Rz = rz; Rx*Rx+Ry*Ry+Rz*Rz < R2max; Rz += L_[2])
		g.get_measures(Rx+g.Rx(),Ry+g.Ry(),Rz+g.Rz(),fm);
	      for(double Rz = rz-L_[2]; Rx*Rx+Ry*Ry+Rz*Rz < R2max; Rz -= L_[2])
		g.get_measures(Rx+g.Rx(),Ry+g.Ry(),Rz+g.Rz(),fm);
	    }
	}
    }   
         
}

void GaussianDensity::get_dmeasures_for_gaussian(int igaussian, double rx, double ry, double rz, double hsd, FundamentalMeasures dfm[5]) const
{
  const Gaussian &g = gaussians_[igaussian];

  // find closest image
  while(rx-g.Rx() > L_[0]/2) rx -= L_[0]; while(rx-g.Rx() < -L_[0]/2) rx += L_[0];
  while(ry-g.Ry() > L_[1]/2) ry -= L_[1]; while(ry-g.Ry() < -L_[1]/2) ry += L_[1];
  while(rz-g.Rz() > L_[2]/2) rz -= L_[2]; while(rz-g.Rz() < -L_[2]/2) rz += L_[2];
	  
  // sum over all contributing images
  for(int imx = -g.get_Nimage_x(); imx <= g.get_Nimage_x(); imx++)
    for(int imy = -g.get_Nimage_y(); imy <= g.get_Nimage_y(); imy++)
      for(int imz = -g.get_Nimage_z(); imz <= g.get_Nimage_z(); imz++)
	{
	  g.get_dmeasures_dX(rx+imx*L_[0],ry+imy*L_[1],rz+imz*L_[2],dfm[0]);
	  g.get_dmeasures_dAlf(rx+imx*L_[0],ry+imy*L_[1],rz+imz*L_[2],dfm[1]);
	  g.get_dmeasures_dR(rx+imx*L_[0],ry+imy*L_[1],rz+imz*L_[2],dfm+2);
	}
}


  double GaussianDensity::FMF(double w0, double r0, vector<double> &x, vector<double> &w, DFT_Vec &dF1) const
  {
    double F = 0;

    int XMAX = 3;
    int YMAX = 3;
    int ZMAX = 3;

    vector<double> dF(5*gaussians_.size(),0.0);

    for(int l=0; l<gaussians_.size(); l++)
      {
	double al     = gaussians_[l].alf();
	double xl     = gaussians_[l].prefactor();
	double dxl_dy = gaussians_[l].dprefactor_dx();
	double dxl_da = gaussians_[l].dprefactor_dalf();
      
	double all = al/2;
	double sall = sqrt(all);

	double ep0 = exp(-all*r0*r0);
	double dFdx = 0;
	double arg = (erf(sall*r0)-M_2_SQRTPI*sall*r0*ep0);
      
	F    += 0.5*w0*xl*xl*arg;
	dFdx += w0*xl*arg;

	double sum1 = 0.0;
	for(int n=0;n<x.size();n++)
	  {
	    arg = w[n]*x[n]*x[n]*exp(-all*x[n]*x[n]);
	    F    += xl*xl*M_2_SQRTPI*sall*all*arg;
	    dFdx +=  2*xl*M_2_SQRTPI*sall*all*arg;
	    sum1 += (3-2*all*x[n]*x[n])*arg;
	  }
      
	dF[5*l+0] += dxl_dy*dFdx;
	dF[5*l+1] += dxl_da*dFdx;
	dF[5*l+1] += 0.25*xl*xl*M_2_SQRTPI*sall*(w0*r0*r0*r0*ep0 + sum1);	      
      
	for(int m=l; m<gaussians_.size(); m++)
	  {
	    double Rx = gaussians_[l].Rx() - gaussians_[m].Rx(); while(Rx > L_[0]/2) Rx -= L_[0]; while(Rx < -L_[0]/2) Rx += L_[0];
	    double Ry = gaussians_[l].Ry() - gaussians_[m].Ry(); while(Ry > L_[1]/2) Ry -= L_[1]; while(Ry < -L_[1]/2) Ry += L_[1];
	    double Rz = gaussians_[l].Rz() - gaussians_[m].Rz(); while(Rz > L_[2]/2) Rz -= L_[2]; while(Rz < -L_[2]/2) Rz += L_[2];

	    for(int ix = -XMAX; ix <= XMAX; ix++)
	      for(int iy = -YMAX; iy <= YMAX; iy++)
		for(int iz = -ZMAX; iz <= ZMAX; iz++)
		  {
		    if(m == l && ix == 0 && iy == 0 && iz == 0) continue;
		    double fac = 1;
		    if(m == l) fac = 0.5;
		  
		    double rx = Rx + ix*L_[0];
		    double ry = Ry + iy*L_[1];
		    double rz = Rz + iz*L_[2];
		  
		    double Rlm = sqrt(rx*rx+ry*ry+rz*rz);
	  
		    double am     = gaussians_[m].alf();
		    double xm     = gaussians_[m].prefactor();
		    double dxm_dy = gaussians_[m].dprefactor_dx();
		    double dxm_da = gaussians_[m].dprefactor_dalf();
		  
		    double alm  = al*am/(al+am);
		    double salm = sqrt(alm);

		    double er = 0.5*erf(salm*(r0-Rlm))+0.5*erf(salm*(r0+Rlm));
		    double em = exp(-alm*(r0-Rlm)*(r0-Rlm));
		    double ep = exp(-alm*(r0+Rlm)*(r0+Rlm));
		  		  
		    F += fac*w0*xl*xm*er;
		    F += fac*w0*xl*xm*(ep-em)*M_2_SQRTPI/(4*salm*Rlm);

		    dF[5*l+0] += fac*dxl_dy*w0*xm*(er+(ep-em)*M_2_SQRTPI/(4*salm*Rlm));
		    dF[5*m+0] += fac*dxm_dy*w0*xl*(er+(ep-em)*M_2_SQRTPI/(4*salm*Rlm));

		    dF[5*l+1] += fac*dxl_da*w0*xm*(er+(ep-em)*M_2_SQRTPI/(4*salm*Rlm));
		    dF[5*m+1] += fac*dxm_da*w0*xl*(er+(ep-em)*M_2_SQRTPI/(4*salm*Rlm));		  
		  
		    dF[5*l+1] += fac*(xl*xm/(al*al))*w0*M_2_SQRTPI*(salm/(8*Rlm))*((1-2*alm*r0*(Rlm-r0))*em-(1+2*alm*r0*(Rlm+r0))*ep);
		    dF[5*m+1] += fac*(xl*xm/(am*am))*w0*M_2_SQRTPI*(salm/(8*Rlm))*((1-2*alm*r0*(Rlm-r0))*em-(1+2*alm*r0*(Rlm+r0))*ep);

		    dF[5*l+2] += fac*w0*rx*xl*xm*(M_2_SQRTPI/(4*Rlm*Rlm*Rlm*salm))*((1-2*alm*Rlm*r0)*em-(1+2*alm*Rlm*r0)*ep);
		    dF[5*m+2] -= fac*w0*rx*xl*xm*(M_2_SQRTPI/(4*Rlm*Rlm*Rlm*salm))*((1-2*alm*Rlm*r0)*em-(1+2*alm*Rlm*r0)*ep);

		    dF[5*l+3] += fac*w0*ry*xl*xm*(M_2_SQRTPI/(4*Rlm*Rlm*Rlm*salm))*((1-2*alm*Rlm*r0)*em-(1+2*alm*Rlm*r0)*ep);
		    dF[5*m+3] -= fac*w0*ry*xl*xm*(M_2_SQRTPI/(4*Rlm*Rlm*Rlm*salm))*((1-2*alm*Rlm*r0)*em-(1+2*alm*Rlm*r0)*ep);

		    dF[5*l+4] += fac*w0*rz*xl*xm*(M_2_SQRTPI/(4*Rlm*Rlm*Rlm*salm))*((1-2*alm*Rlm*r0)*em-(1+2*alm*Rlm*r0)*ep);
		    dF[5*m+4] -= fac*w0*rz*xl*xm*(M_2_SQRTPI/(4*Rlm*Rlm*Rlm*salm))*((1-2*alm*Rlm*r0)*em-(1+2*alm*Rlm*r0)*ep);		  

		    double sum1 = 0;
		    double sum2 = 0;
		    double sum3 = 0;
		    for(int n=0;n<x.size();n++)
		      {
			double enm = exp(-alm*(x[n]-Rlm)*(x[n]-Rlm));
			double enp = exp(-alm*(x[n]+Rlm)*(x[n]+Rlm));

			sum1 += w[n]*x[n]*(enm-enp);
			sum2 += w[n]*x[n]*((1-2*alm*(Rlm-x[n])*(Rlm-x[n]))*enm-(1-2*alm*(Rlm+x[n])*(Rlm+x[n]))*enp);
			sum3 += w[n]*x[n]*((1+2*alm*Rlm*(Rlm+x[n]))*enp-(1+2*alm*Rlm*(Rlm-x[n]))*enm);
		      }
		    F += fac*xl*xm*(M_2_SQRTPI*salm/(2*Rlm))*sum1;
		  
		    dF[5*l+0] += fac*dxl_dy*xm*M_2_SQRTPI*(salm/(2*Rlm))*sum1;
		    dF[5*m+0] += fac*dxm_dy*xl*M_2_SQRTPI*(salm/(2*Rlm))*sum1;
		  
		    dF[5*l+1] += fac*dxl_da*xm*M_2_SQRTPI*(salm/(2*Rlm))*sum1;
		    dF[5*m+1] += fac*dxm_da*xl*M_2_SQRTPI*(salm/(2*Rlm))*sum1;

		    dF[5*l+1] += fac*0.25*(xl*xm/(al*al))*M_2_SQRTPI*alm*(salm/Rlm)*sum2;
		    dF[5*m+1] += fac*0.25*(xl*xm/(am*am))*M_2_SQRTPI*alm*(salm/Rlm)*sum2;

		    dF[5*l+2] += fac*rx*xl*xm*M_2_SQRTPI*(salm/(2*Rlm*Rlm*Rlm))*sum3;
		    dF[5*m+2] -= fac*rx*xl*xm*M_2_SQRTPI*(salm/(2*Rlm*Rlm*Rlm))*sum3;

		    dF[5*l+3] += fac*ry*xl*xm*M_2_SQRTPI*(salm/(2*Rlm*Rlm*Rlm))*sum3;
		    dF[5*m+3] -= fac*ry*xl*xm*M_2_SQRTPI*(salm/(2*Rlm*Rlm*Rlm))*sum3;

		    dF[5*l+4] += fac*rz*xl*xm*M_2_SQRTPI*(salm/(2*Rlm*Rlm*Rlm))*sum3;
		    dF[5*m+4] -= fac*rz*xl*xm*M_2_SQRTPI*(salm/(2*Rlm*Rlm*Rlm))*sum3;		  


		  }
	  }
      
      }

    for(unsigned pos=0;pos<dF.size();pos++)
      dF1.set(pos,dF[pos]);
  
    return F;
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
