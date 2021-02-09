#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>


using namespace std;

#include "Gaussian.h"

void Gaussian::set_parameters(double x, double alf, double Rx, double Ry, double Rz)
{
  x_   = x;
  alf_ = alf;
  Rx_  = Rx;
  Ry_  = Ry;
  Rz_  = Rz;

  update_cache();
}

void Gaussian::update_cache()
{  
  // the vacancy concetration: aliased against a necessary but not sufficient max
  double R    = hsd_/2;

  double z0   = min(Acut_,R);
  double z    = sqrt(alf_)*z0;
  double xmax = 1.0/(erf(z)-M_2_SQRTPI*z*exp(-z*z));

  prefactor_       = xmax*g(x_);
  dprefactor_dx_   = xmax*g1(x_);
  dprefactor_dalf_ = -M_2_SQRTPI*z*z0*z0*exp(-z*z)*xmax*xmax*g(x_);

  norm_ = pow(alf_/M_PI,1.5);

  // spatial limits concerning the influence of this Gaussian:
  // can only be determined once all parameters have been set
  // tol = A*exp(-alf*(r-R)^2) ==> r = R+sqrt(-log(tol/A)/alf)
  const double SMALL_VALUE = 1e-20;
  Rmax_ = 0.5*hsd_+sqrt(-log(SMALL_VALUE/(prefactor_*norm_))/alf_); 
  Nimage_x_ = 1; //(Rmax_ < L[0]/2 ? 0 : 1 + int(-0.5+Rmax_/L[0]));
  Nimage_y_ = 1; //(Rmax_ < L[1]/2 ? 0 : 1 + int(-0.5+Rmax_/L[1]));
  Nimage_z_ = 1; //(Rmax_ < L[2]/2 ? 0 : 1 + int(-0.5+Rmax_/L[2]));

  sq_alf_ = sqrt(alf_);
  exp_A_ = exp(-alf_*Acut_*Acut_);
  erf_A_ = erf(sq_alf_*Acut_);

  //Number of particles
  N_ = prefactor_*(erf_A_-M_2_SQRTPI*sqrt(alf_)*Acut_*exp_A_);

  // ideal-gas stuff
  double fac1 = (erf_A_-M_2_SQRTPI*sq_alf_*Acut_*exp_A_);
  double fac2 = log(prefactor_)+1.5*log(alf_)-2.5;
  double fac3 = M_2_SQRTPI*alf_*sq_alf_*Acut_*Acut_*Acut_*exp_A_;

  double fac = fac1*fac2+fac3;

  double dfac1_da = M_2_SQRTPI*Acut_*Acut_*Acut_*sq_alf_*exp_A_;
  double dfac2_da = 1.5/alf_;
  double dfac3_da = 0.5*M_2_SQRTPI*sq_alf_*Acut_*Acut_*Acut_*exp_A_*(3-2*alf_*Acut_*Acut_);

  f_id_      = prefactor_*(fac1*fac2+fac3);
  
  df_id_dx_  = dprefactor_dx_*(fac1*fac2+fac3+fac1);
  df_id_da_  = dprefactor_dalf_*(fac1*fac2+fac3+fac1);
  df_id_da_ += prefactor_*(dfac1_da*fac2+fac1*dfac2_da+dfac3_da);  
}

void Gaussian::get_parameters(double &x, double &alf, double &Rx, double &Ry, double &Rz) const
{
  alf = alf_;
  Rx  = Rx_;
  Ry  = Ry_;
  Rz  = Rz_;

  double R    = hsd_/2;    
  double z    = sqrt(alf)*min(R,Acut_);  
    //  double z    = sqrt(alf_)*hsd_/2;
  double xmax = 1.0/(erf(z)-M_2_SQRTPI*z*exp(-z*z));

  x = ginv(prefactor_/xmax);
}


double Gaussian::density(double rx, double ry, double rz) const
{
  double dx = fabs(rx-Rx_); 
  double dy = fabs(ry-Ry_); 
  double dz = fabs(rz-Rz_); 

  double r2 = dx*dx+dy*dy+dz*dz;

  if(r2 > Acut_*Acut_) return 0.0;
  
  return (r2 > Rmax_*Rmax_ ? 0.0 : prefactor_*norm_*exp(-alf_*r2));
}

void Gaussian::get_f1f2f3f4(double y, double r2, double &f1, double &f2, double &f3, double &f4) const
{
  double y2 = y*y;
  double y4 = y2*y2;
  double y6= y4*y2;
      
  if(y < 0.5)
    {
      double a1 = 2;
      double a2 = 2.0/3;
      double a3 = 2.0/15;
      double a4 = 2.0/105;
      f1 = a1; f2 = a2; f3 = a3; f4 = a4;
      double f1_old;
      int n = 1;
      do {
	f1_old = f1;

	int n2 = 2*n;
	a1 *= (y2/n2)/(n2+1); f1 += a1;
	a2 *= (y2/n2)/(n2+3); f2 += a2;
	a3 *= (y2/n2)/(n2+5); f3 += a3;
	a4 *= (y2/n2)/(n2+7); f4 += a4;

	n++;
	if(n > 20) throw std::runtime_error("Something went wrong in Gaussian::f1234");
      } while(f1 != f1_old);
      
      double ee = exp(-alf_*(r2+0.25*hsd_*hsd_));
      f1 *= ee;
      f2 *= ee;
      f3 *= ee;
      f4 *= ee;
	
    } else {
    double em = exp(-y-alf_*(r2+0.25*hsd_*hsd_));
    double ep = exp(y-alf_*(r2+0.25*hsd_*hsd_));
    f1 = (ep-em)/y;
    f2 = ((1+y)*em-(1-y)*ep)/(y*y2);
    f3 = ((3-3*y+y*y)*ep-(3+3*y+y*y)*em)/(y*y4);
    f4 = ((-15+15*y-6*y*y+y*y*y)*ep+(15+15*y+6*y*y+y*y*y)*em)/(y*y6);
  }
}


void Gaussian::get_measures(double rx, double ry, double rz, FundamentalMeasures &fm) const
{
  get_measures_internal(rx,ry,rz,prefactor_,fm);
}

void Gaussian::get_dmeasures_dX(double rx, double ry, double rz, FundamentalMeasures &dfm) const
{
  get_measures_internal(rx,ry,rz,dprefactor_dx_,dfm);
}


void Gaussian::get_measures_internal(double rx, double ry, double rz, double coeff, FundamentalMeasures &fm) const
{
  static bool firsttime = true;
  
  double dx = rx-Rx_;
  double dy = ry-Ry_;
  double dz = rz-Rz_;
  
  double R = hsd_/2;

  double r2 = dx*dx+dy*dy+dz*dz;
  double r = sqrt(r2);

  // This is to be consistent with the cutoff (i.e. with Gaussians that are not being used)
  if(r > Acut_+R) return;

  double A = 0.25*sq_alf_*alf_*hsd_*hsd_*M_2_SQRTPI;

  double y = alf_*hsd_*r;

  double e = 0;
  double s = 0;
  double v = 0;
  double T1 = 0;
  double T2 = 0;
  if(r < Acut_-R)
    {
      if(firsttime) {cout << "Region I" << endl; firsttime = false;}
      double f1 = 0;
      double f2 = 0;
      double f3 = 0;
      double f4 = 0;

      get_f1f2f3f4(y,r2,f1,f2,f3,f4);
  
      e  = -A*f1/(alf_*hsd_);
      e += 0.5*(erf(sq_alf_*(R+r))+erf(sq_alf_*(R-r)));      

      s  = A*f1;
      v  = A*f2;
      T1 = A*f2;
      T2 = A*f3;
      
    } else if(fabs(Acut_-R) < r && r < Acut_+R) {

      if(firsttime) {cout << "Region II" << endl; firsttime = false;}
    
    double z = alf_*(r*r+R*R-Acut_*Acut_);
    double C = alf_*sq_alf_*M_2_SQRTPI*R*R;
    double e1 = exp(-alf_*(R-r)*(R-r));
    double e2 = exp_A_;
    double erf_r  = erf(sq_alf_*(R-r));

    //XXX
    
    e  = -e1+(1+alf_*(Acut_-r)*(Acut_-r)-alf_*R*R)*exp_A_;
    e *= (M_2_SQRTPI/(4*sq_alf_*r));
    e += 0.5*(erf_r+erf_A_);

    s  = C*(e1-e2)/y;
    v  = C*((1-y)*e1-(1-z)*e2)/(y*y*y);
    T1 = C*((y-1)*e1-(0.5*z*z-z+1-0.5*y*y)*e2)/(y*y*y);
    T2 = C*((y*y-3*y+3)*e1-(0.5*y*y-1.5*z*z+3*z-3)*e2)/(y*y*y*y*y);
  } else if(R-Acut_ > r) {

    if(firsttime) {cout << "Region III" << endl; firsttime = false;}
    
    e = -M_2_SQRTPI*sq_alf_*Acut_*exp_A_+erf_A_;
  }

  e  *= coeff;
  s  *= coeff;
  v  *= coeff;
  T1 *= coeff;
  T2 *= coeff;

  fm.eta += e;
  
  fm.s0  += s/(hsd_*hsd_);
  fm.s1  += s/hsd_;
  fm.s2  += s;

  dx *= alf_*hsd_;
  dy *= alf_*hsd_;
  dz *= alf_*hsd_;
    
  fm.v1[0] += dx*v/hsd_;
  fm.v1[1] += dy*v/hsd_;
  fm.v1[2] += dz*v/hsd_;

  fm.v2[0] += dx*v;
  fm.v2[1] += dy*v;
  fm.v2[2] += dz*v;
    
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

void Gaussian::get_dmeasures_dR(double rx, double ry, double rz, FundamentalMeasures dfm[3]) const
{
  double dx = rx-Rx_;
  double dy = ry-Ry_;
  double dz = rz-Rz_;
  
  double R = hsd_/2;

  double r2 = dx*dx+dy*dy+dz*dz;
  double r = sqrt(r2);

  // This is to be consistent with the cutoff (i.e. with Gaussians that are not being used)
  if(r > Acut_+R) return;

  double A = 0.25*sq_alf_*alf_*hsd_*hsd_*M_2_SQRTPI;

  double y = alf_*hsd_*r;

  double e  = 0;
  double s  = 0;
  double v1 = 0;
  double v2 = 0;
  double T1 = 0;
  double T2 = 0;
  double T3 = 0;
  if(r < Acut_-R)
    {
      double f1 = 0;
      double f2 = 0;
      double f3 = 0;
      double f4 = 0;

      get_f1f2f3f4(y,r2,f1,f2,f3,f4);

      // -(dn/dr)/(alf*hsd*r)
      e = A*f2;
      
      s = (2*A/hsd_)*(f1-0.5*alf_*hsd_*hsd_*f2);
      v1 = -alf_*hsd_*A*f2;
      v2 = (2*A/hsd_)*(f2-0.5*alf_*hsd_*hsd_*f3);
      T1 = (2*A/hsd_)*f2;
      T2 = alf_*hsd_*A*f3;
      T3 = (2*A/hsd_)*(f3-0.5*alf_*hsd_*hsd_*f4);
    } else if(fabs(Acut_-R) < r && r < Acut_+R) {
    double z  = alf_*(r*r+R*R-Acut_*Acut_);
    double dz = alf_*(2*r);
    double dy = alf_*hsd_;
    double C  = alf_*sq_alf_*M_2_SQRTPI*R*R;
    double e1  = exp(-alf_*(R-r)*(R-r));
    double de1 = 2*alf_*(R-r)*exp(-alf_*(R-r)*(R-r));
    double e2 = exp_A_;
    double erf_r   = erf(sq_alf_*(R-r));
    double derf_r  = -M_2_SQRTPI*sq_alf_*e1;

    //XXX

    e  = -e1+(1+alf_*(Acut_-r)*(Acut_-r)-alf_*R*R)*exp_A_;
    e *= (M_2_SQRTPI/(4*sq_alf_*r));
    e += 0.5*(erf_r+erf_A_);
        
    e  = -(-e1+(1+alf_*(Acut_-r)*(Acut_-r)-alf_*R*R)*exp_A_)/r;
    e += -de1+(-2*alf_*(Acut_-r))*exp_A_;
    e *= (M_2_SQRTPI/(4*sq_alf_*r*r));
    e += 0.5*derf_r/r;

    e /= (-alf_*hsd_);
    
    s   = C*((1-2*alf_*r*(R-r))*e1-e2)/(r*y*y);
    v1  = -2*alf_*R*C*((1-y)*e1-(1-z)*e2)/(y*y*y);
    v2  = -C*(-dy*e1 + (1-y)*de1+dz*e2)/(y*y*y*y);
    v2 += -C*(-3*dy)*((1-y)*e1-(1-z)*e2)/(y*y*y*y*y);
	
    T1  = -C*(dy*e1+(y-1)*de1-(dz*(z-1)-dy*y)*e2)/(y*y*y*y);
    T1 += -C*(-3*dy)*((y-1)*e1-(0.5*z*z-z+1-0.5*y*y)*e2)/(y*y*y*y*y);    
    
    T2 = (2*alf_*R)*C*((y*y-3*y+3)*e1-(0.5*y*y-1.5*z*z+3*z-3)*e2)/(y*y*y*y*y);    

    T1 += T2; // correct for way calculation is organized below.

    T3  = -C*((2*dy*y-3*dy)*e1+(y*y-3*y+3)*de1-(dy*y-3*dz*z+3*dz)*e2)/(y*y*y*y*y*y);
    T3 += -C*(-5*dy)*((y*y-3*y+3)*e1-(0.5*y*y-1.5*z*z+3*z-3)*e2)/(y*y*y*y*y*y*y);

  } else if(R-Acut_ > r) {
    e = 0;
  }

  e  *= prefactor_;
  s  *= prefactor_;
  v1 *= prefactor_;
  v2 *= prefactor_;
  T1 *= prefactor_;
  T2 *= prefactor_;
  T3 *= prefactor_;  
  
  dx *= alf_*hsd_;
  dy *= alf_*hsd_;
  dz *= alf_*hsd_;  

  double yv[] = {dx,dy,dz};
  
  for(int i=0;i<3;i++)
    {
      dfm[i].eta += e*yv[i];
      
      double s1 = s*yv[i];

      dfm[i].s0  += s1/(hsd_*hsd_);
      dfm[i].s1  += s1/hsd_;
      dfm[i].s2  += s1;

      dfm[i].v1[i] += v1/hsd_;
      dfm[i].v2[i] += v1;
      
      for(int j=0;j<3;j++)
	{
	  dfm[i].v1[j] += v2*yv[i]*yv[j]/hsd_;
	  dfm[i].v2[j] += v2*yv[i]*yv[j];

	  dfm[i].T[j][j] += T1*yv[i];
	  
	  dfm[i].T[j][j] -= T2*yv[i];
	  dfm[j].T[i][j] -= T2*yv[i];
	  dfm[j].T[j][i] -= T2*yv[i];

	  for(int k=0;k<3;k++)
	    dfm[i].T[j][k] += T3*yv[i]*yv[j]*yv[k];
	}      
    }
}


void Gaussian::get_dmeasures_dAlf(double rx, double ry, double rz, FundamentalMeasures &dfm) const
{
  double dx = rx-Rx_;
  double dy = ry-Ry_;
  double dz = rz-Rz_;
  
  double R = hsd_/2;

  double r2 = dx*dx+dy*dy+dz*dz;
  double r = sqrt(r2);

  // This is to be consistent with the cutoff (i.e. with Gaussians that are not being used)
  if(r > Acut_+R) return;

  double A  = 0.25*sq_alf_*alf_*hsd_*hsd_*M_2_SQRTPI;
  double A1 = 0.25*sq_alf_*hsd_*hsd_*M_2_SQRTPI;
  double dA = 1.5*A1-(r2+0.25*hsd_*hsd_)*A;
  
  double y = alf_*hsd_*r;

  double e = 0;
  double s = 0;
  double v = 0;
  double T1 = 0;
  double T2 = 0;
  if(r < Acut_-R)
    {
      double f1 = 0;
      double f2 = 0;
      double f3 = 0;
      double f4 = 0;

      get_f1f2f3f4(y,r2,f1,f2,f3,f4);

      e  =  (A1*f1-dA*f1-A1*y*y*f2)/(alf_*hsd_);
      e += ((R+r)*exp(-alf_*(R+r)*(R+r))+(R-r)*exp(-alf_*(R-r)*(R-r)))*M_2_SQRTPI/(4*sq_alf_);
      
      s  = dA*f1+A1*y*y*f2;
      v  = dA*f2+A1*y*y*f3+A1*f2;
      T1 = dA*f2+A1*y*y*f3;
      T2 = dA*f3+A1*y*y*f4+2*A1*f3;
    } else if(fabs(Acut_-R) < r && r < Acut_+R) {

    double z   = alf_*(r*r+R*R-Acut_*Acut_);
    double C   = alf_*sq_alf_*M_2_SQRTPI*R*R;
    double dC  = 1.5*sq_alf_*M_2_SQRTPI*R*R;
    double e1  = exp(-alf_*(R-r)*(R-r));
    double de1 = -(R-r)*(R-r)*e1;
    double e2  = exp_A_;
    double de2 = -Acut_*Acut_*e2;
    double erf_r  = erf(sq_alf_*(R-r));
    double derf_r = (R-r)*e1*M_2_SQRTPI/(2*sq_alf_);
    double derf_A = Acut_*exp_A_*M_2_SQRTPI/(2*sq_alf_);
    double dexp_A = -Acut_*Acut_*exp_A_;

    double f = (M_2_SQRTPI/(4*sq_alf_*r));
    e  = -e1+(1+alf_*(Acut_-r)*(Acut_-r)-alf_*R*R)*exp_A_;
    e *= -f/(2.0*alf_);

    e += f*(-de1+((Acut_-r)*(Acut_-r)-R*R)*exp_A_);
    e += f*(1+alf_*(Acut_-r)*(Acut_-r)-alf_*R*R)*dexp_A;
    e += 0.5*(derf_r+derf_A);
        
    s  = (dC*(e1-e2)/y)+(C*(de1-de2)/y);
    s += -(C*(e1-e2)/y)/alf_;

    v  = dC*((1-y)*e1-(1-z)*e2)/(y*y*y);
    v += C*((-2+y)*e1-(-2+z)*e2)/(y*y*y*alf_);
    v += C*((1-y)*de1-(1-z)*de2)/(y*y*y);

    T1  = dC*((y-1)*e1-(0.5*z*z-z+1-0.5*y*y)*e2)/(y*y*y);
    T1 += C*((-2*y+3)*e1-(-0.5*z*z+2*z-3+0.5*y*y)*e2)/(y*y*y*alf_);
    T1 += C*((y-1)*de1-(0.5*z*z-z+1-0.5*y*y)*de2)/(y*y*y);

    T2  = dC*((y*y-3*y+3)*e1-(0.5*y*y-1.5*z*z+3*z-3)*e2)/(y*y*y*y*y);
    T2 += C*((-y*y+6*y-9)*e1-(-0.5*y*y+1.5*z*z-6*z+9)*e2)/(y*y*y*y*y*alf_);
    T2 += C*((y*y-3*y+3)*de1-(0.5*y*y-1.5*z*z+3*z-3)*de2)/(y*y*y*y*y);
  } else if(R-Acut_ > r) {
    e = M_2_SQRTPI*Acut_*Acut_*Acut_*sq_alf_*exp_A_;
  }

  e  *= prefactor_;
  s  *= prefactor_;
  v  *= prefactor_;
  T1 *= prefactor_;
  T2 *= prefactor_;

  // get contribution of dx/dalf
  get_measures_internal(rx,ry,rz,dprefactor_dalf_,dfm);
  
  dfm.eta += e;
  
  dfm.s0  += s/(hsd_*hsd_);
  dfm.s1  += s/hsd_;
  dfm.s2  += s;

  dx *= alf_*hsd_;
  dy *= alf_*hsd_;
  dz *= alf_*hsd_;
    
  dfm.v1[0] += dx*v/hsd_;
  dfm.v1[1] += dy*v/hsd_;
  dfm.v1[2] += dz*v/hsd_;

  dfm.v2[0] += dx*v;
  dfm.v2[1] += dy*v;
  dfm.v2[2] += dz*v;
    
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

