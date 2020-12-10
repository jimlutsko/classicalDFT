#ifndef __LUTSKO__GAUSSIAN_DENSITY__
#define __LUTSKO__GAUSSIAN_DENSITY__

#include "Density.h"
#include "Fundamental_Measures.h"

class Gaussian
{
 public:
 Gaussian(double n, double alf, double Rx, double Ry, double Rz)
   : n_(n), alf_(alf), Rx_(Rx), Ry_(Ry), Rz_(Rz), norm_(pow(alf_/M_PI,1.5)) {}

 Gaussian()
   : n_(1), alf_(100.0), Rx_(0.0), Ry_(0.0), Rz_(0.0) {}  

  double density(double rx, double ry, double rz) const
  {
    double dx = fabs(rx-Rx_); 
    double dy = fabs(ry-Ry_); 
    double dz = fabs(rz-Rz_); 

    double r2 = dx*dx+dy*dy+dz*dz;

    return n_*norm_*exp(-alf_*r2);
  }
  double density(double r[3]) const { return density(r[0], r[1], r[2]);}

  double R2max(double tol = SMALL_VALUE) const { return -log(tol/(n_*norm_))/alf_;}

  // FMT measures
  void get_measures(double rx, double ry, double rz, double hsd, FundamentalMeasures &fm) const
  {
    double dx = fabs(rx-Rx_); 
    double dy = fabs(ry-Ry_); 
    double dz = fabs(rz-Rz_); 

    double r2 = dx*dx+dy*dy+dz*dz;
    double r = sqrt(r2);
    double R = hsd/2;    

    double sqalf = sqrt(alf_);

    double ep = exp(-alf_*(r+R)*(r+R));
    double em = exp(-alf_*(r-R)*(r-R));
    
    double f = 0.5*(erf(sqalf*(r+R))-erf(sqalf*(r-R)));
    f -= (1/(4*r))*(M_2_SQRTPI/sqalf)*(em-ep);
          
    fm.eta += n_*norm_*f;

    double s = n_*norm_*(hsd*sqalf*M_2_SQRTPI/4)*(em-ep); 
    fm.s0  += s/(hsd*hsd);
    fm.s1  += s/hsd;
    fm.s2  += s;

    double v = n_*norm_*(M_2_SQRTPI/(4*sqalf*r*r))*((1-hsd*r*alf_)*em-(1+hsd*r*alf_)*ep);
    
    fm.v1[0] += (rx/r)*v/hsd;
    fm.v1[1] += (ry/r)*v/hsd;
    fm.v1[2] += (rz/r)*v/hsd;

    fm.v2[0] += (rx/r)*v;
    fm.v2[1] += (ry/r)*v;
    fm.v2[2] += (rz/r)*v;

    double T1 = n_*norm_*(M_2_SQRTPI/(4*hsd*alf_*sqalf))*(1/(r*r*r))*((hsd*r*alf_+1)*ep+(hsd*r*alf_-1)*em);
    double T2 = n_*norm_*(M_2_SQRTPI/(4*hsd*alf_*sqalf))*(1/(r*r*r*r*r))*((3-3*hsd*r*alf_+hsd*hsd*r*r*alf_*alf_)*em-(3+3*hsd*r*alf_+hsd*hsd*r*r*alf_*alf_)*ep);
    
    fm.T[0][0] += T1 + rx*rx*T2;
    fm.T[1][1] += T1 + ry*ry*T2;
    fm.T[2][2] += T1 + rz*rz*T2;

    fm.T[0][1] += rx*ry*T2;
    fm.T[0][2] += rx*rz*T2;

    fm.T[1][0] += ry*rx*T2;
    fm.T[1][2] += ry*rz*T2;

    fm.T[2][0] += rz*rx*T2;
    fm.T[2][1] += rz*ry*T2;        
  }

  
  
  double n_;
  double alf_;
  double Rx_;
  double Ry_;
  double Rz_;
 private:
  double norm_;
};
   



/**
  *  @brief Base class for the density: basically, a wrapper for the array holding the density at each lattice site
  */  

class GaussianDensity : public Density
{
 public:
  GaussianDensity(vector<Gaussian> &g, double dx, double L[])
    : Density(dx, L), gaussians_(g){ initialize_density_field(); }

 GaussianDensity(int num, double dx, double L[])
   : Density(dx, L), gaussians_(num){}

  ~GaussianDensity(){}


  void initialize_density_field(vector<Gaussian> &gaussians, double base_density = 0.0)
  {
    gaussians_.clear();
    for(auto &x: gaussians)
      gaussians_.push_back(x);
    initialize_density_field(base_density);
  }
      
  void initialize_density_field(double base_density = 0.0)
  {
    Density_.zeros();

    DFT_Vec &density = Density_.Real();
    
    density.ShiftBy(-base_density);

    for(auto &g: gaussians_)
      {
	int Nmax = 1+sqrt(g.R2max())/dx_;

	int ix = getIX(g.Rx_);
	int iy = getIY(g.Ry_);
	int iz = getIZ(g.Rz_);

	int jx;
#pragma omp parallel for  private(jx)  schedule(static)	
	for(jx = ix-Nmax; jx < ix+Nmax; jx++)
	  for(int jy = iy-Nmax; jy < iy+Nmax; jy++)
	    for(int jz = iz-Nmax; jz < iz+Nmax; jz++)
	      {
		long pos = get_PBC_Pos(jx, jy, jz);
		set(pos, get(pos) + g.density(getX(jx), getY(jy), getZ(jz)));
	      }
      }
  }

  void get_measures(double rx, double ry, double rz, double hsd, FundamentalMeasures &fm) const
  {
        for(auto &g: gaussians_)
	  g.get_measures(rx,ry,rz,hsd,fm);
  }

 protected:
  vector<Gaussian> gaussians_;
};







#endif // __LUTSKO__GAUSSIAN_DENSITY__
