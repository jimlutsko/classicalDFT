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
  void get_measures(double dx, double dy, double dz, double hsd, FundamentalMeasures &fm) const
  {
    double r2 = dx*dx+dy*dy+dz*dz;
    double r = sqrt(r2);

    double sqalf = sqrt(alf_);

    double A = 0.25*sqalf*alf_*hsd*hsd*M_2_SQRTPI;
    
    double f1 = 0;
    double f2 = 0;
    double f3 = 0;

    double y = alf_*hsd*r;
    
    if(y < 0.01)
      {
	double y2 = y*y;
	double y4 = y2*y2;
	f1 = 2+(y2/3)+(y4/60);
	f2 = (2.0/3)+(y2/15)+(y4/420);
	f3 = (2.0/15)+(y2/105)+(y4/3780);

	double ee = exp(-alf_*(r2+0.25*hsd*hsd));
	f1 *= ee;
	f2 *= ee;
	f3 *= ee;
	
      } else {
      double em = exp(-y-alf_*(r2+0.25*hsd*hsd));
      double ep = exp(y-alf_*(r2+0.25*hsd*hsd));
      f1 = (ep-em)/y;
      f2 = ((1+y)*em-(1-y)*ep)/(y*y*y);
      f3 = ((3-3*y+y*y)*ep-(3+3*y+y*y)*em)/(y*y*y*y*y);
    }

    fm.eta += -n_*(A/(alf_*hsd))*f1+n_*0.5*(erf(sqalf*(r+hsd/2))-erf(sqalf*(r-hsd/2)));

    double s = n_*A*f1;
    fm.s0  += s/(hsd*hsd);
    fm.s1  += s/hsd;
    fm.s2  += s;

    dx *= alf_*hsd;
    dy *= alf_*hsd;
    dz *= alf_*hsd;
    
    double v = n_*A*f2;
    
    fm.v1[0] += dx*v/hsd;
    fm.v1[1] += dy*v/hsd;
    fm.v1[2] += dz*v/hsd;

    fm.v2[0] += dx*v;
    fm.v2[1] += dy*v;
    fm.v2[2] += dz*v;

    double T1 = n_*A*f2;
    double T2 = n_*A*f3;
    
    fm.T[0][0] += T1 + dx*dx*T2;
    fm.T[1][1] += T1 + dy*dy*T2;
    fm.T[2][2] += T1 + dz*dz*T2;

    fm.T[0][1] += dx*dy*T2;
    fm.T[0][2] += dx*dz*T2;

    fm.T[1][0] += dy*dx*T2;
    fm.T[1][2] += dy*dz*T2;

    fm.T[2][0] += dz*dx*T2;
    fm.T[2][1] += dz*dy*T2;

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
		set(pos, get(pos) + g.density(getX(jx),getY(jy),getZ(jz)));
	      }
      }
  }

  void get_measures(long pos, double hsd, FundamentalMeasures &fm) const
  {
    long ix,iy,iz;
    cartesian(pos,ix,iy,iz);
    get_measures(getX(ix), getY(iy), getZ(iz), hsd, fm);
  }
  
  void get_measures(double rx, double ry, double rz, double hsd, FundamentalMeasures &fm) const
  {
    for(auto &g: gaussians_)
      {
	// find closest image
	double dx = rx - g.Rx_; while(dx > L_[0]/2) dx -= L_[0]; while(dx < -L_[0]/2) dx += L_[0];
	double dy = ry - g.Ry_; while(dy > L_[1]/2) dy -= L_[1]; while(dy < -L_[1]/2) dy += L_[1];
	double dz = rz - g.Rz_; while(dz > L_[2]/2) dz -= L_[2]; while(dz < -L_[2]/2) dz += L_[2]; 
	
	// take additional images as necessary: if Rmax > L/2, we need more images
	double Rmax = sqrt(g.R2max());
	int Nimage_x = (Rmax < L_[0]/2 ? 0 : 1 + int(-0.5+Rmax/L_[0]));
	int Nimage_y = (Rmax < L_[1]/2 ? 0 : 1 + int(-0.5+Rmax/L_[1]));
	int Nimage_z = (Rmax < L_[2]/2 ? 0 : 1 + int(-0.5+Rmax/L_[2]));

	
	//	double Rmax = sqrt(g.R2max());
	//	int Nimage_x = std::max(0, int(2*Rmax/L_[0]));
	//	int Nimage_y = std::max(0, int(2*Rmax/L_[1]));
	//	int Nimage_z = std::max(0, int(2*Rmax/L_[2]));

	for(int imx = -Nimage_x; imx <= Nimage_x; imx++)
	  for(int imy = -Nimage_y; imy <= Nimage_y; imy++)
	    for(int imz = -Nimage_z; imz <= Nimage_z; imz++)
	      g.get_measures(dx+imx*L_[0],dy+imy*L_[1],dz+imz*L_[2],hsd,fm);

      }
  }

 protected:
  vector<Gaussian> gaussians_;
};








#endif // __LUTSKO__GAUSSIAN_DENSITY__
