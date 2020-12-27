#ifndef __LUTSKO__GAUSSIAN_DENSITY__
#define __LUTSKO__GAUSSIAN_DENSITY__

#include "Density.h"
#include "Fundamental_Measures.h"



class Gaussian
{
 public:
  // Parameters are always set by calling set_parameters
  Gaussian()  {}
  ~Gaussian() {}

  // Setting the parameters
  // Note that the parameter y is an alias controlling the vacancy concentration
  void set_parameters(double x, double alf, double Rx, double Ry, double Rz, double hsd, double L[]);

  static double get_prefac_alias(double prefac, double alf, double hsd) 
  {
    double z    = sqrt(alf)*hsd/2;
    double xmax = 1.0/(erf(z)-M_2_SQRTPI*z*exp(-z*z));

    Gaussian g;
    return ginv(prefac/xmax);
  }
  
  // maximum range of this gaussian
  double R2max(double tol = SMALL_VALUE) const { return -log(tol/(prefactor_*norm_))/alf_;}
  double Rmax() const { return Rmax_;}

  
  double density(double rx, double ry, double rz) const
  {
    double dx = fabs(rx-Rx_); 
    double dy = fabs(ry-Ry_); 
    double dz = fabs(rz-Rz_); 

    double r2 = dx*dx+dy*dy+dz*dz;

    return prefactor_*norm_*exp(-alf_*r2);
  }
  double density(double r[3]) const { return density(r[0], r[1], r[2]);}


  // FMT measures
  void get_measures(double rx, double ry, double rz, double hsd, FundamentalMeasures &fm) const;
  void get_dmeasures_dAlf(double rx, double ry, double rz, double hsd, FundamentalMeasures &dfm) const;
  void get_dmeasures_dR(double rx, double ry, double rz, double hsd, FundamentalMeasures dfm[3]) const;
  void get_dmeasures_dX(double rx, double ry, double rz, double hsd, FundamentalMeasures &dfm) const;

  int get_Nimage_x() const { return Nimage_x_;}
  int get_Nimage_y() const { return Nimage_y_;}
  int get_Nimage_z() const { return Nimage_z_;}

  double Rx() const { return Rx_;}
  double Ry() const { return Ry_;}
  double Rz() const { return Rz_;}


protected:
  // these are the alias function, its inverse and its derivative
  // Maybe better if these were members that could be over-ridden
  static double    g(double x) { return 1.0-exp(-x*x);}
  static double ginv(double z) { return sqrt(fabs(log(1.0-z)));}
  static double   g1(double x) { return 2*x*exp(-x*x);}

  void get_f1f2f3f4(double y, double hsd, double r2, double &f1, double &f2, double &f3, double &f4) const;
    
protected:   
  //double x_   = 1;
  double alf_ = 1;
  double Rx_  = 0;
  double Ry_  = 0;
  double Rz_  = 0;

  // These are derived quantities we cash
  double norm_       = 1;
  double Rmax_       = 1;
  double prefactor_  = 1;
  double dprefactor_ = 0;

  int Nimage_x_ = 0;
  int Nimage_y_ = 0;
  int Nimage_z_ = 0;
  
  
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

	int ix = getIX(g.Rx());
	int iy = getIY(g.Ry());
	int iz = getIZ(g.Rz());

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
	while(rx-g.Rx() > L_[0]/2) rx -= L_[0]; while(rx-g.Rx() < -L_[0]/2) rx += L_[0];
	while(ry-g.Ry() > L_[1]/2) ry -= L_[1]; while(ry-g.Ry() < -L_[1]/2) ry += L_[1];
	while(rz-g.Rz() > L_[2]/2) rz -= L_[2]; while(rz-g.Rz() < -L_[2]/2) rz += L_[2];
	
	// sum over all contributing images
	for(int imx = -g.get_Nimage_x(); imx <= g.get_Nimage_x(); imx++)
	  for(int imy = -g.get_Nimage_y(); imy <= g.get_Nimage_y(); imy++)
	    for(int imz = -g.get_Nimage_z(); imz <= g.get_Nimage_z(); imz++)
	      g.get_measures(rx+imx*L_[0],ry+imy*L_[1],rz+imz*L_[2],hsd,fm);
      }
  }

 protected:
  vector<Gaussian> gaussians_;
};








#endif // __LUTSKO__GAUSSIAN_DENSITY__
