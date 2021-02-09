#ifndef __LUTSKO__GAUSSIAN__
#define __LUTSKO__GAUSSIAN__

#include "Fundamental_Measures.h"

class FMT_Gaussian_Species;


class Gaussian
{
 public:
  // Parameters are always set by calling set_parameters
  Gaussian(double Acut = 100, double hsd = 1): Acut_(Acut), hsd_(hsd)  {}
  ~Gaussian() {}

  // Setting the parameters
  // Note that the parameter y is an alias controlling the vacancy concentration
  void set_parameters(double x,  double alf,  double Rx,  double Ry,  double Rz);
  void get_parameters(double &x, double &alf, double &Rx, double &Ry, double &Rz) const; 

  void set_hsd(double hsd) { hsd_ = hsd; update_cache();}
  
  double get_prefac_alias(double prefac, double alf) 
  {
    double R    = hsd_/2;    
    double z    = sqrt(alf)*min(R,Acut_);
    double xmax = 1.0/(erf(z)-M_2_SQRTPI*z*exp(-z*z));

    Gaussian g;
    return ginv(prefac/xmax);
  }

  void print()
  {
    cout << setw(20) << prefactor_ << setw(20) << alf_ << setw(20) << Rx_ << setw(20) << Ry_ << setw(20) << Rz_ << setw(20) << Acut_ << setw(20) << hsd_ << endl;
  }
  
  // maximum range of this gaussian
  double Rmax() const { return Rmax_;}

  // Density
  double density(double rx, double ry, double rz) const;
  double density(double r[3]) const { return density(r[0], r[1], r[2]);}

  void zeroForces() { dx_ = dalf_ = dRx_ = dRy_= dRz_ = 0.0;}

  double get_N() const { return N_;};
  double get_ideal_f(double &dfdx, double &dfda) const{ dfdx = df_id_dx_; dfda = df_id_da_; return f_id_;}
  
  // accessors
  double prefactor()       const { return prefactor_;}
  double dprefactor_dx()   const { return dprefactor_dx_;}
  double dprefactor_dalf() const { return dprefactor_dalf_;}
  double alf()             const { return alf_;}
  double A()               const { return Acut_;}
  
  int get_Nimage_x() const { return Nimage_x_;}
  int get_Nimage_y() const { return Nimage_y_;}
  int get_Nimage_z() const { return Nimage_z_;}

  double Rx() const { return Rx_;}
  double Ry() const { return Ry_;}
  double Rz() const { return Rz_;}
  
  // FMT measures
  void get_measures(double rx, double ry, double rz, FundamentalMeasures &fm) const;
  void get_dmeasures_dAlf(double rx, double ry, double rz, FundamentalMeasures &dfm) const;
  void get_dmeasures_dR(double rx, double ry, double rz, FundamentalMeasures dfm[3]) const;
  void get_dmeasures_dX(double rx, double ry, double rz, FundamentalMeasures &dfm) const;

  friend FMT_Gaussian_Species;
  
protected:
  // these are the alias function, its inverse and its derivative
  // Maybe better if these were members that could be over-ridden
  static double    g(double x) { return 1.0-exp(-x*x);}
  static double ginv(double z) { return sqrt(fabs(log(1.0-z)));}
  static double   g1(double x) { return 2*x*exp(-x*x);}

  void get_f1f2f3f4(double y, double r2, double &f1, double &f2, double &f3, double &f4) const;
  void get_measures_internal(double rx, double ry, double rz, double coeff, FundamentalMeasures &fm) const;
  
  void update_cache();
  
protected:   
  double x_   = 1;
  double alf_ = 1;
  double Rx_  = 0;
  double Ry_  = 0;
  double Rz_  = 0;

  double hsd_ = 1;

  // forces
  double dx_   = 0.0;
  double dalf_ = 0.0;
  double dRx_  = 0.0;
  double dRy_  = 0.0;
  double dRz_  = 0.0;
  
  // These are derived quantities we cache
  double norm_            = 1;
  double Rmax_            = 1;
  double prefactor_       = 1;
  double dprefactor_dx_   = 0;
  double dprefactor_dalf_ = 0;
  double sq_alf_          = 0;
  double exp_A_           = 0;
  double erf_A_           = 0;

  int Nimage_x_ = 0;
  int Nimage_y_ = 0;
  int Nimage_z_ = 0;
  
  double Acut_ = 10.0;

  double N_        = 0.0;
  double f_id_     = 0.0;
  double df_id_dx_ = 0.0;
  double df_id_da_ = 0.0;
  
};
   
#endif // __LUTSKO__GAUSSIAN__
