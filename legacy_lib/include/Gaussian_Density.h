#ifndef __LUTSKO__GAUSSIAN_DENSITY__
#define __LUTSKO__GAUSSIAN_DENSITY__

#include "Density.h"
#include "Fundamental_Measures.h"
#include "Gaussian.h"


class GaussianDensity : public Density
{
 public:
  GaussianDensity(double dx, double L[], int num_gaussians = 0, double hsd = 1, double Acut = 4, bool use_discrete = true)
    : Density(dx, L), gaussians_(num_gaussians, Gaussian(Acut,hsd)), use_discrete_(use_discrete)
  {
    Density_.zeros();
    if(use_discrete)
      discrete_gaussian_field.initialize(Nx_, Ny_, Nz_);
  }

  ~GaussianDensity(){}

  bool use_discrete() const {return use_discrete_;}

  void fill_discrete_gaussian_field();  
  double get_discrete_gauss_at(double pos) { return discrete_gaussian_field.Real().get(pos);}

  void set_hsd(double hsd){ for(auto& g: gaussians_) g.set_hsd(hsd);}
  
  
  virtual size_t number_of_parameters() const { return (use_discrete_ ? Ntot_ : 0) + 5*gaussians_.size();}

  size_t number_of_gaussians() const { return gaussians_.size();}
  
  void set_gaussian(int ig, double y, double alf, double Rx, double Ry, double Rz)
  {
    gaussians_[ig].set_parameters(y,alf,Rx,Ry,Rz);
  }

  void get_gaussian(int ig, double &y, double &alf, double &Rx, double &Ry, double &Rz)
  {
    gaussians_[ig].get_parameters(y,alf,Rx,Ry,Rz);
  }

  const Gaussian& get_gaussian(int ig) const { return gaussians_[ig];}
  
  
  virtual double getNumberAtoms() const
  {
    double N = 0;
    for(auto &g: gaussians_) N += g.prefactor();
    return N;
  }

  void get_dN(int ig, double &dN_dx, double &dN_dalf) const;
    
  virtual double getRsquared() const
  {
    double r2 = 0;
    for(auto &g: gaussians_) r2 += g.prefactor()/g.alf();
    return r2/getNumberAtoms();
  }
  
  void get_measures(long pos, double hsd, FundamentalMeasures &fm) const
  {
    long ix,iy,iz;
    cartesian(pos,ix,iy,iz);
    get_measures(getX(ix), getY(iy), getZ(iz), hsd, fm);
  }
  
  void get_measures(double rx, double ry, double rz, double hsd, FundamentalMeasures &fm) const;
  void get_dmeasures_for_gaussian(int igaussian, double rx, double ry, double rz, double hsd, FundamentalMeasures dfm[5]) const;
  double FMF(double w0, double r0, vector<double> &x, vector<double> &w, DFT_Vec &dF) const;  
  
 protected:
  vector<Gaussian> gaussians_;
  bool use_discrete_ = true;
  DFT_FFT discrete_gaussian_field;
};








#endif // __LUTSKO__GAUSSIAN_DENSITY__
