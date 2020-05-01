#ifndef __LUTSKO_VDW__
#define __LUTSKO_VDW__

class VDW1
{
 public:
 VDW1(double d, double a) : d_(d), a_(a){}
  ~VDW1(){}

  // beta * f = beta *F/V
  double helmholtzPerUnitVolume(double x) const
  {
    double e = M_PI*x*d_*d_*d_/6;
    return x*log(x) - x + x*(e*(4.0-3.0*e)/(1-2*e+e*e)) + a_*x*x;
  }

  // beta * Mu = beta * df/dx
  double chemPotential(double x) const
  {
    double e = M_PI*x*d_*d_*d_/6;
    double e1 = 1-e;
    return log(x) + (e*(8-9*e+3*e*e)/(e1*e1*e1)) + 2*a_*x;
  }
  
  // beta * P  = -beta*(f - Mu x)
  double pressure(double x) const
  {
    double e = M_PI*x*d_*d_*d_/6;
    double e1 = 1-e;
    return x*(1+e*(1+e*(1-e)))/(e1*e1*e1) + a_*x*x;
  }

  // beta * dMu/dx
  double dChemPotential(double x) const
  {
    double e = M_PI*x*d_*d_*d_/6;
    double e1 = 1-e;
    return (1.0/x) + (M_PI*d_*d_*d_/6)*((8-2*e)/(e1*e1*e1*e1)) + 2*a_;
  }

  // beta * dPressure/dx
  double dPressure(double x) const
  {
    double e = M_PI*x*d_*d_*d_/6;
    double e1 = 1-e;
    return (1+e*(4+e*(4+e*(-4+e))))/(e1*e1*e1*e1) + 2*a_*x;
  }

  double findDensityFromPressure(double P, double xmin, double xmax) const;
  
  int findCoexistence(double &x1, double &x2) const;

  double findLiquidFromMu(double mu, double mu_coex, double xliq_coex) const ;
  double findLiquidFromMu(double mu, double high_density) const ;
  double findVaporFromMu(double betamu, double maxDensity) const;
  
  void spinodal(double &x1, double &x2) const;
  
  void set_VDW_Parameter(double a) { a_ = a;}
  void set_HardSphere_Diameter(double d) { d_ = d;}

  double get_VDW_Parameter() const { return a_;}
  
 private:
  double d_;
  double a_;
};

#endif // __LUTSKO_VDW__
