#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <utility>
#include <algorithm>

using namespace std;

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_poly.h>

#include "VDW1.h"
#include "poly34.h"

// Solve P= x*(1+e*(1+e*(1-e)))/(e1*e1*e1) + a_*x*x for xmin<x0<x<xmax
/*
double VDW1::findDensityFromPressure(double P, double xmin, double xmax) const
{
  P *= M_PI*d_*d_*d_/6;
  xmin *= M_PI*d_*d_*d_/6;
  xmax *= M_PI*d_*d_*d_/6;

  double ae = a_*6/(M_PI*d_*d_*d_);
  double a[6] = { -P, 1 + 3*P, 1 - 3*P +ae, 1 +P - 3*ae, -1 + 3*ae, -ae};
  double z[10];

  gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc(6);
  int ret =   gsl_poly_complex_solve(a, 6, w, z);
  gsl_poly_complex_workspace_free(w);

  if(ret != GSL_SUCCESS) {cout << " P = " << P << endl; throw std::runtime_error("findDensityFromPressure poly solver failed");}

  double x = -1;
  double xi = 1e30;

  for (int i = 0; i < 5; i++)
    {
      double re = z[2*i];
      double im = abs(z[2*i+1]);

      if(re > xmin && re < xmax)
	if(im < xi) 
	  {x = re; xi = im;}
    }
  if(x < 0 || xi > 1e-6) throw std::runtime_error("VDW1::findDensityFromPressure failed");

  // convert from e to x
  x *= 6/(M_PI*d_*d_*d_);

  return x;
}
*/

 // Solve P= x*(1+e*(1+e*(1-e)))/(e1*e1*e1) + a_*x*x for xmin<x0<x<xmax
double VDW1::findDensityFromPressure(double P, double xmin, double xmax) const
{
  P *= M_PI*d_*d_*d_/6;
  xmin *= M_PI*d_*d_*d_/6;
  xmax *= M_PI*d_*d_*d_/6;

  double ae = a_*6/(M_PI*d_*d_*d_);
  double a[6] = { -P, 1 + 3*P, 1 - 3*P +ae, 1 +P - 3*ae, -1 + 3*ae, -ae};

  double x[5];
  double d = SolveP5(x, (1-3*ae)/ae, (3*ae-1-P)/ae, (3*P-1-ae)/ae, -(1+3*P)/ae, P/ae);

  double y = -1;
  for(int i=0;i<d;i++)
    if((1+x[i]*(4+x[i]*(4+x[i]*(-4+x[i]))))*pow(1-x[i],-4)+2*ae*x[i] >= 0) // non-negative compressibility
      if(x[i] > xmin && x[i] < xmax)
	y = x[i];
  if(y < 0) throw std::runtime_error("VDW1::findDensityFromPressure failed");

  // convert from e to x
  y *= 6/(M_PI*d_*d_*d_);

  return y;
}




// For coexistence we need to solve P1=P2 and mu1=mu2
// where, for the VDW model,
// 
// P = x*(1+e*(1+e*(1-e)))/(e1*e1*e1) + a_*x*x
// mu = log(x) + (e*(8-9*e+3*e*e)/(e1*e1*e1)) + 2*a_*x
//
// The strategy is to first find the spinodal (it is easy)
// Then, to use the fact that the coexisting vapor pressure is between zero and the spinodal
// So for a given vapor x, we compute P and then solve for the corresponding liq x (if it exists)
// and then do bisection on mu

int VDW1::findCoexistence(double &x1, double &x2) const
{
  double x_liq_ceil = 2.0;
  double x_vap_floor = 1e-30;
  
  // Spinodal points
  double xs1,xs2;
  spinodal(xs1,xs2);
  if(xs1 > xs2) { double t = xs1; xs1 = xs2; xs2 = t;}
  
  double Ps1 = pressure(xs1);
  double Ps2 = pressure(xs2);

  // The lowest vapor density corresponds to max(Ps2,0.0)
  double x_vap_min = -1;

  if(Ps2 < 0) {Ps2 = pressure(x_vap_floor); xs2 = findDensityFromPressure(Ps2,xs2,x_liq_ceil); x_vap_min = x_vap_floor;}
  else x_vap_min = findDensityFromPressure(Ps2,0.0,xs1);

  // The max vapor density is either at the spinodal or the largest value below the spinodal for which
  // a liq with the same pressure exists.
  double x_vap_max = xs1;
  double ya        = x_vap_min;
  double y         = x_vap_max;
  
  do{
    try {
      findDensityFromPressure(pressure(y),xs2,x_liq_ceil);
      ya = y;
    } catch(...) {
      x_vap_max = y;
    }
    y = (ya+x_vap_max)/2;
  } while(fabs(y-x_vap_max) > 1e-10);
  x_vap_max = y;

  // Now, we expect that the two densities bracket the coexistence point: is this true?

  double y_liq_a = xs2;
  try {
    y_liq_a = findDensityFromPressure(pressure(x_vap_min),xs2,x_liq_ceil);
  } catch (...) {  }
    
  double y_liq_b = findDensityFromPressure(pressure(x_vap_max),xs2,x_liq_ceil);
  
  double dmu_a = chemPotential(x_vap_min) - chemPotential(y_liq_a);
  double dmu_b = chemPotential(x_vap_max) - chemPotential(y_liq_b);
  if(dmu_a*dmu_b > 0) 
      throw std::runtime_error("Cannot find low enough chem potential in VDW1::findCoexistence");

  // good - just do bisection
  double dmu;
  do {
    double x = (x_vap_min+x_vap_max)/2;
    double y = findDensityFromPressure(pressure(x),xs2,x_liq_ceil);

    dmu = chemPotential(x) - chemPotential(y);

    //    cout << setprecision(12) << x_vap_min << " " << x_vap_max << " " << x << " " << y << " " << dmu_a << " " << dmu_b << " " << dmu << endl;
    
    if(dmu*dmu_a > 0) { dmu_a = dmu; x_vap_min = x;}
    else {dmu_b = dmu; x_vap_max = x;}

  } while(fabs(dmu) > 1e-6);

  x1 = (x_vap_max+x_vap_min)/2;
  x2 = findDensityFromPressure(pressure(x1),xs2,x_liq_ceil);

  return 0;
}


// To find the spinodal, we just need the roots of dP/dx
// and here 
//     double e = M_PI*x*d_*d_*d_/6;
//    double e1 = 1-e;
//    dbetaP/dx (1+e*(4+e*(4+e*(-4+e))))/(e1*e1*e1*e1) + 2*a_*x;
// so we wish to solve
// 1+4e+4e^2-4e^3+e4 + 2a_ x(1-4e+6e^2-4e^3+e^4)
void VDW1::spinodal(double &xs1, double &xs2) const
{
  int i;
  // coefficients of P(x) = a0+a1 x+...
  double ae = a_*6/(M_PI*d_*d_*d_);
  double a[6] = { 1, 4 + 2*ae, 4-8*ae , -4+12*ae, 1-8*ae, 2*ae };
  double z[10];

  gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc(6);
  int ret =   gsl_poly_complex_solve(a, 6, w, z);
  gsl_poly_complex_workspace_free(w);

  if(ret != GSL_SUCCESS) throw std::runtime_error("Determination of spinodal failed 1");

  xs1 = xs2 = -1;

  for (i = 0; i < 5; i++)
    {
      double re = z[2*i];
      double im = z[2*i+1];
      if(re > 0 && fabs(im) < 1e-10)
	{
	  if(xs1 < 0 || (re < xs1)) {xs2 = xs1; xs1 = re;}
	  else if(xs2 < 0 || re < xs2) xs2 = re;
	}
    }
  if(xs1 < 0 || xs2 < 0) throw std::runtime_error("Determination of spinodal failed 2");

  // convert from e to x
  xs1 *= 6/(M_PI*d_*d_*d_);
  xs2 *= 6/(M_PI*d_*d_*d_);
  return;
}

double VDW1::findLiquidFromMu(double mu, double mu_coex, double xliq_coex) const
{
  // Find the liquid that has the correct chemical potential,
  // We know it is between densityLiqCoex and infinity
  double ax = xliq_coex;
  double bx = ax;
  double fa = 0.0;
  double fb = 0.0;
    
  double SuperSaturation = (mu - mu_coex)/fabs(mu_coex);

  int sign = (SuperSaturation > 0 ? -1 : 1);
    
  while(sign*fb > sign*SuperSaturation)
    {
      ax = bx;
      bx += -sign*0.01;
      fb = (chemPotential(bx)-mu_coex)/fabs(mu_coex);
    }
  double densityLiq;
  do {
    densityLiq = (ax+bx)/2;
    double fc = (chemPotential(densityLiq)-mu_coex)/fabs(mu_coex);
    if(sign*fc < sign*SuperSaturation){bx = densityLiq; fb = fc;}
    else { ax = densityLiq; fa = fc;}
  } while(fabs(fa-fb) > 1e-6*fabs(SuperSaturation));
  return densityLiq;
}


double VDW1::findLiquidFromMu(double mu, double high_density) const
{
  // Find the liquid that has the correct chemical potential,
  // We know it is between densityLiqCoex and infinity
  double ax = high_density;
  double bx = ax;
  double fa = chemPotential(ax);
  double fb = 0.0;
        
  do {
    bx = ax; fb = fa;
    ax -= 0.01;
    fa = chemPotential(ax);
  } while((fa-mu)*(fb-mu) > 0);

  double cx;
  do {
    cx = (ax+bx)/2;
    double fc = chemPotential(cx);

    if(fa > mu)
      {
	if(fc > mu) { fa = fc; ax = cx;}
	else {fb = fc; bx = cx;}
      } else {
      	if(fc > mu) { fb = fc; bx = cx;}
	else {fa = fc; ax = cx;}
    }
  } while(fabs(fa-fb) > 1e-6);
  return cx;
}

double VDW1::findVaporFromMu(double betamu, double maxDensity) const
{
  // first, find a lower bound:
  // this is a guess for ideal gas
  double x1 = exp(betamu)/2;
  while(chemPotential(x1) > betamu) x1 /= 2;

  // now we need an upper bound: is there a spinodal?
  double xs1 = -1;
  double xs2 = -1;
  double x2  = -1;
  try {
    spinodal(xs1,xs2);
    x2 = min(xs1,xs2);
  } catch(...) {
    // no problem - there is no spinodal.
    x2 = maxDensity;
  }

  if(chemPotential(x2) < betamu) return -1;

  // do bisection
  double f1 = chemPotential(x1);
  double f2 = chemPotential(x2);

  while(fabs(x1-x2) > 1e-8)
    {
      double x = (x1+x2)/2;
      double f = chemPotential(x);
      if(f > betamu) x2 = x;
      else x1 = x;
    }
  return (x1+x2)/2;  
}


