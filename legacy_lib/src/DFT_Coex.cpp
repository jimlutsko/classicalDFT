#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

using namespace std;

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_gamma.h>

#ifdef USE_OMP
#include <omp.h>
#endif

#include "DFT.h"

// Spinodal is defined by the stationary points of the pressure. Here, I will do
// the very simplest search for a stationary point.

static vector<double> xv(1);

static double pressure(double x, const DFT* dft)
{
  xv[0] = x;
  return -dft->Omega(xv);
}

static double chempot(double x, const DFT* dft)
{
  xv[0] = x;
  return dft->Mu(xv,0);
}

void DFT::findSpinodal(double xmax, double dx, double &xs1, double &xs2, double tol) const
{
  // 1. Bracket the maximum
  double x = 2*dx;
  double p0 = pressure(dx,this);
  double p  = pressure(2*dx,this);
  double dp = p-p0;

  if(dp < 0) throw std::runtime_error("0. Could not get started in DFT::findSpinodal");

  while(dp > 0 && x < xmax-dx)
    {
      dp = -p;
      x   += dx;
      p = pressure(x,this);
      dp += p;
    }
  if(x >= xmax-dx) throw std::runtime_error("1. Xmax exceeded in DFT::findSpinodal");
  
  // 2. Refine
  // first point is bracket between x and x-2*dx
  double a = x-2*dx;
  double b = x;
  double fa = pressure(a,this);
  double fb = pressure(b,this);

  double r = (3-sqrt(5))/2;
  double u = a+r*(b-a);
  double v = b-r*(b-a);

  double fu = pressure(u,this);
  double fv = pressure(v,this);


  // sanity check
  double fsanity = pressure(a+dx,this);
  if(fsanity < fa || fsanity < fb) throw std::runtime_error("Failed sanity check 1 in DFT::findSpinodal");
  
  do {
    if(fu > fv) { b = v; fb = fv; v = u; fv = fu; u = a+r*(b-a); fu = pressure(u,this);}
    else { a = u; fa = fu; u = v; fu = fv; v = b-r*(b-a);  fv = pressure(v,this);}
  } while(b-a > 1e-8);
  xs1 = (a+b)/2;

  // 3. Bracket the minimum.
  x = xs1;
  p = pressure(x,this);
  
  do {
    dp = -p;
    x   += dx;
    p = pressure(x,this);
    dp += p;
  } while(dp < 0 && x < xmax-dx);
  if(x >= xmax-dx) throw std::runtime_error("2. Xmax exceeded in DFT::findSpinodal");

  //4. Refine
  // first point is bracket between x and x-2*dx
  a = x-2*dx;
  b = x;
  fa = pressure(a,this);
  fb = pressure(b,this);

  // sanity check
  fsanity = pressure(a+dx,this);
  if(fsanity > fa || fsanity > fb) throw std::runtime_error("Failed sanity check 2 in DFT::findSpinodal");

  r = (3-sqrt(5))/2;
  u = a+r*(b-a);
  v = b-r*(b-a);

  fu = pressure(u,this);
  fv = pressure(v,this);
  
  do {
    if(fu < fv) { b = v; fb = fv; v = u; fv = fu; u = a+r*(b-a); fu = pressure(u,this);}
    else { a = u; fa = fu; u = v; fu = fv; v = b-r*(b-a);  fv = pressure(v,this);}
  } while(b-a > tol);
  xs2 = (a+b)/2;
  
  return;
}

double  DFT::find_density_from_mu(double mu, double xmin, double xmax, double tol) const
{
  double mu1 = chempot(xmin,this);
  double mu2 = chempot(xmax,this);

  if(mu1 > mu2) { double tmp = xmin; xmin = xmax; xmax = tmp; tmp = mu1; mu1 = mu2; mu2 = tmp;}
  
  if(mu2 < mu || mu1 > mu)
    {
      cout << mu1 << " " << mu << " " << mu2 << endl;
      throw std::runtime_error("DFT::find_density_from_mu failed");
    }
  
  do {
    double x = (xmin+xmax)/2;
    double mux = chempot(x,this);
    if(mux > mu) {xmax = x;  mu2 = mux;}
    else {xmin = x; mu1 = mux;}
  } while(mu2-mu1 > tol);

  return (xmin+xmax)/2;
}

void DFT::findCoex(double xmax, double dx, double &x1, double &x2, double tol) const
{
  // 0. Get limits of vap and liq regions: (0,xs1) and (xs2,xmax)
  double xs1,xs2;
  findSpinodal(xmax,dx,xs1,xs2,tol);

  // 1. Bracket  
  double xvap = xs1;
  double xliq = find_density_from_mu(chempot(xvap,this), xs2,xmax-dx, tol);
  double dp1 = pressure(xvap,this) - pressure(xliq,this);
  double dp2;
  do {
    dp2 = dp1;
    xvap -= dx;
    xliq = find_density_from_mu(chempot(xvap,this), xs2,xmax-dx,tol);      
    dp1 =  pressure(xvap,this) - pressure(xliq,this);
  } while(xvap > dx && ((dp1<0) == (dp2<0)));

  if(xvap <= dx) throw std::runtime_error("DFT::findCoex failed 1");

  
  // 2. Bracket is (xvap,xvap+dx) : refine

  double y1 = xvap;
  double y2 = xvap+dx;  
  do {
    double y = (y1+y2)/2;
    xliq = find_density_from_mu(chempot(y,this), xs2,xmax-dx,tol);      
    double dp =  pressure(y,this) - pressure(xliq,this);    
    if((dp<0) == (dp1<0)) { y1 = y; dp1 = dp;}
    else { y2 = y; dp2 = dp;}
  } while(fabs(y2-y1) > tol);

  x1 = (y1+y2)/2;
  x2 = xliq;
}
