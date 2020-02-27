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

#include "poly34.h"

#ifdef USE_OMP
#include <omp.h>
#endif

#include "DFT.h"

double DFT::Mu(const vector<double> &x, int species) const
{
  double mu = log(x[species]);

  if(fmt_)
    mu += fmt_->BulkMuex(x, allSpecies_, species);

  for(auto &interaction: Interactions_)
    mu += interaction->Mu(x,species);
  return mu;
}

double DFT::Omega(const vector<double> &x) const
{
  double omega = Fhelmholtz(x);
  for(int i=0;i<allSpecies_.size();i++)
    {
      omega -= x[i]*Mu(x,i);
    }
  
  return omega;
}

double DFT::Fhelmholtz(const vector<double> &x) const
{
  double F = 0.0;

  for(auto &y: x)
    F += y*log(y)-y;

  if(fmt_)
    F += fmt_->BulkFex(x, allSpecies_);

  for(auto &interaction: Interactions_)
    F += interaction->Fhelmholtz(x);  
  return F;  
}  

double DFT::XLiq_From_Mu(double mu, double high_density) const
{
  if(allSpecies_.size() > 1) throw std::runtime_error("Xliq_From_Mu only implemented for single component systems");

  vector<double> x(1);

  // Find the liquid that has the correct chemical potential,
  // Start at teh high density and go dozn until the chemical potential is bracketed.
  double ax = high_density;
  double bx = ax;
  double fa = 0.0;
  double fb = 0.0;

  x[0] = ax; fa = Mu(x,0);
  
  do {
    bx = ax; fb = fa;
    ax -= 0.01;
    if(ax < 0) throw std::runtime_error("No liquid found");

    x[0] = ax; fa = Mu(x,0);
  } while((fa-mu)*(fb-mu) > 0);

  double cx;
  do {
    cx = (ax+bx)/2; x[0] = cx;
    double fc = Mu(x,0);

    if(fa > mu)
      {
	if(fc > mu) { fa = fc; ax = cx;}
	else {fb = fc; bx = cx;}
      } else {
      	if(fc > mu) { fb = fc; bx = cx;}
	else {fa = fc; ax = cx;}
    }
  } while(fabs(fa-fb) > 1e-8 + 1e-6*fabs(fa+fb));
  return cx;
}

double DFT::XVap_From_Mu(double mu, double maxDensity) const
{
  if(allSpecies_.size() > 1) throw std::runtime_error("XVap_From_Mu only implemented for single component systems");

  vector<double> x(1);
  
  // first, find a lower bound:
  // this is a guess for ideal gas
  x[0] = exp(mu)/2;
  while(Mu(x,0) > mu) x[0] /= 2;
  double x1 = x[0];
  
  // now we need an upper bound: is there a spinodal?
  double xs1 = -1;
  double xs2 = -1;
  double x2  = -1;
  try {
    spinodal(xs1,xs2);
    cout << "Spinodal: " << x1 << " " << x2 << endl;
    x2 = min(xs1,xs2);
  } catch(...) {
    // no problem - there is no spinodal.
    x2 = maxDensity;
  }

  x[0] = x2;
  if(Mu(x,0) < mu) return -1;

  // do bisection
  x[0] = x1; double f1 = Mu(x,0);
  x[0] = x2; double f2 = Mu(x,0);
  double f = f1;
  
  while(fabs(f-mu) > 1e-8*fabs(f+mu))
    {
      x[0] = (x1+x2)/2;
      f = Mu(x,0);
      if(f > mu) x2 = x[0];
      else x1 = x[0];
    }
  return x[0];
}

// To find the spinodal, we just need the roots of dP/dx
// and here 
//     double e = M_PI*x*d_*d_*d_/6;
//    double e1 = 1-e;
//    dbetaP/dx (1+e*(4+e*(4+e*(-4+e))))/(e1*e1*e1*e1) + a_*x;
// so we wish to solve
// 1+4e+4e^2-4e^3+e4 + 2a_ x(1-4e+6e^2-4e^3+e^4)
void DFT::spinodal(double &xs1, double &xs2) const
{
  if(allSpecies_.size() != 1)   throw std::runtime_error("DFT::spinodal only implemented for a single component systems");
  if(Interactions_.size() != 1) throw std::runtime_error("DFT::spinodal only implemented for a single attractive interaction");


  double d = allSpecies_[0]->getHSD();
  double ae = Interactions_[0]->getVDWParameter()*6/(M_PI*d*d*d);
  double a[6] = { 1, 4 + ae, 4-4*ae , -4+6*ae, 1-4*ae, ae };
  double x[5];
  int nroots = SolveP5(x, a[4]/a[5], a[3]/a[5], a[2]/a[5], a[1]/a[5], a[0]/a[5]);

  xs1 = xs2 = -1;
  for(int i=0;i<nroots;i++)
    if(x[i] > 0 && x[i] < 0.74) // i.e. < close packing limit
      {
	if(xs2 < 0) xs2 = x[i];
	else {
	  if(xs1 < 0) xs1 = min(x[i],xs2);
	  else xs1 = min(xs1, x[i]);
	  xs2 = max(xs2, x[i]);
	}
      }
  /*  
  vector<double> x(1);

  // coefficients of P(x) = a0+a1 x+...
  double d = allSpecies_[0]->getHSD();
  double ae = Interactions_[0]->getVDWParameter()*6/(M_PI*d*d*d);
  double a[6] = { 1, 4 + ae, 4-4*ae , -4+6*ae, 1-4*ae, ae };
  double z[10];

  gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc(6);
  int ret = gsl_poly_complex_solve(a, 6, w, z);
  gsl_poly_complex_workspace_free(w);

  if(ret != GSL_SUCCESS) throw std::runtime_error("Determination of spinodal failed 1");

  xs1 = xs2 = -1;

  for(int i = 0; i < 5; i++)
    {
      double re = z[2*i];
      double im = z[2*i+1];
      cout << re << " " << im << endl;
      if(re > 0 && fabs(im) < 1e-10)
	{
	  if(xs1 < 0 || (re < xs1)) {xs2 = xs1; xs1 = re;}
	  else if(xs2 < 0 || re < xs2) xs2 = re;
	}
    }
  */
  if(xs1 < 0 || xs2 < 0) throw std::runtime_error("Determination of spinodal failed 2");

  // convert from e to x
  xs1 *= 6/(M_PI*d*d*d);
  xs2 *= 6/(M_PI*d*d*d);
  return;
}


double DFT::calculateFreeEnergyAndDerivatives(bool onlyFex)
{
  for(auto &species : allSpecies_)  
    species->zeroForce();
  
  for(auto &s: allSpecies_) s->beginForceCalculation();

  double F = calculateFreeEnergyAndDerivatives_internal_(onlyFex);

  for(auto &s: allSpecies_)
    s->endForceCalculation();

  return F;
}

#pragma omp declare reduction(SummationPlus: Summation: omp_out += omp_in) 

double DFT::calculateFreeEnergyAndDerivatives_internal_(bool onlyFex)
{
  Summation F;
  double f = 0;
  // Ideal gas contribution  
  if(!onlyFex) 
    for(auto &species : allSpecies_)
      {
	const Density& density = species->getDensity();
	double dV = density.dV();
	long Ntot = density.Ntot();
	long pos;
#pragma omp parallel for				\
  shared(species, dV)				\
  private(pos)						\
  schedule(auto)				\
  reduction(SummationPlus:F)
	for(pos=0;pos<Ntot;pos++)
	  {
	    double d0 = density.getDensity(pos);
	    F += (d0*log(d0)-d0)*dV;
	    species->addToForce(pos,log(d0)*dV);
	  }
      }
  cout << "Fid = " << F << endl;
  // Hard-sphere contribution
  if(fmt_)
    {    
      try{
	f = fmt_->calculateFreeEnergyAndDerivatives(allSpecies_);
	F += f;
      } catch( Eta_Too_Large_Exception &e) {
	throw e;
      }
    }
  cout << "Fhs = " << f << endl;  
  //< Mean field contribution to F and dF
  // Need the following only if the fmt object is not called
  if(!fmt_)
    for(auto &species : allSpecies_)
      species->doFFT();

  f = 0;
  for(auto &interaction: DFT::Interactions_)    
    f += interaction->getInteractionEnergyAndForces();
  F += f;

  cout << "Fmf = " << f << endl;
  // External field + chemical potential
  f = 0;
  for(auto &species : allSpecies_)
    f += species->externalField(true); // bCalcForces = true: obsolete?
  F += f;
  cout << "Fxt = " << f << endl;  
  return F.sum();  
}
