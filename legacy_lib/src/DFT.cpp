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

#include "poly34.h"

#include "Poly.h"

#ifdef USE_OMP
#include <omp.h>
#endif

#include "DFT.h"

double DFT::Mu(const vector<double> &x, int species) const
{
  Summation mu;

  mu += log(x[species]);

  if(fmt_)
    mu += fmt_->BulkMuex(x, allSpecies_, species);

  for(auto &interaction: Interactions_)
    mu += interaction->Mu(x,species);
  
  return mu.sum();
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

  double V = lattice().getVolume();

  for(auto &y: x)
    F += y*log(y)-y;

  double Fhs = 0.0;
  if(fmt_)
    Fhs += fmt_->BulkFex(x, allSpecies_);

  F += Fhs;

  double Fmf = 0.0;
  for(auto &interaction: Interactions_)
    Fmf += interaction->Fhelmholtz(x);

  F += Fmf;

  return F;  
}  

double DFT::XLiq_From_Mu(double mu, double high_density) const
{
  if(allSpecies_.size() > 1) throw std::runtime_error("Xliq_From_Mu only implemented for single component systems");

  vector<double> x(1);

  // Find the liquid that has the correct chemical potential,
  // Start at the high density and go dozn until the chemical potential is bracketed.
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
  } while(fabs(fa-fb) > 1e-12 + 1e-12*fabs(fa+fb));
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
//    dbetaP/dx (1+e*(4+e*(4+e*(-4+e))))/(e1*e1*e1*e1) + a_*x 
// so we wish to solve
// 1+4e+4e^2-4e^3+e4 + a_ x(1-4e+6e^2-4e^3+e^4)
void DFT::spinodal(double &xs1, double &xs2) const
{
  if(allSpecies_.size() != 1)   throw std::runtime_error("DFT::spinodal only implemented for a single component systems");
  if(Interactions_.size() != 1) throw std::runtime_error("DFT::spinodal only implemented for a single attractive interaction");

  //  double a[6] = { 1, 4 + ae, 4-4*ae , -4+6*ae, 1-4*ae, ae };
  //  int nroots = SolveP5(x, a[4]/a[5], a[3]/a[5], a[2]/a[5], a[1]/a[5], a[0]/a[5]);

  xs1 = xs2 = -1;
  
  double d = allSpecies_[0]->getHSD();
  double fac = 6/(M_PI*d*d*d);
  double ae = Interactions_[0]->getVDWParameter()*fac;

  vector<double> a;
  a.push_back(1);
  a.push_back( 4 + ae);
  a.push_back( 4 - 4*ae);
  a.push_back(-4 + 6*ae);
  a.push_back( 1 - 4*ae);
  a.push_back(       ae);
    
  double roots[2*(a.size()-1)];
  gsl_poly_complex_workspace * w= gsl_poly_complex_workspace_alloc(a.size());
  gsl_poly_complex_solve(a.data(), a.size(), w, roots);
  gsl_poly_complex_workspace_free(w);

  for(int i=0;i<a.size()-1;i++)
    if(fabs(gsl_poly_eval(a.data(),a.size(),roots[2*i])) < 1e-8) // only check real roots
      if(roots[2*i] > 0 && roots[2*i] < (M_PI/3)*M_SQRT1_2) // i.e. pi/(3sqrt(2)) =  close packing limit
	{
	  if(xs2 < 0) xs2 = roots[2*i];
	  else {
	    if(xs1 < 0) xs1 = min(roots[2*i],xs2);
	    else xs1 = min(xs1, roots[2*i]);
	    xs2 = max(xs2, roots[2*i]);
	  }
	}

  if(xs1 < 0 || xs2 < 0) throw std::runtime_error("Determination of spinodal failed 2");

  // convert from e to x
  xs1 *= 6/(M_PI*d*d*d);
  xs2 *= 6/(M_PI*d*d*d);
  return;
}

// Solve  P = x1*(1+e+e*e-e*e*e)*pow(1.0-e,-3) + 0.5*a*x1*x1 + 2b*x1*x1*x1
// or     (P*pi*d^3/6)(1.0-e)^3 = e*(1+e+e*e-e*e*e) + 0.5*ae*e*e*(1-e)^3 + 2*be*e*e*e*(1-e)^3
double DFT::XLiq_from_P(double P) const
{
  if(allSpecies_.size() != 1)   throw std::runtime_error("DFT::spinodal only implemented for a single component systems");
  if(Interactions_.size() != 1) throw std::runtime_error("DFT::spinodal only implemented for a single attractive interaction");

  double d = allSpecies_[0]->getHSD();
  double fac = 6/(M_PI*d*d*d);
  double ae = (Interactions_.size() == 0 ? 0.0 : Interactions_[0]->getVDWParameter()*fac);
  
  P *= M_PI*d*d*d/6;

  double x = -1;
  
  /*  
  double roots[5];
  int nroots = SolveP5(roots, -3+2/ae,(3*ae-2*P-2)/ae, (-ae+6*P-2)/ae, (-6*P-2)/ae, 2*P/ae);
  
  for(int i=0;i<nroots;i++)
    if(roots[i] > 0 && roots[i] < (M_PI/3)*M_SQRT1_2) // i.e. pi/(3sqrt(2)) =  close packing limit
      x = max(x, roots[i]);
  */
  /*
  const int order = 6;
  double a[] = { -P, 3*P+1, -3*P+1+0.5*ae, P+1-1.5*ae+2*be, -1+1.5*ae-6*be, -0.5*ae+6*be, -2*be};  
  double roots[order*2];
  gsl_poly_complex_workspace * w= gsl_poly_complex_workspace_alloc(order+1);
  gsl_poly_complex_solve(a, order+1, w, roots);
  gsl_poly_complex_workspace_free (w);
  int nroots = order;

  for(int i=0;i<nroots;i++)
    if(fabs(gsl_poly_eval(a,1+order,roots[2*i])) < 1e-12) // only check real roots
  */

  vector<double> a;
  a.push_back(-P);
  a.push_back(  3*P+1);
  a.push_back( -3*P+1+0.5*ae);
  a.push_back(  P+1-1.5*ae);
  a.push_back(  -1+1.5*ae);
  a.push_back(   -0.5*ae);
    
  double roots[2*(a.size()-1)];
  gsl_poly_complex_workspace * w= gsl_poly_complex_workspace_alloc(a.size());
  gsl_poly_complex_solve(a.data(), a.size(), w, roots);
  gsl_poly_complex_workspace_free(w);

  for(int i=0;i<a.size()-1;i++)
    if(fabs(gsl_poly_eval(a.data(),a.size(),roots[2*i])) < 1e-12) // only check real roots      
      if(roots[2*i] > 0 && roots[2*i] < (M_PI/3)*M_SQRT1_2) // i.e. pi/(3sqrt(2)) =  close packing limit
        x = max(x, roots[2*i]);

  if(x < 0) throw std::runtime_error("DFT::Xliq_From_P failed");

  // convert from e to x
  x *= 6/(M_PI*d*d*d);
  return x;
}

void DFT::liq_vap_coex(double &xs1, double &xs2, double &x1, double &x2) const
{
  xs1 = xs2 = x1 = x2 = -1;

  spinodal(xs1,xs2);
  if(xs1 < 0 || xs2 < 0) throw std::runtime_error("No liq-vap coexistence found");

  if(xs1 > xs2) swap(xs1,xs2);
  
  x1 = xs1;  
  vector<double> v1(1); v1[0] = x1;
  double Mu1 = Mu(v1,0);

  vector<double> v2(1);
  v2[0] = max(XLiq_from_P(-Omega(v1)),xs2);
  double Mu2 = Mu(v2,0);

  if(Mu2 > Mu1) throw std::runtime_error("DFT::liq_vap_coex failed at point 1");

  while(Mu2 < Mu1){ v1[0] /= 2; Mu1 = Mu(v1,0);}

  // the value is between v1[0] and v1[0]*2
  double a = v1[0];
  double b = 2*a;
  while(fabs(b-a) > 1e-12 + 1e-12*fabs(a+b))
    {
      v1[0] = (a+b)/2;
      v2[0] = max(XLiq_from_P(-Omega(v1)),xs2);      
      Mu1 = Mu(v1,0);
      Mu2 = Mu(v2,0);
      if(Mu1 < Mu2) a = v1[0];
      else b = v1[0];
    }
  x1 = v1[0];
  x2 = v2[0];
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
  F_id_ = F;

  // Hard-sphere contribution
  if(fmt_)
    {    
      try{
	F_hs_ = fmt_->calculateFreeEnergyAndDerivatives(allSpecies_);
	F += F_hs_;
      } catch( Eta_Too_Large_Exception &e) {
	throw e;
      }
    } else {   // Need the following only if the fmt object is not called
    for(auto &species : allSpecies_)
      species->doFFT();
  }

  //< Mean field contribution to F and dF
  F_mf_ = 0;
  for(auto &interaction: DFT::Interactions_)    
    F_mf_ += interaction->getInteractionEnergyAndForces();
  F += F_mf_;

  // External field + chemical potential
  F_ext_ = 0;
  for(auto &species : allSpecies_)
    F_ext_ += species->externalField(true); // bCalcForces = true: obsolete?
  F += F_ext_;

  return F.sum();  
}

  /**
   *   @brief  Computes (d2F/dn_i dn_j) v_j: Currently implemented only for a SINGLE SPECIES
   *
   *  
   *   @return nothing
   */  
void DFT::second_derivative(vector<DFT_FFT> &v, vector<DFT_Vec> &d2F, bool noFID)
{
  double dV = allSpecies_[0]->getDensity().dV();
    
  // ideal gas contribution: v_i/n_i
  if(noFID == false)
    for(int s=0;s<allSpecies_.size();s++)
      for(unsigned i=0;i<v[s].cReal().size();i++)
	d2F[s].set(i, dV*v[s].cReal().get(i)/allSpecies_[s]->getDensity().getDensity(i));
  
  // Hard-sphere
  if(fmt_)
    {    
      try{
	fmt_->add_second_derivative(v,d2F, allSpecies_);
      } catch( Eta_Too_Large_Exception &e) {
	throw e;
      }
    }
  
  if(!fmt_)
    for(auto &species : allSpecies_)
      species->doFFT();  
  
  // Mean field
  for(auto &interaction: DFT::Interactions_)    
    interaction->add_second_derivative(v,d2F);  
  
}






namespace dft_util {
/**
 *   @brief  This will fit the value of thethe vDW coefficient to the supplied data for a single species. The input densities are supplied in dimensionless form, 
 *           so they imply some length scale  l. The output is the hsd/l and the vdw parmeter in the form  beta a/l^3. Normally, for a potential with length scale sigma and energy scale epsilon, the
 *           calculated vdw parameter is some number times beta epsilon sigma^3 so one could match the two either by adjusting epsilon or sigma or both. 
 *  
 *   @param  data is the array of pairs of densities and z_factor = P/(n kT) for which the pressure is being supplied
 *   @param hsd: this is the hard-sphere diameter.
 *   @returns the vDW parameter in the dimensionless form beta avdw/l^3. 
 */   
double fit_avdw_to_data(vector<pair<double, double> > &data, double hsd)
{
  double s1 = 0.0;
  double s2 = 0.0;

  for(auto &p: data)
    {
      double e = p.first*M_PI*hsd*hsd*hsd/6;
      double Z = (1+e*(1+e*(1-e)))/pow(1-e,3);
      s1 += (p.second - Z)*p.first;
      s2 += p.first*p.first;
    }
  return 2*s1/s2;
}

static double objective(vector<pair<double, double> > &data, double hsd)
{
  double sum = 0.0;
  double a   = fit_avdw_to_data(data,hsd);

  for(auto &p: data)
    {      
      double e = p.first*M_PI*hsd*hsd*hsd/6;
      double Z = (1+e*(1+e*(1-e)))/pow(1-e,3);
      Z += 0.5*a*p.first;
      Z -= p.second;

      sum += Z*p.first*(2+2*e-e*e)*pow(1-e,-4);
    }
  return sum; 
}

/**
 *   @brief  This will fit the value of the hsd and the vDW coefficient to the supplied data for a single species. The input densities are supplied in dimensionless form, 
 *           so they imply some length scale  l. The output is the hsd/l and the vdw parmeter in the form  beta a/l^3. Normally, for a potential with length scale sigma and energy scale epsilon, the
 *           calculated vdw parmaeter is some number times beta epsilon sigma^3 so one could match the two either by adjusting epsilon or sigma or both. 
 *  
 *   @param  data is the array of pairs of densities and z_factor = P/(n kT) for which the pressure is being supplied
 *   @param hsd: this is an in/out parameter. The supplied value is used as an initial guess in the numerical search. Afterwards, it contains the determined value.
 *   @param aVDW this is an output parameter. It is actually  beta aVDW/l^3 . 
 *   @param tol is the tolerence of the fit. 
 */   

double fit_to_data(std::vector<std::pair<double, double> > &data, double &avdw, double tol)
{
  double a = 0;
  double b = 0;
  double fa = 0;
  double fb = 0;

  for(double d = 0.5; d < 2.0; d += 0.1)
    {
      a = b;
      b = d;
      fa = fb;
      fb = objective(data, d);
      if(fa*fb < 0) break;
    }

  if(fa*fb > 0) throw std::runtime_error("No crossing found in fit_to_data");

  while(fabs(a-b) > tol*(a+b))
    {
      double c = (a+b)/2;
      double fc = objective(data,c);
      if(fa*fc > 0) {a = c; fa = fc;}
      else { b = c; fb = fc;}
    }

  avdw = fit_avdw_to_data(data,(a+b)/2);  
  return (a+b)/2;
}
}  
