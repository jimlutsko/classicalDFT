 /* This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
 * To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
 *
 * Author: James F. Lutsko
 * www.lutsko.com
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

#include <gsl/gsl_integration.h>

using namespace std;

#ifdef USE_OMP
#include <omp.h>
#endif

#include "Potential1.h"
#include "myColor.h"

BOOST_SERIALIZATION_ASSUME_ABSTRACT(Potential1)

BOOST_CLASS_EXPORT(Potential1)
BOOST_CLASS_EXPORT(LJ)
BOOST_CLASS_EXPORT(tWF)
BOOST_CLASS_EXPORT(WHDF)

Potential1::Potential1(double sigma, double eps, double rcut)
  : sigma_(sigma), eps_(eps),rcut_(rcut), shift_(0), rmin_(0), Vmin_(0) {}

Potential1::Potential1() {}

double Potential1::getRcut() const { return rcut_;} ///< Returns the cutoff
double Potential1::getR0()   const { return r0_;} ///< Returns the position at which potential is zero
  
void Potential1::set_WCA_limit(double r) { r_att_min_ = r;} ///< Distance from origin that attractive continuate extends
void Potential1::setBH() { bhFlag_ = true; r_att_min_ = getR0();} ///< Use the BH split of the potential
  
double Potential1::V(double r)   const { return vr(r)-shift_;}    ///< The cut and shifted potential at point r
double Potential1::V2(double r2) const { return vr2(r2)-shift_;}  ///< The cut and shifted potential at point r calculated from r^2

double Potential1::Watt(double r) const {return Watt2(r*r);} ///< The attractive tail 
double Potential1::Watt2(double r2) const ///< Attractive part calcuated from r2
{
  double ret = 0.0;
  if(r2 < rcut_*rcut_)  // zero outside of cutoff      
    if(bhFlag_) ret =  (r2 < getR0()*getR0() ? 0.0 : V2(r2)); // BH case 
    else if(r2 < r_att_min_*r_att_min_) ret = 0.0; // Our "generalized" WCA
    else if (r2 < rmin_*rmin_) ret = Vmin_; // WCA continuation inside rmin
    else ret = V2(r2); // Just the potential outsize of rmin
  return ret;
}  

double Potential1::V0(double r) const ///< The repulsive part of the potential
{
  if(bhFlag_) return (r < getR0() ? V(r) : 0.0);
  return  (r < rmin_ ? V(r)-Vmin_ : 0.0);
}

double Potential1::getHSD(double kT) const ///< Calculates HSD by numeric integration. 
{
  kT_ = kT;
  double hc = getHardCore();

  double rlimit = rmin_;
  if(bhFlag_) rlimit = getR0();
    
  Integrator<const Potential1> II(this, &Potential1::dBH_Kernal, 1e-4, 1e-6);
  return hc + II.integrateFinite(hc,rlimit);
}

double Potential1::getVDW_Parameter(double kT) const ///< Calculates the vDW parameter by numeric integration
{
  kT_ = kT;
  Integrator<const Potential1> II(this, &Potential1::vdw_Kernal, 1e-6, 1e-8);
  if(bhFlag_) return (2*M_PI/kT)*II.integrateFinite(getR0(),rcut_);
  return (2*M_PI/kT)*(II.integrateFinite(r_att_min_,rmin_)+II.integrateFinite(rmin_,rcut_));
}

/////////////////////////////////////////////
///////  LJ potential

LJ::LJ(double sigma, double eps, double rcut) : Potential1(sigma, eps, rcut)
{
  shift_ = (rcut_ < 0 ? 0.0 : vr(rcut_));
  rmin_  = getRmin(); 
  Vmin_  = V(rmin_);
  r0_    = pow(0.5*sqrt(1+shift_)+0.5,-1.0/6.0);          
}

LJ::LJ() : Potential1() {}
  
double LJ::getRmin() const { return pow(2.0,1.0/6.0)*sigma_;}

double LJ::getHardCore() const { return 0.0;}

string LJ::getIdentifier() const
{
    stringstream ss;
    ss << "LJ_" << sigma_ << "_" << eps_ << "_" << rcut_ << r_att_min_ << "_" << bhFlag_;
    return ss.str();    
}
  
double LJ::vr(double r) const
{
  double y = sigma_/r;
  double y3 = y*y*y;
  double y6 = y3*y3;
  return 4*eps_*(y6*y6-y6);
}

double LJ::vr2(double r2) const
{
  double y2 = sigma_*sigma_/r2;
  double y6 = y2*y2*y2;
  return 4*eps_*(y6*y6-y6);
}

/////////////////////////////////////////////
/////// tWF potential

tWF::tWF(double sigma, double eps, double rcut, double alpha) : Potential1(sigma, eps, rcut), alpha_(alpha)
{
  shift_ = (rcut_ <= 0.0 ? 0.0 : vr(rcut_));
  rmin_  = getRmin();
  Vmin_  = V(rmin_);
  r0_    = sqrt(1+pow(25*sqrt(1+shift_)+25,-1.0/3.0));
}

tWF::tWF() : Potential1() {}

double tWF::getRmin() const { return sigma_*sqrt(1+pow(2.0/alpha_,1.0/3.0));}

double tWF::getHardCore() const { return sigma_;}

string tWF::getIdentifier() const
{
  stringstream ss;
  ss << "TWF_" << sigma_ << "_" << eps_ << "_" << alpha_ << "_" << rcut_ << r_att_min_ << "_" << bhFlag_;
  return ss.str();    
}

double tWF::vr(double r) const
{
  if(r < sigma_) return 1e50;
  
  double s = r/sigma_;
  double y = 1.0/(s*s-1);
  double y3 = y*y*y;
  
  return (4*eps_/(alpha_*alpha_))*(y3*y3-alpha_*y3);
}

double tWF::vr2(double r2) const
{
  if(r2 < sigma_*sigma_) return 1e50;
  
  double s2 = r2/(sigma_*sigma_);
  double y = 1.0/(s2-1);
  double y3 = y*y*y;
  
  return (4*eps_/(alpha_*alpha_))*(y3*y3-alpha_*y3);
}  

/////////////////////////////////////////////
/////// WHDF potential

WHDF::WHDF(double sigma, double eps, double rcut) : Potential1(sigma, eps, rcut)
{
  eps_rescaled_ = eps_*2*pow(rcut/sigma,2)*pow(2*((rcut/sigma)*(rcut/sigma)-1)/3.0,-3.0);
  
  shift_ = 0.0; // cutoff is built in
  rmin_  = getRmin(); 
  Vmin_  = V(rmin_);
  r0_    = 1.0;
}

WHDF::WHDF() : Potential1() {}
  
double WHDF::getRmin()       const { return rcut_*pow((1.0+2.0*(rcut_/sigma_)*(rcut_/sigma_))/3.0,-0.5);}
double WHDF::getHardCore()   const { return 0.0;}
string WHDF::getIdentifier() const
{
  stringstream ss;
  ss << "WHDF_" << sigma_ << "_" << eps_ << "_" << rcut_ << r_att_min_ << "_" << bhFlag_;
  return ss.str();    
}

double WHDF::vr(double r) const
{
  double y = sigma_/r;
  double z = rcut_/r;
  return eps_rescaled_*(y*y-1)*(z*z-1)*(z*z-1);
}

double WHDF::vr2(double r2) const
{
  double y2 = sigma_*sigma_/r2;
  double z2 = rcut_*rcut_/r2;
  return eps_rescaled_*(y2-1)*(z2-1)*(z2-1);
}  
