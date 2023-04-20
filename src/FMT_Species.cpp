#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>

#include <gsl/gsl_integration.h>

using namespace std;

#ifdef USE_MPI
#include <mpi.h>  
#endif

#include "Species.h"
#include "myColor.h"

const double dmin = SMALL_VALUE;
const double SMALL_VALUE_FOR_J_INTEGRAL = 1e-14; //value does not affect results

// TODO: Declare constant --> should be passed as an argument to G functions
double SMALL_VALUE_FOR_WORKAROUND_NAN_IN_G_FUNCTION = 1e-8;

#ifdef DEBUG_FMT_WEIGHTS
  bool USE_OLD_WORKAROUND_NAN_IN_G_FUNCTION = false;
  bool SET_J_INTEGRAL_TO_ZERO_WHEN_V_ON_BOUNDARY = false;
#else
  #if defined FMT_WEIGHTS_BEFORE_JUN_2021
    const bool USE_OLD_WORKAROUND_NAN_IN_G_FUNCTION = true; //small value will be 1e-15
    const bool SET_J_INTEGRAL_TO_ZERO_WHEN_V_ON_BOUNDARY = false;
  #elif defined FMT_WEIGHTS_BEFORE_JAN_2023
    const bool USE_OLD_WORKAROUND_NAN_IN_G_FUNCTION = true; //small value will be 1e-15
    const bool SET_J_INTEGRAL_TO_ZERO_WHEN_V_ON_BOUNDARY = false;
  #else
    const bool USE_OLD_WORKAROUND_NAN_IN_G_FUNCTION = false;
    const bool SET_J_INTEGRAL_TO_ZERO_WHEN_V_ON_BOUNDARY = false; // true should guarantee no NAN
  #endif
#endif


// These impose density limits based on the fact that in the most extreme case,
// if the density at some point is dmax and all other points it is dmin, then
// eta at that point is (4*pi/3)*hsd*hsd*hsd*dmin + (dmax-dmin)*dV.
// So the largest possible value of dmax, based on eta < 1, is
//    dmax = dmin + (eta-etamin)/dV
// Note that the limits used here are really only valid for FMT_Species
// but for now, those are really the only ones used.

void FMT_Species::set_density_from_alias(const DFT_Vec &x)
{
  //Species::set_density_from_alias(x);
  
  long pos;
  const double etamin = dmin*(4*M_PI*getHSD()*getHSD()*getHSD()/3);
  const double c = (1.0-etamin)/density_->dV();
  //const double c = (0.99-etamin)/density_->dV();
  
  #ifdef USE_OMP    
  #pragma omp parallel for  private(pos)  schedule(static)
  #endif  
  for(pos=0;pos<x.size();pos++)
  {
    double y = x.get(pos);
    //double z = dmin +c*(1-exp(-y*y));
    double z = dmin +c*y*y/(1+y*y);
    density_->set(pos,z);
  }
}

void FMT_Species::get_density_alias(DFT_Vec &x) const
{
  //Species::get_density_alias(x);
  x.zeros(density_->size());
   
  long pos;
  const double etamin = dmin*(4*M_PI*getHSD()*getHSD()*getHSD()/3);
  const double c = (1.0-etamin)/density_->dV();
  //const double c = (0.99-etamin)/density_->dV();

  #ifdef USE_OMP
  #pragma omp parallel for  private(pos)  schedule(static)
  #endif
  for(pos=0;pos<x.size();pos++)
  {
    double z = (density_->get(pos) - dmin)/c;
    //double y = sqrt(fabs(log(1.0/(1-z))));
    double y = sqrt(fabs(z/(1-z)));
    x.set(pos, y);          
  }
}

void FMT_Species::convert_to_alias_deriv(DFT_Vec &dF_dRho) const
{
  DFT_Vec x; get_density_alias(x);
  convert_to_alias_deriv(x, dF_dRho);
}

void FMT_Species::convert_to_alias_deriv(DFT_Vec &x, DFT_Vec &dF_dRho) const
{
  //Species::convert_to_alias_deriv(x,dF_dRho);  
  
  long pos;
  const double etamin = dmin*(4*M_PI*getHSD()*getHSD()*getHSD()/3);
  const double c = (1.0-etamin)/density_->dV();
  //  const double c = (0.99-etamin)/density_->dV();  

  #ifdef USE_OMP
  #pragma omp parallel for  private(pos)  schedule(static)
  #endif
  for(pos=0;pos<x.size();pos++)
  {
    double y = x.get(pos);
    double df = dF_dRho.get(pos);
    //dF_dRho.set(pos, df*(2*c*y*exp(-y*y)));
    dF_dRho.set(pos, df*(2*c*y/((1+y*y)*(1+y*y))));
  }
}

void FMT_Species::convert_to_alias_increment(DFT_Vec &dRho) const
{
  DFT_Vec x; get_density_alias(x);
  convert_to_alias_increment(x, dRho);
}

void FMT_Species::convert_to_alias_increment(DFT_Vec &x, DFT_Vec &dRho) const
{
  //Species::convert_to_alias_increment(x,dF_dRho);  
  
  long pos;
  const double etamin = dmin*(4*M_PI*getHSD()*getHSD()*getHSD()/3);
  const double c = (1.0-etamin)/density_->dV();
  //  const double c = (0.99-etamin)/density_->dV();  

  #ifdef USE_OMP
  #pragma omp parallel for  private(pos)  schedule(static)
  #endif
  for(pos=0;pos<x.size();pos++)
  {
    double y = x.get(pos);
    double drho = dRho.get(pos);
    //dRho.set(pos, drho/(2*c*y*exp(-y*y)));
    dRho.set(pos, drho*(1+y*y)*(1+y*y)/(2*c*y));
  }
}

void FMT_Species::get_second_derivatives_of_density_wrt_alias(DFT_Vec &d2Rhodx2) const 
{
  d2Rhodx2.zeros(density_->size());
  
  long pos;
  const double etamin = dmin*(4*M_PI*getHSD()*getHSD()*getHSD()/3);
  const double c = (1.0-etamin)/density_->dV();
  //  const double c = (0.99-etamin)/density_->dV();  

  #ifdef USE_OMP
  #pragma omp parallel for  private(pos)  schedule(static)
  #endif
  for(pos=0;pos<density_->size();pos++)
  {
    double z = (density_->get(pos) - dmin)/c;
    d2Rhodx2.set(pos, 2*c*(1+2*z)*(1-z)*(1-z));
  }
}

// Trivial alias
/*
void FMT_Species::set_density_from_alias(const DFT_Vec &x) {density_->set(x);}
void FMT_Species::get_density_alias(DFT_Vec &x) const {x = density_->get_density_real();}
void FMT_Species::convert_to_alias_deriv(DFT_Vec &x, DFT_Vec &dF_dRho) const {;}
void FMT_Species::convert_to_alias_increment(DFT_Vec &x, DFT_Vec &dF_dRho) const {;}
void FMT_Species::convert_to_alias_deriv(DFT_Vec &dF_dRho) const {;}
void FMT_Species::convert_to_alias_increment(DFT_Vec &dF_dRho) const {;}
void FMT_Species::get_second_derivatives_of_density_wrt_alias(DFT_Vec &d2Rhodx2) const {d2Rhodx2.zeros(density_->size()); d2Rhodx2.add(1.0);}
*/


// Identical transformation if jacobian is diagonal (local alias)
void FMT_Species::convert_to_density_increment(DFT_Vec &dRho) const {convert_to_alias_deriv(dRho);}
void FMT_Species::convert_to_density_increment(DFT_Vec &x, DFT_Vec &dRho) const {convert_to_alias_deriv(x, dRho);}


FMT_Species::FMT_Species(Density& density, double hsd, double mu, bool verbose, int seq): Species(density,mu,seq), hsd_(hsd), fmt_weighted_densities(11)
{
  verbose_ = verbose;
  Initialize();
}

void FMT_Species::Initialize()
{
  #ifdef DEBUG_FMT_WEIGHTS
    Check();
  #endif
  
  long Nx = density_->Nx();
  long Ny = density_->Ny();
  long Nz = density_->Nz();

  verbose_ = false;
  
  for(FMT_Weighted_Density &d: fmt_weighted_densities)
    d.initialize(Nx, Ny, Nz);

  #if defined FMT_WEIGHTS_BEFORE_JUN_2021
    // The old function "generateWeights_before_jun_2021" could be removed 
    // as the new one does reproduce the same calculation to machine precision
    // by using here the arguments on the commented line below.
    generateWeights_before_jun_2021(hsd_, fmt_weighted_densities);
    //generateWeights(hsd_, fmt_weighted_densities, density_->getDX(), false);
  #elif defined FMT_WEIGHTS_BEFORE_JAN_2023
    // The old function "generateWeights_before_jan_2023" could be removed 
    // as the new one does reproduce the same calculation to machine precision
    // by using here the arguments on the commented line below.
    generateWeights_before_jan_2023(hsd_, fmt_weighted_densities);
    //generateWeights(hsd_, fmt_weighted_densities, 1.0, false);
  #else
    generateWeights(hsd_, fmt_weighted_densities);
  #endif

  for(FMT_Weighted_Density &d: fmt_weighted_densities)
    d.transformWeights();
}


#ifdef DEBUG_FMT_WEIGHTS
void FMT_Species::Check()
{
  long Nx = density_->Nx();
  long Ny = density_->Ny();
  long Nz = density_->Nz();
  
  // New code version
  vector<FMT_Weighted_Density> fmt_weighted_densities_new(11);
  for(FMT_Weighted_Density &d: fmt_weighted_densities_new) d.initialize(Nx, Ny, Nz);
  generateWeights(hsd_, fmt_weighted_densities_new);
  
  // New code with parameters set to reproduce old calculations
  USE_OLD_WORKAROUND_NAN_IN_G_FUNCTION = true;
  
  vector<FMT_Weighted_Density> fmt_weighted_densities_new_as_before_jun_2021(11);
  for(FMT_Weighted_Density &d: fmt_weighted_densities_new_as_before_jun_2021) d.initialize(Nx, Ny, Nz);
  generateWeights(hsd_, fmt_weighted_densities_new_as_before_jun_2021, density_->getDX(), false);
  
  vector<FMT_Weighted_Density> fmt_weighted_densities_new_as_before_jan_2023(11);
  for(FMT_Weighted_Density &d: fmt_weighted_densities_new_as_before_jan_2023) d.initialize(Nx, Ny, Nz);
  generateWeights(hsd_, fmt_weighted_densities_new_as_before_jan_2023, 1.0, false);
  
  USE_OLD_WORKAROUND_NAN_IN_G_FUNCTION = false;
  
  //////////////////////////////////////////////////////////////////////////////
  // This can go if we decide to remove the old functions
  // "generateWeights_before_jun_2021" and "generateWeights_before_jan_2023"
  
  // Old code versions
  #ifdef FMT_WEIGHTS_BEFORE_JUN_2021
    USE_OLD_WORKAROUND_NAN_IN_G_FUNCTION = true;
    
    vector<FMT_Weighted_Density> fmt_weighted_densities_before_jun_2021(11);
    for(FMT_Weighted_Density &d: fmt_weighted_densities_before_jun_2021) d.initialize(Nx, Ny, Nz);
    generateWeights_before_jun_2021(hsd_, fmt_weighted_densities_before_jun_2021);
    
    USE_OLD_WORKAROUND_NAN_IN_G_FUNCTION = false;
    
    cout << endl;
    cout << scientific << setprecision(2);
    cout << "===========================================================" << endl;
    cout << "  Testing FMT weights calculation: new vs old (<jun 2021)  " << endl;
    cout << endl;
    
    Check_Print(fmt_weighted_densities_new, fmt_weighted_densities_before_jun_2021); cout << endl;
    
    cout << endl;
    cout << scientific << setprecision(2);
    cout << "===========================================================" << endl;
    cout << "  Testing FMT weights calculation: new vs old (<jun 2021)  " << endl;
    cout << "  Reproducing old calculation using new code scaled, the   " << endl;
    cout << "  old G-function regulator (1e-15) and keep contribution   " << endl;
    cout << endl;
    
    Check_Print(fmt_weighted_densities_new_as_before_jun_2021, fmt_weighted_densities_before_jun_2021); cout << endl;
  #endif
  
  #ifdef FMT_WEIGHTS_BEFORE_JAN_2023
    USE_OLD_WORKAROUND_NAN_IN_G_FUNCTION = true;
    
    vector<FMT_Weighted_Density> fmt_weighted_densities_before_jan_2023(11);
    for(FMT_Weighted_Density &d: fmt_weighted_densities_before_jan_2023) d.initialize(Nx, Ny, Nz);
    generateWeights_before_jan_2023(hsd_, fmt_weighted_densities_before_jan_2023);
    
    USE_OLD_WORKAROUND_NAN_IN_G_FUNCTION = false;
    
    cout << endl;
    cout << scientific << setprecision(2);
    cout << "===========================================================" << endl;
    cout << "  Testing FMT weights calculation: new vs old (<jan 2023)  " << endl;
    cout << endl;
    
    Check_Print(fmt_weighted_densities_new, fmt_weighted_densities_before_jan_2023); cout << endl;
    
    cout << endl;
    cout << scientific << setprecision(2);
    cout << "===========================================================" << endl;
    cout << "  Testing FMT weights calculation: new vs old (<jan 2023)  " << endl;
    cout << "  Reproducing old calculation using new code unscaled, the " << endl;
    cout << "  old G-function regulator (1e-15) and keep contribution   " << endl;
    cout << endl;
    
    Check_Print(fmt_weighted_densities_new_as_before_jan_2023, fmt_weighted_densities_before_jan_2023); cout << endl;
  #endif
  //////////////////////////////////////////////////////////////////////////////
  
  // Modification of G-function regulator
  USE_OLD_WORKAROUND_NAN_IN_G_FUNCTION = true;
  
  vector<FMT_Weighted_Density> fmt_weighted_densities_new_code_but_old_regulator(11);
  for(FMT_Weighted_Density &d: fmt_weighted_densities_new_code_but_old_regulator) d.initialize(Nx, Ny, Nz);
  generateWeights(hsd_, fmt_weighted_densities_new_code_but_old_regulator);
  
  USE_OLD_WORKAROUND_NAN_IN_G_FUNCTION = false;
  
  cout << endl;
  cout << scientific << setprecision(2);
  cout << "=======================================================================" << endl;
  cout << "  Testing FMT weights calculation: old vs new (handling of regulator)  " << endl;
  cout << "  old: regulator R += " << setw(8) << 1e-15 << ", not subtracting contribution to G  " << endl;
  cout << "  new: regulator R *= (1+" << setw(8) << SMALL_VALUE_FOR_WORKAROUND_NAN_IN_G_FUNCTION << "), subtracting contribution to G  " << endl;
  cout << "  This is what changed before and after January 2023" << endl;
  cout << endl;
  
  Check_Print(fmt_weighted_densities_new, fmt_weighted_densities_new_code_but_old_regulator); cout << endl;
  
  // Scaling under old G-function regulator
  cout << endl;
  cout << scientific << setprecision(2);
  cout << "==========================================================" << endl;
  cout << "  Testing FMT weights calculation: /=1 vs /=dx (scaling)  " << endl;
  cout << "  Using old regulator R += " << setw(8) << 1e-15 << ", not subtracting contribution to G  " << endl;
  cout << "  This is what changed before and after June 2021" << endl;
  cout << endl;
  
  Check_Print(fmt_weighted_densities_new_as_before_jan_2023, fmt_weighted_densities_new_as_before_jun_2021); cout << endl;
  
  // Scaling under new G-function regulator
  vector<FMT_Weighted_Density> fmt_weighted_densities_new_unscaled(11);
  for(FMT_Weighted_Density &d: fmt_weighted_densities_new_unscaled) d.initialize(Nx, Ny, Nz);
  generateWeights(hsd_, fmt_weighted_densities_new_unscaled, 1.0, true);
  
  vector<FMT_Weighted_Density> fmt_weighted_densities_new_scaled(11);
  for(FMT_Weighted_Density &d: fmt_weighted_densities_new_scaled) d.initialize(Nx, Ny, Nz);
  generateWeights(hsd_, fmt_weighted_densities_new_scaled, density_->getDX(), true);
  
  cout << endl;
  cout << scientific << setprecision(2);
  cout << "==========================================================" << endl;
  cout << "  Testing FMT weights calculation: /=1 vs /=dx (scaling)  " << endl;
  cout << "  Using new regulator R *= (1+" << setw(8) << SMALL_VALUE_FOR_WORKAROUND_NAN_IN_G_FUNCTION << "), subtracting contribution to G  " << endl;
  cout << endl;
  
  Check_Print(fmt_weighted_densities_new_unscaled, fmt_weighted_densities_new_scaled); cout << endl;
  
  // Contribution of the G-function regulator
  SMALL_VALUE_FOR_WORKAROUND_NAN_IN_G_FUNCTION *= 2;
  
  vector<FMT_Weighted_Density> fmt_weighted_densities_new_Gx2(11);
  for(FMT_Weighted_Density &d: fmt_weighted_densities_new_Gx2) d.initialize(Nx, Ny, Nz);
  generateWeights(hsd_, fmt_weighted_densities_new_Gx2, 1.0, false);
  
  SMALL_VALUE_FOR_WORKAROUND_NAN_IN_G_FUNCTION /= 2;
  
  cout << endl;
  cout << scientific << setprecision(2);
  cout << "=====================================================================" << endl;
  cout << "  Testing FMT weights calculation: new vs new (G: small value x 2)  " << endl;
  cout << "  Using new regulator R *= (1+" << setw(8) << SMALL_VALUE_FOR_WORKAROUND_NAN_IN_G_FUNCTION << "), scale=1 and NOT subtracting contribution to G  " << endl;
  cout << endl;
  
  Check_Print(fmt_weighted_densities_new, fmt_weighted_densities_new_Gx2); cout << endl;
  
  // Contribution of the G-function regulator (after subtraction of 1st order perturbation)
  SMALL_VALUE_FOR_WORKAROUND_NAN_IN_G_FUNCTION *= 2;
  
  fmt_weighted_densities_new_Gx2.clear();
  fmt_weighted_densities_new_Gx2 = vector<FMT_Weighted_Density>(11);
  for(FMT_Weighted_Density &d: fmt_weighted_densities_new_Gx2) d.initialize(Nx, Ny, Nz);
  generateWeights(hsd_, fmt_weighted_densities_new_Gx2, 1.0, true);
  
  SMALL_VALUE_FOR_WORKAROUND_NAN_IN_G_FUNCTION /= 2;
  
  cout << endl;
  cout << scientific << setprecision(2);
  cout << "=====================================================================" << endl;
  cout << "  Testing FMT weights calculation: new vs new (G: small value x 2)  " << endl;
  cout << "  Using new regulator R *= (1+" << setw(8) << SMALL_VALUE_FOR_WORKAROUND_NAN_IN_G_FUNCTION << "), scale=1 and subtracting contribution to G  " << endl;
  cout << endl;
  
  Check_Print(fmt_weighted_densities_new, fmt_weighted_densities_new_Gx2); cout << endl;
  
  // Setting or not setting J(X,V) to zero when V=0
  const bool default_policy = SET_J_INTEGRAL_TO_ZERO_WHEN_V_ON_BOUNDARY;
  SET_J_INTEGRAL_TO_ZERO_WHEN_V_ON_BOUNDARY = true;
  
  vector<FMT_Weighted_Density> fmt_weighted_densities_new_J_zero(11);
  for(FMT_Weighted_Density &d: fmt_weighted_densities_new_J_zero) d.initialize(Nx, Ny, Nz);
  generateWeights(hsd_, fmt_weighted_densities_new_J_zero, 1.0, true);
  
  SET_J_INTEGRAL_TO_ZERO_WHEN_V_ON_BOUNDARY = default_policy;
  
  cout << endl;
  cout << scientific << setprecision(2);
  cout << "========================================================================" << endl;
  cout << "  Testing FMT weights calculation: new with vs without J(X,0) set to 0  " << endl;
  cout << "  Using new regulator R *= (1+" << setw(8) << SMALL_VALUE_FOR_WORKAROUND_NAN_IN_G_FUNCTION << "), scale=1 and subtracting contribution to G  " << endl;
  cout << endl;
  
  Check_Print(fmt_weighted_densities_new, fmt_weighted_densities_new_J_zero); cout << endl;
}


void FMT_Species::Check_Print(const vector<FMT_Weighted_Density> &fmt_weighted_densities_old, 
                              const vector<FMT_Weighted_Density> &fmt_weighted_densities_new)
{
  for(int i=0; i<11; i++)
  {
    const FMT_Weighted_Density &d_new = fmt_weighted_densities_new[i];
    const FMT_Weighted_Density &d_old = fmt_weighted_densities_old[i];
    
    if      (i== EI()    ) cout << "eta   ";
    else if (i== SI()    ) cout << "s     ";
    else if (i== VI(0)   ) cout << "vx    ";
    else if (i== VI(1)   ) cout << "vy    ";
    else if (i== VI(2)   ) cout << "vz    ";
    else if (i== TI(0,0) ) cout << "txx   ";
    else if (i== TI(1,1) ) cout << "tyy   ";
    else if (i== TI(2,2) ) cout << "tzz   ";
    else if (i== TI(0,1) ) cout << "txy   ";
    else if (i== TI(0,2) ) cout << "txz   ";
    else if (i== TI(1,2) ) cout << "tyz   ";
    
    double norm = 0.0;
    double norm_diff = 0.0;
    
    for (long p=0; p<density_->Ntot(); p++)
    {
      norm += pow(d_new.getWeight(p),2);
      norm_diff += pow(d_new.getWeight(p)-d_old.getWeight(p),2);
    }
    
    norm      = sqrt(norm      /density_->Ntot());
    norm_diff = sqrt(norm_diff /density_->Ntot());
    
    cout << "    norm = " << setw(10) << norm;
    cout << "    norm_diff = " << setw(10) << norm_diff;
    cout << "    norm_diff/norm = " << setw(10) << norm_diff/norm;
    cout << endl;
  }
}
#endif // DEBUG_FMT_WEIGHTS

/*
double my_asin(double x)
{
  if (fabs(x)-1>SMALL_VALUE_FOR_J_INTEGRAL) 
  {
    cout << "|x|-1 = " << fabs(x)-1 << endl;
    throw runtime_error("FMT_Species::my_asin was given argument |x|>1");
  }
  else
  {
    x = max(-1.0,min(1.0,x));
  }
  
  //return asin(x);
  
  // Trying to improve the evaluation of asin(x) when x near 1.
  // Is it really why the evaluations of J(X,V) cause problems?
  
  double y = 1-x;
  double z = 1+x;
  
  if      (fabs(y)<SMALL_VALUE_FOR_J_INTEGRAL) return  M_PI/2;
  else if (fabs(z)<SMALL_VALUE_FOR_J_INTEGRAL) return -M_PI/2;
  else return asin(x);
  
  
  //if (y<SMALL_VALUE_FOR_ASIN)
  //{
  //  double phi = sqrt(4*y-2*y*y);
  //  cout << "Asin:    y = " << y << "    PI/2-asin(x) = " << M_PI/2-asin(x)  << "    PI/2-approx = " << phi << endl;
  //  return M_PI/2-phi;
  //}
  //else if (z<SMALL_VALUE_FOR_ASIN)
  //{
  //  double phi = sqrt(-4*z-2*z*z);
  //  cout << "Asin:    z = " << z << "    PI/2+asin(x) = " << M_PI/2+asin(x)  << "    PI/2+approx = " << phi << endl;
  //  return -M_PI/2+phi;
  //}
  //else
  //{
  //  return asin(x);
  //}
  
  //if (x>1/sqrt(2)) // I expect more accuracy with acos --> no
  //{
  //  double y = 1-x;
  //  return acos(sqrt(2*y-y*y));
  //}
  //if (x<-1/sqrt(2)) // I expect more accuracy with acos --> no
  //{
  //  double y = 1+x;
  //  return acos(sqrt(-2*y-y*y));
  //}
  //else
  //{
  //  return asin(x);
  //}
}
*/


static void getI(double X, double A, double I[])
{
  //double a = asin(X/A);
  double a = asin(max(-1.0,min(1.0,X/A)));
  //double a = my_asin(X/A);
  double b = sqrt(fabs(A*A-X*X));
  
  I[0] = 0.5*X*b+0.5*A*A*a;
  I[1] = -(1.0/3)*(A*A-X*X)*b;
  I[2] = 0.125*X*(2*X*X-A*A)*b+0.125*A*A*A*A*a;
  I[3] = -(1.0/15)*(A*A-X*X)*(3*X*X+2*A*A)*b;
  I[4] = (1.0/48)*(8*X*X*X*X-2*A*A*X*X-3*A*A*A*A)*X*b+0.0625*A*A*A*A*A*A*a;
}
  
static void getJ(double X, double V, double R, double J[])
{
  // Aren't they supposed to be equal to zero when V=0? 
  // Why is there a R^4 term in the SI of the PRE 2020?
  if (fabs(V/R)<SMALL_VALUE_FOR_J_INTEGRAL && SET_J_INTEGRAL_TO_ZERO_WHEN_V_ON_BOUNDARY)
  {
    for (int i=0; i<5; i++) J[i]=0.0; //J[2] = (M_PI/8)*R*R*R*R;
    return;
  }
  
  double aV  = asin(max(-1.0,min(1.0,V/max(SMALL_VALUE_FOR_J_INTEGRAL,sqrt(R*R-X*X)))));
  double aX  = asin(max(-1.0,min(1.0,X/sqrt(R*R-V*V))));
  double aVX = asin(max(-1.0,min(1.0,X*V/max(SMALL_VALUE_FOR_J_INTEGRAL,sqrt((R*R-V*V)*(R*R-X*X))))));
  
  //double aV  = my_asin(V/sqrt(R*R-X*X));
  //double aX  = my_asin(X/sqrt(R*R-V*V));
  //double aVX = my_asin(V*X/sqrt((R*R-V*V)*(R*R-X*X)));
  
  double b   = sqrt(fabs(R*R-V*V-X*X));
  
  J[0] = X*aV+V*aX-R*aVX;
  J[1] = -0.5*V*b+0.5*(X*X-R*R)*aV;
  J[2] = -(1.0/6)*V*(V*V-3*R*R)*aX+(1.0/3)*X*X*X*aV-(1.0/3)*R*R*R*aVX-(1.0/6)*V*X*b;
  J[3] = 0.25*(X*X*X*X-R*R*R*R)*aV+(M_PI/8)*R*R*R*R-(V/12.0)*(5*R*R+X*X-2*V*V)*b;
  J[4] = 0.025*X*V*(3*V*V-2*X*X-7*R*R)*b+0.025*V*(3*V*V*V*V-10*V*V*R*R+15*R*R*R*R)*aX+0.2*X*X*X*X*X*aV-0.2*R*R*R*R*R*aVX;
  
  if(std::isnan(J[1]) || std::isnan(J[2]) || std::isnan(J[3]) || std::isnan(J[4]))
    {
      cout << scientific << setprecision(2);
      cout << X << " " << V << " " << R << " | " << b << " " << R*R-V*V << " " << R*R-X*X << " " << R*R - V*V -X*X << endl;
      throw std::runtime_error("done");
    }
  
  // Trying to identify where inaccuracies come from
  //J[0] = J[1] = J[2] = J[3] = J[4] = 0; //ok (the inaccuracies come only from the evaluations of J(X,V))
  //J[0] = J[2] = J[4] = 0; //err
  //J[1] = J[3] = 0; //err
  //       J[1] = J[2] = J[3] = J[4] = 0; //ok  except eta,s
  //J[0] =        J[2] = J[3] = J[4] = 0; //err except txx,tyy,tzz
  //J[0] = J[1]        = J[3] = J[4] = 0; //err except s
  //J[0] = J[1] = J[2] =        J[4] = 0; //err except s,vx,vy,vz
  //J[0] = J[1] = J[2] = J[3] =        0; //ok  except txy,txz,tyz
  //---> basically the calculations are fine unless there is a J(X,V) involved
}


double G_eta(double R, double X, double Vy, double Vz, double Tx, double Ty, double Tz)
{
  double Iy[5];
  double Iz[5];
  double Jy[5];
  double Jz[5];

  // A temporary workaround to avoid NAN.
  if (USE_OLD_WORKAROUND_NAN_IN_G_FUNCTION) R += 1e-15;
  else R *= (1+SMALL_VALUE_FOR_WORKAROUND_NAN_IN_G_FUNCTION);
  
  getI(X, sqrt(R*R-Vy*Vy), Iy);
  getI(X, sqrt(R*R-Vz*Vz), Iz);
  getJ(X,Vy,R,Jy);
  getJ(X,Vz,R,Jz);
  
  double A = Ty*Tz*Vy*Vz+0.25*Vy*Vy*Vz*Vz+0.125*R*R*R*R
    -(1.0/6)*(Ty*Vy*Vy*Vy+Tz*Vz*Vz*Vz) -0.5*Tz*Vz*Vy*Vy-0.5*Ty*Vy*Vz*Vz
    +0.125*Vy*Vy*Vy*Vy+0.125*Vz*Vz*Vz*Vz+0.5*R*R*(Ty*Vy+Tz*Vz-0.5*Vy*Vy-0.5*Vz*Vz+0.5*M_PI*Ty*Tz);
  double B = 0.5*Ty*Vy+0.5*Tz*Vz-0.25*Vy*Vy-0.25*Vz*Vz+(M_PI/4)*Ty*Tz+0.25*R*R;
  double C = 0.5*Ty*Vy+(1.0/3)*(R*R-Vy*Vy);
  double D = 0.5*Tz*Vz+(1.0/3)*(R*R-Vz*Vz);

  double g = Tx*A*X-0.5*A*X*X-(1.0/3)*B*Tx*X*X*X+0.25*X*X*X*X*B
    +0.025*Tx*X*X*X*X*X-(1.0/48)*X*X*X*X*X*X;

  if(std::isnan(g))
    throw std::runtime_error("Here 1");
  
  g -= 0.5*Ty*Tz*(Tx*R*R*Jy[0]-R*R*Jy[1]-Tx*Jy[2]+Jy[3]);
  if(std::isnan(g))
    throw std::runtime_error("Here 2");  
  g -= 0.5*Ty*Tz*(Tx*R*R*Jz[0]-R*R*Jz[1]-Tx*Jz[2]+Jz[3]);
  if(std::isnan(g))
    {
      cout << Jz[0] << " " << Jz[1] << " " << Jz[2] << " " << Jz[3] << endl;
      throw std::runtime_error("Here 3");
    }
  g -= (Tx*Tz*C*Iy[0]-Tz*C*Iy[1]-(1.0/3)*Tx*Tz*Iy[2]+(1.0/3)*Tz*Iy[3]);
  if(std::isnan(g))
    throw std::runtime_error("Here 4");  
  g -= (Tx*Ty*D*Iz[0]-Ty*D*Iz[1]-(1.0/3)*Tx*Ty*Iz[2]+(1.0/3)*Ty*Iz[3]);


  
  return g;
}

double G_s(double R, double X, double Vy, double Vz, double Tx, double Ty, double Tz)
{
  double Iy[5];
  double Iz[5];
  double Jy[5];
  double Jz[5];

  // A temporary workaround to avoid NAN.
  if (USE_OLD_WORKAROUND_NAN_IN_G_FUNCTION) R += 1e-15;
  else R *= (1+SMALL_VALUE_FOR_WORKAROUND_NAN_IN_G_FUNCTION);
  
  getI(X, sqrt(R*R-Vy*Vy), Iy);
  getI(X, sqrt(R*R-Vz*Vz), Iz);
  getJ(X,Vy,R,Jy);
  getJ(X,Vz,R,Jz);

  double A = Ty*Vy+Tz*Vz-0.5*Vy*Vy-0.5*Vz*Vz+0.5*R*R+(M_PI/2)*Ty*Tz;
  
  double g = 0.125*R*X*X*X*X-(1.0/6)*R*Tx*X*X*X-0.5*R*A*X*X+R*A*Tx*X;

  g -= R*Tz*(Tx*Iy[0]-Iy[1]);
  g -= R*Ty*(Tx*Iz[0]-Iz[1]);

  g -= R*Ty*Tz*(Tx*Jy[0]-Jy[1]);
  g -= R*Ty*Tz*(Tx*Jz[0]-Jz[1]);

  return g;
}

double G_vx(double R, double X, double Vy, double Vz, double Tx, double Ty, double Tz)
{
  double Iy[5];
  double Iz[5];
  double Jy[5];
  double Jz[5];

  // A temporary workaround to avoid NAN.
  if (USE_OLD_WORKAROUND_NAN_IN_G_FUNCTION) R += 1e-15;
  else R *= (1+SMALL_VALUE_FOR_WORKAROUND_NAN_IN_G_FUNCTION);  

  getI(X, sqrt(R*R-Vy*Vy), Iy);
  getI(X, sqrt(R*R-Vz*Vz), Iz);
  getJ(X,Vy,R,Jy);
  getJ(X,Vz,R,Jz);

  double A = Ty*Vy+Tz*Vz-0.5*Vy*Vy-0.5*Vz*Vz+0.5*R*R+(M_PI/2)*Ty*Tz;
  
  double g = 0.1*R*X*X*X*X*X-0.125*R*Tx*X*X*X*X-(1.0/3)*R*A*X*X*X+0.5*R*A*Tx*X*X;

  g -= R*Tz*(Tx*Iy[1]-Iy[2]);
  g -= R*Ty*(Tx*Iz[1]-Iz[2]);

  g -= R*Ty*Tz*(Tx*Jy[1]-Jy[2]);
  g -= R*Ty*Tz*(Tx*Jz[1]-Jz[2]);

  
  return g;
}

double G_txx(double R, double X, double Vy, double Vz, double Tx, double Ty, double Tz)
{
  double Iy[5];
  double Iz[5];
  double Jy[5];
  double Jz[5];

  // A temporary workaround to avoid NAN.
  if (USE_OLD_WORKAROUND_NAN_IN_G_FUNCTION) R += 1e-15;
  else R *= (1+SMALL_VALUE_FOR_WORKAROUND_NAN_IN_G_FUNCTION);
  
  getI(X, sqrt(R*R-Vy*Vy), Iy);
  getI(X, sqrt(R*R-Vz*Vz), Iz);
  getJ(X,Vy,R,Jy);
  getJ(X,Vz,R,Jz);

  double A = Ty*Vy+Tz*Vz-0.5*Vy*Vy-0.5*Vz*Vz+0.5*R*R+(M_PI/2)*Ty*Tz;
  
  double g = (1.0/12)*R*X*X*X*X*X*X-0.1*R*Tx*X*X*X*X*X-0.25*R*A*X*X*X*X+(1.0/3)*R*A*Tx*X*X*X;

  g -= R*Tz*(Tx*Iy[2]-Iy[3]);
  g -= R*Ty*(Tx*Iz[2]-Iz[3]);

  g -= R*Ty*Tz*(Tx*Jy[2]-Jy[3]);
  g -= R*Ty*Tz*(Tx*Jz[2]-Jz[3]);

  return g;
}

double G_txy(double R, double X, double Vy, double Vz, double Tx, double Ty, double Tz)
{
  double Iy[5];
  double Iz[5];
  double Jy[5];
  double Jz[5];

  // A temporary workaround to avoid NAN.
  if (USE_OLD_WORKAROUND_NAN_IN_G_FUNCTION) R += 1e-15;
  else R *= (1+SMALL_VALUE_FOR_WORKAROUND_NAN_IN_G_FUNCTION);
  
  getI(X, sqrt(R*R-Vy*Vy), Iy);
  getI(X, sqrt(R*R-Vz*Vz), Iz);
  getJ(X,Vy,R,Jy);
  getJ(X,Vz,R,Jz);

  double A = (1.0/6)*(3*Ty*R*R+2*Vy*Vy*Vy-3*Ty*Vy*Vy-3*Ty*Vz*Vz+6*Ty*Tz*Vz+1.5*M_PI*R*R*Tz);
  double B = 0.25*M_PI*Tz+0.5*Ty;
  double C = -(1.0/6)*(2*R*R-2*Vz*Vz+3*Tz*Vz);
  
  double g = -0.5*R*A*Tx*X*X+(1.0/3)*R*A*X*X*X+0.25*R*B*Tx*X*X*X*X-0.2*R*B*X*X*X*X*X;

  g += 0.5*R*Tz*(Tx*R*R*Jy[1]-R*R*Jy[2]-Tx*Jy[3]+Jy[4]);
  g += 0.5*R*Tz*(Tx*R*R*Jz[1]-R*R*Jz[2]-Tx*Jz[3]+Jz[4]);

  g += 0.5*R*Tz*(2*Ty-Vy)*(Tx*Iy[1]-Iy[2]);
  g += R*(-Tx*C*Iz[1]+C*Iz[2]-(1.0/3)*Tx*Iz[3]+(1.0/3)*Iz[4]);

  return g;
}


#ifdef FMT_WEIGHTS_BEFORE_JUN_2021 // could be removed following instructions in FMT_Species::Initialize()
void FMT_Species::generateWeights_before_jun_2021(double hsd, vector<FMT_Weighted_Density> &fmt_weights)
{
  cout << "OK: Generating weights" << endl;
  
  double dx = density_->getDX();
  double dy = density_->getDY();
  double dz = density_->getDZ();

  double dV = dx*dy*dz;

  
  double dS = dx*dx;

  double hsr = hsd/2; // the hard-sphere radius
  
  // This saves having to a code a useless special case
  if(hsd < min(min(dx,dy),dz))
    throw std::runtime_error("hsd is less than the lattice spacing ... aborting");
  
  int Sx_max = 2+int(hsr/dx);
  int Sy_max = 2+int(hsr/dy);
  int Sz_max = 2+int(hsr/dz);

  long pmax = (Sx_max+1)*(Sy_max+1)*(Sz_max+1);
  
  int I[2] = {-1,1};

  cout << endl;
  cout << myColor::GREEN;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << "/////  Generating weights using analytic formulae" << endl;
  cout << myColor::RESET << endl;

  long counter = 0;

  int numWeights = fmt_weights.size();

  //REGARDING dx,dy,dz: In the following, I allow for different spacings in the different directions. In order to preserve - exactly - the previous results, I essentially
  // scale everything by dx. A more rationale implementation would involve changes as noted below in comments marked with !!!. 

  //!!! do not scale by dx  
  double R = hsr/dx; 
  dy /= dx; 
  dz /= dx;
  dx = 1;
  
  for(int Sx = 0; Sx <= Sx_max; Sx++)
    for(int Sy = 0; Sy <= Sy_max; Sy++)
      for(int Sz = 0; Sz <= Sz_max; Sz++)
	{
	  counter++;
	  if(counter%1000 == 0) {if(counter > 0) cout << '\r'; cout << "\t" << int(double(counter)*100.0/pmax) << "% finished: " << counter << " out of " << pmax; cout.flush();}

	  double R2_min = dx*dx*(Sx-(Sx == 0 ? 0 : 1))*(Sx-(Sx == 0 ? 0 : 1))+dy*dy*(Sy-(Sy == 0 ? 0 : 1))*(Sy-(Sy == 0 ? 0 : 1))+dz*dz*(Sz-(Sz == 0 ? 0 : 1))*(Sz-(Sz == 0 ? 0 : 1));
	  double R2_max = dx*dx*(Sx+1)*(Sx+1)+dy*dy*(Sy+1)*(Sy+1)+dz*dz*(Sz+1)*(Sz+1);
	  
	  double w_eta = 0.0;
	  double w_s   = 0.0;
	  double w_v[3] = {0.0,0.0,0.0};
	  double w_T[3][3] = {{0.0,0.0,0.0},
			      {0.0,0.0,0.0},
			      {0.0,0.0,0.0}};

	  // The weight for point Sx,Sy,Sz has contributions from all adjoining cells
	  // The furthest corner is Sx+1,Sy+1,Sz+1 and if this less than hsr*hsr, the volume weights are 1, and the surface weights are zero
	  // else, the nearest corner is Sx-1,Sy-1,Sz-1 (unless Sx, Sy or Sz = 0) and if this is less than hsr*hsr, then the boundary is between these limits and we must compute
	  // else, all hsd boundary is less than the nearest corner and all weights are zero.

	  //Note: Special cases of I-sums, e.g. when one or more components of S are zero, are handled in the called functions.
	  
	  if(R*R > R2_max) {w_eta = dV;}
	  else if(R*R > R2_min)
	    for(int ax:I)
	      {
		int ix = ax;
		double Tx = dx*(Sx+ix);
		int px = (Tx < 0 ? -1 : 1);
		if(Tx < 0) {ix = 1; Tx = dx;}		   
		for(int ay:I)
		  {
		    int iy = ay;
		    double Ty = dy*(Sy+iy);
		    int py = (Ty < 0 ? -1 : 1);
		    if(Ty < 0) {iy = 1; Ty = dy;}	       		    
		    for(int az:I)
		      {
			int iz = az;
			double Tz = dz*(Sz+iz);
			int pz = (Tz < 0 ? -1 : 1);
			if(Tz < 0) {iz = 1; Tz = dz;}
  
			int vx[2] = {0,ix};
			int vy[2] = {0,iy};
			int vz[2] = {0,iz};

			double j = 0.0;

			for(int jx = 0; jx < 2; jx++)
			  for(int jy = 0; jy < 2; jy++)
			    for(int jz = 0; jz < 2; jz++)
			      {
				double Vx = Tx - dx*vx[jx];
				double Vy = Ty - dy*vy[jy];
				double Vz = Tz - dz*vz[jz];
			    
				int sgn = 1;
				if(jx == 1) sgn *= -1;
				if(jy == 1) sgn *= -1;
				if(jz == 1) sgn *= -1;
			    
				// test up to machine precision
				if((R*R - (Vx*Vx+Vy*Vy+Vz*Vz)) < std::nextafter(0.d,1.d)) continue;

				//!!! Replace dV and dS in the following by (1/dV) in all cases
				
				w_eta += dV*sgn*(G_eta(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_eta(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				if(numWeights > 1)
				  {
				    w_s   += dS*sgn*(G_s(R,Vx,Vy,Vz,Tx,Ty,Tz)   - G_s(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));

				    w_v[0] += dS*px*(1.0/R)*sgn*(G_vx(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_vx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_v[1] += dS*py*(1.0/R)*sgn*(G_vx(R,Vy,Vx,Vz,Ty,Tx,Tz) - G_vx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz));
				    w_v[2] += dS*pz*(1.0/R)*sgn*(G_vx(R,Vz,Vy,Vx,Tz,Ty,Tx) - G_vx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx));

				    w_T[0][0] += dS*(1.0/(R*R))*sgn*(G_txx(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_txx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_T[1][1] += dS*(1.0/(R*R))*sgn*(G_txx(R,Vy,Vx,Vz,Ty,Tx,Tz) - G_txx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz));
				    w_T[2][2] += dS*(1.0/(R*R))*sgn*(G_txx(R,Vz,Vy,Vx,Tz,Ty,Tx) - G_txx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx));			    

				    w_T[0][1] += dS*px*py*(1.0/(R*R))*sgn*(G_txy(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_T[0][2] += dS*px*pz*(1.0/(R*R))*sgn*(G_txy(R,Vx,Vz,Vy,Tx,Tz,Ty) - G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vz,Vy,Tx,Tz,Ty));
				    w_T[1][2] += dS*py*pz*(1.0/(R*R))*sgn*(G_txy(R,Vy,Vz,Vx,Ty,Tz,Tx) - G_txy(R,sqrt(R*R-Vz*Vz-Vx*Vx),Vz,Vx,Ty,Tz,Tx));
				  }
			      }
		      }
		  }
	      }	
	  // Add in for all octants of the sphere: take account of parity of vector and tensor quantities
	  for(int ix = 0; ix < (Sx == 0 ? 1 : 2); ix++)
	    for(int iy = 0; iy < (Sy == 0 ? 1 : 2); iy++)
	      for(int iz = 0; iz < (Sz == 0 ? 1 : 2); iz++)
		{		  
		  long pos = density_->get_PBC_Pos((1-2*ix)*Sx,(1-2*iy)*Sy,(1-2*iz)*Sz);	  
		  fmt_weights[EI()].addToWeight(pos,w_eta);
		  if(isnan(fmt_weights[EI()].getWeight(pos)))
		    {
		      cout << ix << " " << iy << " " << iz << " " << Sx << " " << Sy << " " << Sz << endl;
		      throw std::runtime_error("Found NAN");
		    }
		  if(numWeights > 1)
		    {
		      fmt_weights[SI()].addToWeight(pos,w_s);
		      for(int iv = 0;iv < 3;iv++)
			{
			  fmt_weights[VI(iv)].addToWeight(pos,(iv == 0 ? (1-2*ix) : (iv == 1 ? (1-2*iy) : (1-2*iz)))*w_v[iv]);
			  for(int it=iv;it<3;it++)
			    fmt_weights[TI(iv,it)].addToWeight(pos,(iv == 0 ? (1-2*ix) : (iv == 1 ? (1-2*iy) : (1-2*iz)))*(it == 0 ? (1-2*ix) : (it == 1 ? (1-2*iy) : (1-2*iz)))*w_T[iv][it]);
			}
		    }
		}		  
	}

  cout << endl;
  cout << myColor::GREEN;
  cout << "/////  Finished.  " << endl;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << myColor::RESET << endl;
}
#endif //FMT_WEIGHTS_BEFORE_JUN_2021

#ifdef FMT_WEIGHTS_BEFORE_JAN_2023 // could be removed following instructions in FMT_Species::Initialize()
void FMT_Species::generateWeights_before_jan_2023(double hsd, vector<FMT_Weighted_Density> &fmt_weights)
{
  double dx = density_->getDX();
  double dy = density_->getDY();
  double dz = density_->getDZ();

  double dV = dx*dy*dz;

  double R = hsd/2; // the hard-sphere radius
  
  // This saves having to a code a useless special case
  if(hsd < min(min(dx,dy),dz))
    throw std::runtime_error("hsd is less than the lattice spacing ... aborting");
  
  int Sx_max = 2+int(R/dx);
  int Sy_max = 2+int(R/dy);
  int Sz_max = 2+int(R/dz);

  long pmax = (Sx_max+1)*(Sy_max+1)*(Sz_max+1);
  
  int I[2] = {-1,1};

  if(verbose_) cout << endl;
  if(verbose_) cout << myColor::YELLOW;
  if(verbose_) cout << "/////  Generating FMT weights using analytic formulae" << endl;
  if(verbose_) cout << myColor::RESET;

  long counter = 0;

  int numWeights = fmt_weights.size();
  
  for(int Sx = 0; Sx <= Sx_max; Sx++)
    for(int Sy = 0; Sy <= Sy_max; Sy++)
      for(int Sz = 0; Sz <= Sz_max; Sz++)
	{
	  counter++;
	  if(counter%1000 == 0) {if(counter > 0) if(verbose_) cout << '\r'; if(verbose_) cout << "\t" << int(double(counter)*100.0/pmax) << "% finished: " << counter << " out of " << pmax; if(verbose_) cout.flush();}

	  double R2_min = dx*dx*(Sx-(Sx == 0 ? 0 : 1))*(Sx-(Sx == 0 ? 0 : 1))+dy*dy*(Sy-(Sy == 0 ? 0 : 1))*(Sy-(Sy == 0 ? 0 : 1))+dz*dz*(Sz-(Sz == 0 ? 0 : 1))*(Sz-(Sz == 0 ? 0 : 1));
	  double R2_max = dx*dx*(Sx+1)*(Sx+1)+dy*dy*(Sy+1)*(Sy+1)+dz*dz*(Sz+1)*(Sz+1);
	  
	  double w_eta = 0.0;
	  double w_s   = 0.0;
	  double w_v[3] = {0.0,0.0,0.0};
	  double w_T[3][3] = {{0.0,0.0,0.0},
			      {0.0,0.0,0.0},
			      {0.0,0.0,0.0}};

	  // The weight for point Sx,Sy,Sz has contributions from all adjoining cells
	  // The furthest corner is Sx+1,Sy+1,Sz+1 and if this less than hsr*hsr, the volume weights are 1, and the surface weights are zero
	  // else, the nearest corner is Sx-1,Sy-1,Sz-1 (unless Sx, Sy or Sz = 0) and if this is less than hsr*hsr, then the boundary is between these limits and we must compute
	  // else, all hsd boundary is less than the nearest corner and all weights are zero.

	  //Note: Special cases of I-sums, e.g. when one or more components of S are zero, are handled in the called functions.
	  
	  if(R*R > R2_max) {w_eta = dV;}
	  else if(R*R > R2_min)
	    for(int ax:I)
	      {
		int ix = ax;
		double Tx = dx*(Sx+ix);
		int px = (Tx < 0 ? -1 : 1);
		if(Tx < 0) {ix = 1; Tx = dx;}		   
		for(int ay:I)
		  {
		    int iy = ay;
		    double Ty = dy*(Sy+iy);
		    int py = (Ty < 0 ? -1 : 1);
		    if(Ty < 0) {iy = 1; Ty = dy;}	       		    
		    for(int az:I)
		      {
			int iz = az;
			double Tz = dz*(Sz+iz);
			int pz = (Tz < 0 ? -1 : 1);
			if(Tz < 0) {iz = 1; Tz = dz;}
  
			int vx[2] = {0,ix};
			int vy[2] = {0,iy};
			int vz[2] = {0,iz};

			double j = 0.0;

			for(int jx = 0; jx < 2; jx++)
			  for(int jy = 0; jy < 2; jy++)
			    for(int jz = 0; jz < 2; jz++)
			      {
				double Vx = Tx - dx*vx[jx];
				double Vy = Ty - dy*vy[jy];
				double Vz = Tz - dz*vz[jz];
			    
				int sgn = 1;
				if(jx == 1) sgn *= -1;
				if(jy == 1) sgn *= -1;
				if(jz == 1) sgn *= -1;
			    
				// test up to machine precision
				if((R*R - (Vx*Vx+Vy*Vy+Vz*Vz)) < std::nextafter(0.d,1.d)) continue;
				
				w_eta += (1/dV)*sgn*(G_eta(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_eta(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				if(numWeights > 1)
				  {
				    w_s   += (1/dV)*sgn*(G_s(R,Vx,Vy,Vz,Tx,Ty,Tz)   - G_s(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));

				    w_v[0] += (1/dV)*px*(1.0/R)*sgn*(G_vx(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_vx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_v[1] += (1/dV)*py*(1.0/R)*sgn*(G_vx(R,Vy,Vx,Vz,Ty,Tx,Tz) - G_vx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz));
				    w_v[2] += (1/dV)*pz*(1.0/R)*sgn*(G_vx(R,Vz,Vy,Vx,Tz,Ty,Tx) - G_vx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx));

				    w_T[0][0] += (1/dV)*(1.0/(R*R))*sgn*(G_txx(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_txx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_T[1][1] += (1/dV)*(1.0/(R*R))*sgn*(G_txx(R,Vy,Vx,Vz,Ty,Tx,Tz) - G_txx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz));
				    w_T[2][2] += (1/dV)*(1.0/(R*R))*sgn*(G_txx(R,Vz,Vy,Vx,Tz,Ty,Tx) - G_txx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx));			    

				    w_T[0][1] += (1/dV)*px*py*(1.0/(R*R))*sgn*(G_txy(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_T[0][2] += (1/dV)*px*pz*(1.0/(R*R))*sgn*(G_txy(R,Vx,Vz,Vy,Tx,Tz,Ty) - G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vz,Vy,Tx,Tz,Ty));
				    w_T[1][2] += (1/dV)*py*pz*(1.0/(R*R))*sgn*(G_txy(R,Vy,Vz,Vx,Ty,Tz,Tx) - G_txy(R,sqrt(R*R-Vz*Vz-Vx*Vx),Vz,Vx,Ty,Tz,Tx));
				  }
			      }
		      }
		  }
	      }
	  
	  // Add in for all octants of the sphere: take account of parity of vector and tensor quantities
	  for(int ix = 0; ix < (Sx == 0 ? 1 : 2); ix++)
	    for(int iy = 0; iy < (Sy == 0 ? 1 : 2); iy++)
	      for(int iz = 0; iz < (Sz == 0 ? 1 : 2); iz++)
		{		  
		  long pos = density_->get_PBC_Pos((1-2*ix)*Sx,(1-2*iy)*Sy,(1-2*iz)*Sz);	  
		  fmt_weights[EI()].addToWeight(pos,w_eta);
		  if(std::isnan(fmt_weights[EI()].getWeight(pos)))
		    {
		      if(verbose_) cout << ix << " " << iy << " " << iz << " " << Sx << " " << Sy << " " << Sz << endl;
		      throw std::runtime_error("Found NAN");
		    }
		  if(numWeights > 1)
		    {
		      fmt_weights[SI()].addToWeight(pos,w_s);
		      for(int iv = 0;iv < 3;iv++)
			{
			  fmt_weights[VI(iv)].addToWeight(pos,(iv == 0 ? (1-2*ix) : (iv == 1 ? (1-2*iy) : (1-2*iz)))*w_v[iv]);
			  for(int it=iv;it<3;it++)
			    fmt_weights[TI(iv,it)].addToWeight(pos,(iv == 0 ? (1-2*ix) : (iv == 1 ? (1-2*iy) : (1-2*iz)))*(it == 0 ? (1-2*ix) : (it == 1 ? (1-2*iy) : (1-2*iz)))*w_T[iv][it]);
			}
		    }
		}		  
	}

  if(verbose_) cout << '\r'; if(verbose_) cout << ""; if(verbose_) cout.flush();
}
#endif //FMT_WEIGHTS_BEFORE_JAN_2023

//New FMT weight calculation (after january 2023)
void FMT_Species::generateWeights(double hsd, vector<FMT_Weighted_Density> &fmt_weights, double scale, bool subtract_regulator_contribution)
{
  #ifndef DEBUG_FMT_WEIGHTS
    #if   defined FMT_WEIGHTS_BEFORE_JUN_2021
      throw runtime_error("FMT_Species::generateWeights should not be called in this context");
    #elif defined FMT_WEIGHTS_BEFORE_JAN_2023
      throw runtime_error("FMT_Species::generateWeights should not be called in this context");
    #endif
  #endif
  
  double dx = density_->getDX();
  double dy = density_->getDY();
  double dz = density_->getDZ();

  double dV = dx*dy*dz;

  double R = hsd/2; // the hard-sphere radius
  
  // This saves having to a code a useless special case
  if(hsd < min(min(dx,dy),dz))
    throw std::runtime_error("hsd is less than the lattice spacing ... aborting");
  
  int Sx_max = 2+int(R/dx);
  int Sy_max = 2+int(R/dy);
  int Sz_max = 2+int(R/dz);

  long pmax = (Sx_max+1)*(Sy_max+1)*(Sz_max+1);
  
  int I[2] = {-1,1};

  if(verbose_) cout << endl;
  if(verbose_) cout << myColor::YELLOW;
  if(verbose_) cout << "/////  Generating FMT weights using analytic formulae" << endl;
  if(verbose_) cout << myColor::RESET;

  long counter = 0;
  long num_contributions = 0;

  int numWeights = fmt_weights.size();
  
  ////////////////////////////  
  // Scaling here gives the same results
  double _dx = dx; double _dy = dy; double _dz = dz;
  double  _R =  R; double _dV = dV;
  /*
  dx /= scale; 
  dy /= scale; 
  dz /= scale;
  
  R  /= scale;
  dV /= pow(scale,3);
  */
  ////////////////////////////
  
  for(int Sx = 0; Sx <= Sx_max; Sx++)
    for(int Sy = 0; Sy <= Sy_max; Sy++)
      for(int Sz = 0; Sz <= Sz_max; Sz++)
	{
	  ////////////////////////////
	  // Scaling here gives the same results
	  dx = _dx; dy = _dy; dz = _dz;
	  R  =  _R; dV = _dV;
	  /*
	  dx /= scale;
	  dy /= scale;
	  dz /= scale;
	  
	  R  /= scale;
	  dV /= pow(scale,3);
	  */
	  ////////////////////////////
	  
	  counter++;
	  if(counter%1000 == 0) {if(counter > 0) if(verbose_) cout << '\r'; if(verbose_) cout << "\t" << int(double(counter)*100.0/pmax) << "% finished: " << counter << " out of " << pmax; if(verbose_) cout.flush();}

	  double R2_min = dx*dx*(Sx-(Sx == 0 ? 0 : 1))*(Sx-(Sx == 0 ? 0 : 1))+dy*dy*(Sy-(Sy == 0 ? 0 : 1))*(Sy-(Sy == 0 ? 0 : 1))+dz*dz*(Sz-(Sz == 0 ? 0 : 1))*(Sz-(Sz == 0 ? 0 : 1));
	  double R2_max = dx*dx*(Sx+1)*(Sx+1)+dy*dy*(Sy+1)*(Sy+1)+dz*dz*(Sz+1)*(Sz+1);
	  /*
	  double w_eta = 0.0;
	  double w_s   = 0.0;
	  double w_v[3] = {0.0,0.0,0.0};
	  double w_T[3][3] = {{0.0,0.0,0.0},
			      {0.0,0.0,0.0},
			      {0.0,0.0,0.0}};
    */
    Summation w_eta, w_s, w_vx, w_vy, w_vz, w_txx, w_tyy, w_tzz, w_txy, w_txz, w_tyz;
    
	  // The weight for point Sx,Sy,Sz has contributions from all adjoining cells
	  // The furthest corner is Sx+1,Sy+1,Sz+1 and if this less than hsr*hsr, the volume weights are 1, and the surface weights are zero
	  // else, the nearest corner is Sx-1,Sy-1,Sz-1 (unless Sx, Sy or Sz = 0) and if this is less than hsr*hsr, then the boundary is between these limits and we must compute
	  // else, all hsd boundary is less than the nearest corner and all weights are zero.

	  //Note: Special cases of I-sums, e.g. when one or more components of S are zero, are handled in the called functions.
	  
	  if(R*R > R2_max) {w_eta += dV;}
	  //if(R*R > R2_max) {w_eta += dV*pow(scale,3);} // use this line if you scale above
	  else if(R*R > R2_min)
	  {
	    ////////////////////////////
	    // Scaling here gives the same results
	    
	    dx /= scale;
	    dy /= scale;
	    dz /= scale;
	    
	    R  /= scale;
	    dV /= pow(scale,3);
	    
	    ////////////////////////////
      
	    for(int ax:I)
	      {
		int ix = ax;
		double Tx = dx*(Sx+ix);
		int px = (Tx < -dx/10 ? -1 : 1);
		if(Tx < -dx/10) {ix = 1; Tx = dx;}		   
		for(int ay:I)
		  {
		    int iy = ay;
		    double Ty = dy*(Sy+iy);
		    int py = (Ty < -dy/10 ? -1 : 1);
		    if(Ty < -dy/10) {iy = 1; Ty = dy;}	       		    
		    for(int az:I)
		      {
			int iz = az;
			double Tz = dz*(Sz+iz);
			int pz = (Tz < -dz/10 ? -1 : 1);
			if(Tz < -dz/10) {iz = 1; Tz = dz;}
  
			int vx[2] = {0,ix};
			int vy[2] = {0,iy};
			int vz[2] = {0,iz};

			for(int jx = 0; jx < 2; jx++)
			  for(int jy = 0; jy < 2; jy++)
			    for(int jz = 0; jz < 2; jz++)
			      {
				double Vx = Tx - dx*vx[jx];
				double Vy = Ty - dy*vy[jy];
				double Vz = Tz - dz*vz[jz];
			    
				int sgn = 1;
				if(jx == 1) sgn *= -1;
				if(jy == 1) sgn *= -1;
				if(jz == 1) sgn *= -1;
			    
				// test up to machine precision
				if((R*R - (Vx*Vx+Vy*Vy+Vz*Vz)) < std::nextafter(0.d,1.d)) continue;
				
				////////////////////////
				/*if (fabs(dx-dy)>0 || fabs(dx-dz)>0)
				{
				  cout << endl; cout << scientific << setprecision(2);
				  cout << "dx-dy = " << dx-dy << endl;
				  cout << "dx-dz = " << dx-dz << endl;
				  throw runtime_error("dx!=dy!=dz");
				}*/
				////////////////////////
				
				////////////////////////
				// Trying to identify where inaccuracies come from
				/*if (fabs(Vx/R)<SMALL_VALUE_FOR_J_INTEGRAL ||    // problems with vy,vz and tyy,tzz,tyz
				    fabs(Vy/R)<SMALL_VALUE_FOR_J_INTEGRAL ||    // problems with all but vy and tyy,tyz
				    fabs(Vz/R)<SMALL_VALUE_FOR_J_INTEGRAL       // problems with all but vz and tzz
				//    ) continue;
				    ); else continue;*/
				////////////////////////
				
				num_contributions++;
				
				// Scaling here gives the DIFFERENT results if you are using the old regulator
				/*
				double _Vx = Vx; double _Vy = Vy; double _Vz = Vz;
				double _Tx = Tx; double _Ty = Ty; double _Tz = Tz;
				
				R  /= scale; dV /= pow(scale,3);
				Vx /= scale; Vy /= scale; Vz /= scale;
				Tx /= scale; Ty /= scale; Tz /= scale;
				*/
				/*
				w_eta += pow(scale,3)*(1/dV)*sgn*(G_eta(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_eta(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				if(numWeights > 1)
				  {
				    w_s   += pow(scale,2)*(1/dV)*sgn*(G_s(R,Vx,Vy,Vz,Tx,Ty,Tz)   - G_s(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));

				    w_v[0] += pow(scale,2)*(1/dV)*px*(1.0/R)*sgn*(G_vx(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_vx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_v[1] += pow(scale,2)*(1/dV)*py*(1.0/R)*sgn*(G_vx(R,Vy,Vx,Vz,Ty,Tx,Tz) - G_vx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz));
				    w_v[2] += pow(scale,2)*(1/dV)*pz*(1.0/R)*sgn*(G_vx(R,Vz,Vy,Vx,Tz,Ty,Tx) - G_vx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx));

				    w_T[0][0] += pow(scale,2)*(1/dV)*(1.0/(R*R))*sgn*(G_txx(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_txx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_T[1][1] += pow(scale,2)*(1/dV)*(1.0/(R*R))*sgn*(G_txx(R,Vy,Vx,Vz,Ty,Tx,Tz) - G_txx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz));
				    w_T[2][2] += pow(scale,2)*(1/dV)*(1.0/(R*R))*sgn*(G_txx(R,Vz,Vy,Vx,Tz,Ty,Tx) - G_txx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx));			    

				    w_T[0][1] += pow(scale,2)*(1/dV)*px*py*(1.0/(R*R))*sgn*(G_txy(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_T[0][2] += pow(scale,2)*(1/dV)*px*pz*(1.0/(R*R))*sgn*(G_txy(R,Vx,Vz,Vy,Tx,Tz,Ty) - G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vz,Vy,Tx,Tz,Ty));
				    w_T[1][2] += pow(scale,2)*(1/dV)*py*pz*(1.0/(R*R))*sgn*(G_txy(R,Vy,Vz,Vx,Ty,Tz,Tx) - G_txy(R,sqrt(R*R-Vz*Vz-Vx*Vx),Vz,Vx,Ty,Tz,Tx));
				  }
				*/
				
				if (!subtract_regulator_contribution)
				{
				w_eta += pow(scale,3)*(1/dV)*sgn*G_eta(R,Vx,Vy,Vz,Tx,Ty,Tz);//- G_eta(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				w_eta -= pow(scale,3)*(1/dV)*sgn*G_eta(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz);
				if(numWeights > 1)
				  {
				    w_s   += pow(scale,2)*(1/dV)*sgn*G_s(R,Vx,Vy,Vz,Tx,Ty,Tz);//- G_s(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_s   -= pow(scale,2)*(1/dV)*sgn*G_s(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz);

				    w_vx  += pow(scale,2)*(1/dV)*px*(1.0/R)*sgn*G_vx(R,Vx,Vy,Vz,Tx,Ty,Tz);// - G_vx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_vx  -= pow(scale,2)*(1/dV)*px*(1.0/R)*sgn*G_vx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz);
				    w_vy  += pow(scale,2)*(1/dV)*py*(1.0/R)*sgn*G_vx(R,Vy,Vx,Vz,Ty,Tx,Tz);// - G_vx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz));
				    w_vy  -= pow(scale,2)*(1/dV)*py*(1.0/R)*sgn*G_vx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz);
				    w_vz  += pow(scale,2)*(1/dV)*pz*(1.0/R)*sgn*G_vx(R,Vz,Vy,Vx,Tz,Ty,Tx);// - G_vx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx));
				    w_vz  -= pow(scale,2)*(1/dV)*pz*(1.0/R)*sgn*G_vx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx);

				    w_txx += pow(scale,2)*(1/dV)*(1.0/(R*R))*sgn*G_txx(R,Vx,Vy,Vz,Tx,Ty,Tz);// - G_txx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_txx -= pow(scale,2)*(1/dV)*(1.0/(R*R))*sgn*G_txx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz);
				    w_tyy += pow(scale,2)*(1/dV)*(1.0/(R*R))*sgn*G_txx(R,Vy,Vx,Vz,Ty,Tx,Tz);// - G_txx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz));
				    w_tyy -= pow(scale,2)*(1/dV)*(1.0/(R*R))*sgn*G_txx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz);
				    w_tzz += pow(scale,2)*(1/dV)*(1.0/(R*R))*sgn*G_txx(R,Vz,Vy,Vx,Tz,Ty,Tx);// - G_txx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx));			    
				    w_tzz -= pow(scale,2)*(1/dV)*(1.0/(R*R))*sgn*G_txx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx);

				    w_txy += pow(scale,2)*(1/dV)*px*py*(1.0/(R*R))*sgn*G_txy(R,Vx,Vy,Vz,Tx,Ty,Tz);// - G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_txy -= pow(scale,2)*(1/dV)*px*py*(1.0/(R*R))*sgn*G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz);
				    w_txz += pow(scale,2)*(1/dV)*px*pz*(1.0/(R*R))*sgn*G_txy(R,Vx,Vz,Vy,Tx,Tz,Ty);// - G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vz,Vy,Tx,Tz,Ty));
				    w_txz -= pow(scale,2)*(1/dV)*px*pz*(1.0/(R*R))*sgn*G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vz,Vy,Tx,Tz,Ty);
				    w_tyz += pow(scale,2)*(1/dV)*py*pz*(1.0/(R*R))*sgn*G_txy(R,Vy,Vz,Vx,Ty,Tz,Tx);// - G_txy(R,sqrt(R*R-Vz*Vz-Vx*Vx),Vz,Vx,Ty,Tz,Tx));
				    w_tyz -= pow(scale,2)*(1/dV)*py*pz*(1.0/(R*R))*sgn*G_txy(R,sqrt(R*R-Vz*Vz-Vx*Vx),Vz,Vx,Ty,Tz,Tx);
				  }
				}
				else
				{
				// Here we compute twice the weight increment. Due to the small regulator
				// in G-functions (SMALL_VALUE_FOR_WORKAROUND_NAN_IN_G_FUNCTION), the 
				// expected result dW is off by a small amount dw. We can subtract this
				// perturbation by multiplying our result dW+dw by a factor of 2 and then
				// subtract the result dW+2*dw of a second calculation with a regulator
				// twice as large. This effectively removes the perturbation introduced
				// by the regulator (at first order).
				
				double prefac = 2;
				w_eta += prefac*pow(scale,3)*(1/dV)*sgn*G_eta(R,Vx,Vy,Vz,Tx,Ty,Tz);//- G_eta(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				w_eta -= prefac*pow(scale,3)*(1/dV)*sgn*G_eta(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz);
				if(numWeights > 1)
				  {
				    w_s   += prefac*pow(scale,2)*(1/dV)*sgn*G_s(R,Vx,Vy,Vz,Tx,Ty,Tz);//- G_s(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_s   -= prefac*pow(scale,2)*(1/dV)*sgn*G_s(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz);

				    w_vx  += prefac*pow(scale,2)*(1/dV)*px*(1.0/R)*sgn*G_vx(R,Vx,Vy,Vz,Tx,Ty,Tz);// - G_vx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_vx  -= prefac*pow(scale,2)*(1/dV)*px*(1.0/R)*sgn*G_vx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz);
				    w_vy  += prefac*pow(scale,2)*(1/dV)*py*(1.0/R)*sgn*G_vx(R,Vy,Vx,Vz,Ty,Tx,Tz);// - G_vx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz));
				    w_vy  -= prefac*pow(scale,2)*(1/dV)*py*(1.0/R)*sgn*G_vx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz);
				    w_vz  += prefac*pow(scale,2)*(1/dV)*pz*(1.0/R)*sgn*G_vx(R,Vz,Vy,Vx,Tz,Ty,Tx);// - G_vx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx));
				    w_vz  -= prefac*pow(scale,2)*(1/dV)*pz*(1.0/R)*sgn*G_vx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx);

				    w_txx += prefac*pow(scale,2)*(1/dV)*(1.0/(R*R))*sgn*G_txx(R,Vx,Vy,Vz,Tx,Ty,Tz);// - G_txx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_txx -= prefac*pow(scale,2)*(1/dV)*(1.0/(R*R))*sgn*G_txx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz);
				    w_tyy += prefac*pow(scale,2)*(1/dV)*(1.0/(R*R))*sgn*G_txx(R,Vy,Vx,Vz,Ty,Tx,Tz);// - G_txx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz));
				    w_tyy -= prefac*pow(scale,2)*(1/dV)*(1.0/(R*R))*sgn*G_txx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz);
				    w_tzz += prefac*pow(scale,2)*(1/dV)*(1.0/(R*R))*sgn*G_txx(R,Vz,Vy,Vx,Tz,Ty,Tx);// - G_txx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx));			    
				    w_tzz -= prefac*pow(scale,2)*(1/dV)*(1.0/(R*R))*sgn*G_txx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx);

				    w_txy += prefac*pow(scale,2)*(1/dV)*px*py*(1.0/(R*R))*sgn*G_txy(R,Vx,Vy,Vz,Tx,Ty,Tz);// - G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_txy -= prefac*pow(scale,2)*(1/dV)*px*py*(1.0/(R*R))*sgn*G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz);
				    w_txz += prefac*pow(scale,2)*(1/dV)*px*pz*(1.0/(R*R))*sgn*G_txy(R,Vx,Vz,Vy,Tx,Tz,Ty);// - G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vz,Vy,Tx,Tz,Ty));
				    w_txz -= prefac*pow(scale,2)*(1/dV)*px*pz*(1.0/(R*R))*sgn*G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vz,Vy,Tx,Tz,Ty);
				    w_tyz += prefac*pow(scale,2)*(1/dV)*py*pz*(1.0/(R*R))*sgn*G_txy(R,Vy,Vz,Vx,Ty,Tz,Tx);// - G_txy(R,sqrt(R*R-Vz*Vz-Vx*Vx),Vz,Vx,Ty,Tz,Tx));
				    w_tyz -= prefac*pow(scale,2)*(1/dV)*py*pz*(1.0/(R*R))*sgn*G_txy(R,sqrt(R*R-Vz*Vz-Vx*Vx),Vz,Vx,Ty,Tz,Tx);
				  }
				
				SMALL_VALUE_FOR_WORKAROUND_NAN_IN_G_FUNCTION *= 2;
				prefac = -1;
				w_eta += prefac*pow(scale,3)*(1/dV)*sgn*G_eta(R,Vx,Vy,Vz,Tx,Ty,Tz);//- G_eta(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				w_eta -= prefac*pow(scale,3)*(1/dV)*sgn*G_eta(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz);
				if(numWeights > 1)
				  {
				    w_s   += prefac*pow(scale,2)*(1/dV)*sgn*G_s(R,Vx,Vy,Vz,Tx,Ty,Tz);//- G_s(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_s   -= prefac*pow(scale,2)*(1/dV)*sgn*G_s(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz);

				    w_vx  += prefac*pow(scale,2)*(1/dV)*px*(1.0/R)*sgn*G_vx(R,Vx,Vy,Vz,Tx,Ty,Tz);// - G_vx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_vx  -= prefac*pow(scale,2)*(1/dV)*px*(1.0/R)*sgn*G_vx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz);
				    w_vy  += prefac*pow(scale,2)*(1/dV)*py*(1.0/R)*sgn*G_vx(R,Vy,Vx,Vz,Ty,Tx,Tz);// - G_vx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz));
				    w_vy  -= prefac*pow(scale,2)*(1/dV)*py*(1.0/R)*sgn*G_vx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz);
				    w_vz  += prefac*pow(scale,2)*(1/dV)*pz*(1.0/R)*sgn*G_vx(R,Vz,Vy,Vx,Tz,Ty,Tx);// - G_vx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx));
				    w_vz  -= prefac*pow(scale,2)*(1/dV)*pz*(1.0/R)*sgn*G_vx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx);

				    w_txx += prefac*pow(scale,2)*(1/dV)*(1.0/(R*R))*sgn*G_txx(R,Vx,Vy,Vz,Tx,Ty,Tz);// - G_txx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_txx -= prefac*pow(scale,2)*(1/dV)*(1.0/(R*R))*sgn*G_txx(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz);
				    w_tyy += prefac*pow(scale,2)*(1/dV)*(1.0/(R*R))*sgn*G_txx(R,Vy,Vx,Vz,Ty,Tx,Tz);// - G_txx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz));
				    w_tyy -= prefac*pow(scale,2)*(1/dV)*(1.0/(R*R))*sgn*G_txx(R,sqrt(R*R-Vx*Vx-Vz*Vz),Vx,Vz,Ty,Tx,Tz);
				    w_tzz += prefac*pow(scale,2)*(1/dV)*(1.0/(R*R))*sgn*G_txx(R,Vz,Vy,Vx,Tz,Ty,Tx);// - G_txx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx));			    
				    w_tzz -= prefac*pow(scale,2)*(1/dV)*(1.0/(R*R))*sgn*G_txx(R,sqrt(R*R-Vy*Vy-Vx*Vx),Vy,Vx,Tz,Ty,Tx);

				    w_txy += prefac*pow(scale,2)*(1/dV)*px*py*(1.0/(R*R))*sgn*G_txy(R,Vx,Vy,Vz,Tx,Ty,Tz);// - G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
				    w_txy -= prefac*pow(scale,2)*(1/dV)*px*py*(1.0/(R*R))*sgn*G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz);
				    w_txz += prefac*pow(scale,2)*(1/dV)*px*pz*(1.0/(R*R))*sgn*G_txy(R,Vx,Vz,Vy,Tx,Tz,Ty);// - G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vz,Vy,Tx,Tz,Ty));
				    w_txz -= prefac*pow(scale,2)*(1/dV)*px*pz*(1.0/(R*R))*sgn*G_txy(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vz,Vy,Tx,Tz,Ty);
				    w_tyz += prefac*pow(scale,2)*(1/dV)*py*pz*(1.0/(R*R))*sgn*G_txy(R,Vy,Vz,Vx,Ty,Tz,Tx);// - G_txy(R,sqrt(R*R-Vz*Vz-Vx*Vx),Vz,Vx,Ty,Tz,Tx));
				    w_tyz -= prefac*pow(scale,2)*(1/dV)*py*pz*(1.0/(R*R))*sgn*G_txy(R,sqrt(R*R-Vz*Vz-Vx*Vx),Vz,Vx,Ty,Tz,Tx);
				  }
				SMALL_VALUE_FOR_WORKAROUND_NAN_IN_G_FUNCTION /= 2;
				}
				/*
				R  = _R;  dV = _dV;
				Vx = _Vx; Vy = _Vy; Vz = _Vz;
				Tx = _Tx; Ty = _Ty; Tz = _Tz;
				*/
			      }
		      }
		  }
	      }
	      }
	  
	  double w_v[3] = {w_vx.sum(),w_vy.sum(),w_vz.sum()};
	  double w_T[3][3] = {{w_txx.sum(),w_txy.sum(),w_txz.sum()},
	                      {        0.0,w_tyy.sum(),w_tyz.sum()},
		                    {        0.0,        0.0,w_tzz.sum()}};
	  
	  // Add in for all octants of the sphere: take account of parity of vector and tensor quantities
	  for(int ix = 0; ix < (Sx == 0 ? 1 : 2); ix++)
	    for(int iy = 0; iy < (Sy == 0 ? 1 : 2); iy++)
	      for(int iz = 0; iz < (Sz == 0 ? 1 : 2); iz++)
		{		  
		  long pos = density_->get_PBC_Pos((1-2*ix)*Sx,(1-2*iy)*Sy,(1-2*iz)*Sz);	  
		  fmt_weights[EI()].addToWeight(pos,w_eta.sum());
		  //fmt_weights[EI()].addToWeight(pos,w_eta);
		  if(std::isnan(fmt_weights[EI()].getWeight(pos)))
		    {
		      if(verbose_) cout << ix << " " << iy << " " << iz << " " << Sx << " " << Sy << " " << Sz << endl;
		      throw std::runtime_error("Found NAN");
		    }
		  if(numWeights > 1)
		    {
		      fmt_weights[SI()].addToWeight(pos,w_s.sum());
		      //fmt_weights[SI()].addToWeight(pos,w_s);
		      for(int iv = 0;iv < 3;iv++)
			{
			  fmt_weights[VI(iv)].addToWeight(pos,(iv == 0 ? (1-2*ix) : (iv == 1 ? (1-2*iy) : (1-2*iz)))*w_v[iv]);
			  for(int it=iv;it<3;it++)
			    fmt_weights[TI(iv,it)].addToWeight(pos,(iv == 0 ? (1-2*ix) : (iv == 1 ? (1-2*iy) : (1-2*iz)))*(it == 0 ? (1-2*ix) : (it == 1 ? (1-2*iy) : (1-2*iz)))*w_T[iv][it]);
			}
		    }
		}		  
	}

  if(verbose_)
    cout << '\r'; if(verbose_) cout << ""; if(verbose_) cout.flush();
}


FMT_Species_EOS::FMT_Species_EOS(double D_EOS, Density& density, double hsd, double mu, int seq)
  : FMT_Species(density,hsd,mu,seq), eos_weighted_density_(1), D_EOS_(D_EOS)
										    
{
  long Nx = density_->Nx();
  long Ny = density_->Ny();
  long Nz = density_->Nz();

  eos_weighted_density_[0].initialize(Nx, Ny, Nz);
  generateWeights(D_EOS*hsd, eos_weighted_density_);
  eos_weighted_density_[0].transformWeights();  
}
/*  //NOW OUTDATED
// This is an exact copy of FMT_Species::generateWeights. It is necessary because we want to do
// the same but with a different hsd. Obviously, these could (and probably SHOULD) be combined
// with some additional logic to only do the work required. For now, this is a quick solution
// with the only drawback being that essentially the same code must be maintained in two different
// places (e.g. if any errors are found). Since the code is fairly compact, I do not see this as a major problem.

void FMT_Species_EOS::generate_additional_Weight()
{
  cout << "Generating additional FMT eos weight" << endl;
  
  double dx = density_->getDX();
  double dy = density_->getDY();
  double dz = density_->getDZ();

  double dV = dx*dy*dz;
  double dS = dx*dy;

  double hsr = hsd_; // twide the hard-sphere radius
  
  // This saves having to a code a useless special case
  if(2*hsd_ < dx)
    throw std::runtime_error("2*hsd is less than the lattice spacing ... aborting");
  
  int Sx_max = 2+int(hsr/dx);
  int Sy_max = 2+int(hsr/dy);
  int Sz_max = 2+int(hsr/dz);

  long pmax = (Sx_max+1)*(Sy_max+1)*(Sz_max+1);
  
  int I[2] = {-1,1};

  cout << endl;
  cout << myColor::YELLOW;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << "/////  Generating FMT eos weight using analytic formulae" << endl;
  cout << myColor::RESET << endl;

  long counter = 0;
  
  for(int Sx = 0; Sx <= Sx_max; Sx++)
    for(int Sy = 0; Sy <= Sy_max; Sy++)
      for(int Sz = 0; Sz <= Sz_max; Sz++)
	{
	  counter++;
	  if(counter%1000 == 0) {if(counter > 0) cout << '\r'; cout << "\t" << int(double(counter)*100.0/pmax) << "% finished: " << counter << " out of " << pmax; cout.flush();}

	  double R2_min = (Sx-(Sx == 0 ? 0 : 1))*(Sx-(Sx == 0 ? 0 : 1))+(Sy-(Sy == 0 ? 0 : 1))*(Sy-(Sy == 0 ? 0 : 1))+(Sz-(Sz == 0 ? 0 : 1))*(Sz-(Sz == 0 ? 0 : 1));
	  double R2_max = (Sx+1)*(Sx+1)+(Sy+1)*(Sy+1)+(Sz+1)*(Sz+1);
	  
	  double w_eta = 0.0;

	  // The weight for point Sx,Sy,Sz has contributions from all adjoining cells
	  // The furthest corner is Sx+1,Sy+1,Sz+1 and if this less than hsr*hsr, the volume weights are 1, and the surface weights are zero
	  // else, the nearest corner is Sx-1,Sy-1,Sz-1 (unless Sx, Sy or Sz = 0) and if this is less than hsr*hsr, then the boundary is between these limits and we must compute
	  // else, all hsd boundary is less than the nearest corner and all weights are zero.

	  double R = hsr/dx;

	  //Note: Special cases of I-sums, e.g. when one or more components of S are zero, are handled in the called functions.
	  
	  if(R*R > R2_max) {w_eta = dV;}
	  else if(R*R > R2_min)
	    for(int ax:I)
	      {
		int ix = ax;
		double Tx = Sx+ix;
		int px = (Tx < 0 ? -1 : 1);
		if(Tx < 0) Tx = ix = 1;		   
		for(int ay:I)
		  {
		    int iy = ay;
		    double Ty = Sy+iy;
		    int py = (Ty < 0 ? -1 : 1);
		    if(Ty < 0) {Ty = 1; iy = 1;}	       		    
		    for(int az:I)
		      {
			int iz = az;
			double Tz = Sz+iz;
			int pz = (Tz < 0 ? -1 : 1);
			if(Tz < 0) Tz = iz = 1;
  
			int vx[2] = {0,ix};
			int vy[2] = {0,iy};
			int vz[2] = {0,iz};

			double j = 0.0;

			for(int jx = 0; jx < 2; jx++)
			  for(int jy = 0; jy < 2; jy++)
			    for(int jz = 0; jz < 2; jz++)
			      {
				double Vx = Tx - vx[jx];
				double Vy = Ty - vy[jy];
				double Vz = Tz - vz[jz];
			    
				int sgn = 1;
				if(jx == 1) sgn *= -1;
				if(jy == 1) sgn *= -1;
				if(jz == 1) sgn *= -1;
			    
				// test up to machine precision
				if((R*R - (Vx*Vx+Vy*Vy+Vz*Vz)) < std::nextafter(0.d,1.d)) continue;

				w_eta += dV*sgn*(G_eta(R,Vx,Vy,Vz,Tx,Ty,Tz) - G_eta(R,sqrt(R*R-Vy*Vy-Vz*Vz),Vy,Vz,Tx,Ty,Tz));
			      }
		      }
		  }
	      }	
	  
	  // Add in for all octants of the sphere: take account of parity of vector and tensor quantities
	  for(int ix = 0; ix < (Sx == 0 ? 1 : 2); ix++)
	    for(int iy = 0; iy < (Sy == 0 ? 1 : 2); iy++)
	      for(int iz = 0; iz < (Sz == 0 ? 1 : 2); iz++)
		{		  
		  long pos = density_->get_PBC_Pos((1-2*ix)*Sx,(1-2*iy)*Sy,(1-2*iz)*Sz);	  
		  eos_weighted_density_->addToWeight(pos,w_eta);
		  if(std::isnan(eos_weighted_density_->getWeight(pos)))
		    {
		      cout << ix << " " << iy << " " << iz << " " << Sx << " " << Sy << " " << Sz << endl;
		      throw std::runtime_error("Found NAN");
		    }
		}		  
	}

  cout << endl;
  cout << myColor::YELLOW;
  cout << "/////  Finished.  " << endl;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << myColor::RESET << endl;
}
*/


////////////////////////////
// AO model
FMT_AO_Species:: FMT_AO_Species(Density& density, double hsd, double Rp, double reservoir_density, double mu, int seq)
  : Rp_(Rp), reservoir_density_(reservoir_density), FMT_Species(density,hsd,mu,seq), fmt_weighted_densitiesAO_(11)
{
  // get the polymer species weights
  long Nx = density_->Nx();
  long Ny = density_->Ny();
  long Nz = density_->Nz();

  for(FMT_Weighted_Density &d: fmt_weighted_densitiesAO_)
    d.initialize(Nx, Ny, Nz);
  
  #if defined FMT_WEIGHTS_BEFORE_JUN_2021
    // The old function "generateWeights_before_jun_2021" could be removed 
    // as the new one does reproduce the same calculation to machine precision
    // by using here the arguments on the commented line below.
    generateWeights_before_jun_2021(2*Rp_, fmt_weighted_densitiesAO_);
    //generateWeights(2*Rp_, fmt_weighted_densitiesAO_, density_->getDX(), false);
  #elif defined FMT_WEIGHTS_BEFORE_JAN_2023
    // The old function "generateWeights_before_jan_2023" could be removed 
    // as the new one does reproduce the same calculation to machine precision
    // by using here the arguments on the commented line below.
    generateWeights_before_jan_2023(2*Rp_, fmt_weighted_densitiesAO_);
    //generateWeights(2*Rp_, fmt_weighted_densitiesAO_, 1.0, false);
  #else
    generateWeights(2*Rp_, fmt_weighted_densitiesAO_);
  #endif
  
  for(FMT_Weighted_Density &d: fmt_weighted_densitiesAO_)
    d.transformWeights();

  PSI_.initialize(Nx,Ny,Nz);
}


// The job of this function is to take the information concerning dF/dn_{a}(pos) (stored in DPHI) and to construct the dPhi_dn for partial measures that are stored in fmt_weighted_densities.
// This is done via the call to FMT_Species::set_fundamental_measure_derivatives. 
void FMT_AO_Species::set_fundamental_measure_derivatives(FundamentalMeasures &DPHI, long pos, bool needsTensor)
{
  FMT_Species::set_fundamental_measure_derivatives(DPHI,pos,needsTensor);

  double hsdp = 2*Rp_;
  
  double dPhi_dEta = DPHI.eta;
  double dPhi_dS = (DPHI.s0/(hsdp*hsdp)) + DPHI.s1*((1/hsdp)) + DPHI.s2;
  double dPhi_dV[3] = {DPHI.v2[0] + DPHI.v1[0]/hsdp,
			  DPHI.v2[1] + DPHI.v1[1]/hsdp,
			  DPHI.v2[2] + DPHI.v1[2]/hsdp};
  
  fmt_weighted_densitiesAO_[EI()].Set_dPhi(pos,dPhi_dEta);
  fmt_weighted_densitiesAO_[SI()].Set_dPhi(pos,dPhi_dS);

  //note parity factor for the vector
  for(int j=0;j<3;j++)
    {
      fmt_weighted_densitiesAO_[VI(j)].Set_dPhi(pos, -dPhi_dV[j]);	
      if(needsTensor)
	for(int k=j;k<3;k++)
	  fmt_weighted_densitiesAO_[TI(j,k)].Set_dPhi(pos,(j == k ? 1 : 2)*DPHI.T[j][k]); // taking account that we only use half the entries
    }
}

// This is where we finish the construction of Psi(r) and
// evaluate the contribution to the free energy
double FMT_AO_Species::free_energy_post_process(bool needsTensor)  
{  
  PSI_.Four().zeros();  
  
  // The call to add_to_dPhi s does the fft of dPhi/dn_{a} for each fm and adds to array PSI
  int number_of_weights = (needsTensor ? fmt_weighted_densities.size() : 5);

  int i1;  
  //#pragma omp parallel for   private(i1) schedule(static)	 
  for(i1=0;i1<number_of_weights;i1++)
    fmt_weighted_densitiesAO_[i1].add_weight_schur_dPhi_to_arg(PSI_.Four()); 

  PSI_.do_fourier_2_real();

  // Here, we compute the contribution to the free energy and begin to reuse PSI in preparation for computing the forces
  double F = 0;

  long Ntot = PSI_.cReal().size();
  long i;
#ifdef USE_OMP      
#pragma omp parallel for private(i) schedule(static) reduction(+:F)
#endif  
  for(long i=0;i<Ntot;i++)
    {
      double val = exp(-PSI_.cReal().get(i));
      PSI_.Real().set(i,val);
      F += val;
    }
  F *= -reservoir_density_;

  // This is the "standard" shift of the free energy. The Ntot eventually gets multiplied by dV to become V. 
  F += reservoir_density_*density_->Ntot();
  
  // This is to prepare for the force calculation which comes later. Strictly speaking, it is not a part of the free energy calculation. 
  // Upsilon requires convoluting the (now exponentiated) PSI with the various weights
  // We do this here because the next step is carried out in FMT::calculateFreeEnergyAndDerivatives and so PSI_ is not accessible.
  PSI_.do_real_2_fourier(); // do FFT

  for(auto &x: fmt_weighted_densitiesAO_)
    x.convoluteWith(PSI_.cFour()); // point-wise multiplication of FFT's and call to fourier_2_real (norm factor was included in definition of weights)
  
  return F;
}
