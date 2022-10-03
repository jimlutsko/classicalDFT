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

int Species::SequenceNumber_ = 0;


// These functions just alias the density so that it is always non-zero and in fact bigger than SMALL_VALUE.
/*
void Species::set_density_from_alias(const DFT_Vec &x)
{
  long pos;
  const double dmin = SMALL_VALUE;
  
#pragma omp parallel for  private(pos)  schedule(static)
  for(pos=0;pos<x.size();pos++)
    density_.set(pos,dmin+x.get(pos)*x.get(pos));
}
  
void Species::get_density_alias(DFT_Vec &x) const
{
  long pos;
  const double dmin = SMALL_VALUE;
  
#pragma omp parallel for  private(pos)  schedule(static)				    
  for(pos=0;pos<x.size();pos++)
    x.set(pos, sqrt(std::max(0.0, density_.get(pos)-1e-20)));    
}

void Species::convert_to_alias_deriv(DFT_Vec &x, DFT_Vec &dF_dRho) const
{
  dF_dRho.Schur(x,dF_dRho);
  dF_dRho.MultBy(2.0);
}
*/

// This is another alias which maps the zero of the density 
// (or rathera minimal value SMALL_VALUE) to -inf and the
// maximal value physically possible (1/dV) to +inf
void Species::set_density_from_alias(const DFT_Vec &x)
{
  long pos;
  const double rho_max = 1/density_.dV();
  
  #pragma omp parallel for  private(pos)  schedule(static)
  for(pos=0;pos<x.size();pos++)
    density_.set(pos, SMALL_VALUE + 0.5*rho_max*tanh(2*x.get(pos)/rho_max) );
}
  
void Species::get_density_alias(DFT_Vec &x) const
{
  long pos;
  const double rho_max = 1/density_.dV();
  
  #pragma omp parallel for  private(pos)  schedule(static)				    
  for(pos=0;pos<x.size();pos++)
    x.set(pos, 0.5*rho_max*atanh( 2*(density_.get(pos)-SMALL_VALUE)/rho_max -1 ) );    
}

void Species::convert_to_alias_deriv(DFT_Vec &x, DFT_Vec &dF_dRho) const
{
  dF_dRho.zeros(x.size());
  
  long pos;
  const double rho_max = 1/density_.dV();
  
  #pragma omp parallel for  private(pos)  schedule(static)				    
  for(pos=0;pos<x.size();pos++)
    dF_dRho.set(pos, dF_dRho.get(pos) / pow( cosh(2*x.get(pos)/rho_max) ,2) );
}


double Species::externalField(bool bCalcForces)
{
  double dV = density_.dV();
  double Fx = density_.get_field_dot_density()*dV - density_.getNumberAtoms()*mu_;
  if(bCalcForces)
    {
      dF_.IncrementBy_Scaled_Vector(density_.get_external_field(),dV);
      dF_.add(-mu_*dV);
    }
  return Fx;
}


void Species::beginForceCalculation()
{
  if(fixedMass_ > 0.0)
    {
      density_.scale_to(fixedMass_/density_.getNumberAtoms());
      mu_ = 0.0;
    }
}

double Species::endForceCalculation()
{
  if(fixedBackground_ && fixedMass_ > 0.0)
    throw std::runtime_error("Cannot have both fixed background and fixed mass .... aborting");
        
  if(fixedBackground_)
    {
      for(long pos = 0; pos < density_.get_Nboundary(); pos++)
	dF_.set(density_.boundary_pos_2_pos(pos),0.0);	
    }

  if(homgeneousBoundary_)
    {
	
      double average_border_force = 0;

      for(long pos = 0; pos < density_.get_Nboundary(); pos++)
	average_border_force += dF_.get(density_.boundary_pos_2_pos(pos));

      average_border_force /= density_.get_Nboundary();

      for(long pos = 0; pos < density_.get_Nboundary(); pos++)
	dF_.set(density_.boundary_pos_2_pos(pos),average_border_force);
    }    

    
  if(fixedMass_ > 0.0)
    {
      mu_ = 0.0;

      double Mtarget = fixedMass_;

      for(long p=0;p<density_.Ntot();p++)
	mu_ += dF_.get(p)*density_.get(p);
      mu_ /= Mtarget; //fixedMass_;
      for(long p=0;p<density_.Ntot();p++)
	dF_.set(p, dF_.get(p)-mu_*density_.dV());
    }
  return 0;
}
