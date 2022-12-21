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

BOOST_CLASS_EXPORT(Species)
BOOST_CLASS_EXPORT(FMT_Species)





int Species::SequenceNumber_ = 0;


// These functions just alias the density so that it is always non-zero and in fact bigger than SMALL_VALUE.
void Species::set_density_from_alias(const DFT_Vec &x)
{
  long pos;
  const double dmin = SMALL_VALUE;
  
#pragma omp parallel for  private(pos)  schedule(static)
  for(pos=0;pos<x.size();pos++)
    density_->set(pos,dmin+x.get(pos)*x.get(pos));
}
  
void Species::get_density_alias(DFT_Vec &x) const
{
  long pos;
  const double dmin = SMALL_VALUE;
  
#pragma omp parallel for  private(pos)  schedule(static)				    
  for(pos=0;pos<x.size();pos++)
    x.set(pos, sqrt(std::max(0.0, density_->get(pos)-1e-20)));    
}

void Species::convert_to_alias_deriv(DFT_Vec &x, DFT_Vec &dF_dRho) const
{
  dF_dRho.Schur(x,dF_dRho);
  dF_dRho.MultBy(2.0);
}

/*
double Species::externalField(bool bCalcForces)
{
  double dV = density_->dV();
  double Fx = density_->get_field_dot_density()*dV - density_->getNumberAtoms()*mu_;
  if(bCalcForces)
    {
      dF_.IncrementBy_Scaled_Vector(density_->get_external_field(),dV);
      dF_.add(-mu_*dV);
    }
  return Fx;
}
*/

double Species::evaluate_contribution_chemical_potential()
{
  double dV = density_->dV();
  dF_.add(-mu_*dV);
  return  - density_->get_mass()*mu_;
}

double Species::evaluate_external_field(const External_Field &f)
{
  double dV = density_->dV();

  dF_.IncrementBy_Scaled_Vector(f.get_field(),dV);
  //  dF_.add(-mu_*dV);

  return f.get_field().dotWith(density_->get_density_real())*dV; // - density_->getNumberAtoms()*mu_;
}



void Species::beginForceCalculation()
{
  if(fixedMass_ > 0.0)
    {
      //      density_->scale_to(fixedMass_/density_->getNumberAtoms());
      *density_ *= (fixedMass_/density_->getNumberAtoms());
      mu_ = 0.0;
    }
}

double Species::endForceCalculation()
{
  //  if(fixedBackground_ && fixedMass_ > 0.0)
  //    throw std::runtime_error("Cannot have both fixed background and fixed mass .... aborting");
        
  if(fixedBackground_)
    {
      for(long pos = 0; pos < density_->get_Nboundary(); pos++)
	dF_.set(density_->boundary_pos_2_pos(pos),0.0);	
    }

  if(homogeneousBoundary_)
    {
	
      double average_border_force = 0;

      for(long pos = 0; pos < density_->get_Nboundary(); pos++)
	average_border_force += dF_.get(density_->boundary_pos_2_pos(pos));

      average_border_force /= density_->get_Nboundary();

      for(long pos = 0; pos < density_->get_Nboundary(); pos++)
	dF_.set(density_->boundary_pos_2_pos(pos),average_border_force);
    }    

    
  if(fixedMass_ > 0.0)
    {
      mu_ = 0.0;

      double Mtarget = fixedMass_;

      for(long p=0;p<density_->Ntot();p++)
	//	if(!fixedBackground_ || !(density_->is_boundary_point(p)))
	  mu_ += dF_.get(p)*density_->get(p);
      mu_ /= Mtarget; //fixedMass_;
      for(long p=0;p<density_->Ntot();p++)
	//	if(!fixedBackground_ || !(density_->is_boundary_point(p)))
	  dF_.set(p, dF_.get(p)-mu_*density_->dV());
    }

  // 
  //  if(fixedBackground_ && fixedMass_ > 0.0) // need to do this again for case of both constraints ...
  //    {
  //      for(long pos = 0; pos < density_->get_Nboundary(); pos++)
  //	dF_.set(density_->boundary_pos_2_pos(pos),0.0);	
  //    }

  
  return 0;
}


