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
#include "FMT.h"


FMT_Gaussian_Species::FMT_Gaussian_Species(GaussianDensity& density, double hsd, double mu, int seq)
  : FMT_Species(density,hsd,mu,seq)
{
  density.set_hsd(hsd);
}

void FMT_Gaussian_Species::set_density_from_alias(const DFT_Vec &x)
{
  GaussianDensity &density = *(static_cast<GaussianDensity*>(&density_));

  long pos = 0;
  if(density.use_discrete())
    {
      Species::set_density_from_alias(x);
      pos += density.Ntot();
    }

  for(int ig=0; ig<density.number_of_gaussians(); ig++, pos += 5)  
    density.set_gaussian(ig, x.get(pos), x.get(pos+1), x.get(pos+2), x.get(pos+3), x.get(pos+4));
}

  
void FMT_Gaussian_Species::get_density_alias(DFT_Vec &x) const
{
  GaussianDensity &density = *(static_cast<GaussianDensity*>(&density_));

  long pos = 0;
  if(density.use_discrete())
    {
      FMT_Species::get_density_alias(x);
      pos += density.Ntot();
    }

  for(int ig=0; ig<density.number_of_gaussians(); ig++, pos += 5)  
    {
      double y,b,Rx,Ry,Rz;

      density.get_gaussian(ig, y,b,Rx,Ry,Rz);

      x.set(pos,   y);
      x.set(pos+1, b);
      x.set(pos+2, Rx);
      x.set(pos+3, Ry);
      x.set(pos+4, Rz);
    }        
}

double FMT_Gaussian_Species::calculateFreeEnergyAndDerivatives_IdealGas_()
{
  GaussianDensity &density = *(static_cast<GaussianDensity*>(&density_));

  double F = 0.0;
  const double dV = density.dV();
  long pos = 0;

  if(density.use_discrete())
    {
#pragma omp parallel for  private(pos)  shared(dV) schedule(static)				
      for(pos = 0; pos < density.Ntot(); pos++)
	{
	  double dg = density.get_discrete_gauss_at(pos);
	  double dd = density.get(pos);

	  double ld = log(dd+dg);
	  
	  F += ((dd+dg)*ld - dg*log(dg) - dd)*dV;

	  dF_.IncrementBy(pos,ld*dV);
	}
    }
  
  for(int ig=0; ig<density.number_of_gaussians(); ig++, pos += 5)
    {
      const Gaussian &g = density.get_gaussian(ig);
      
      double x = g.prefactor();
      double alf = g.alf();

      double dfdx = 0;
      double dfdb = 0;
      
      F += g.get_ideal_f(dfdx, dfdb);
      // flatten 
      dF_.IncrementBy(pos+0, dfdx);
      dF_.IncrementBy(pos+1, dfdb);      
    }

  return F;
}

void FMT_Gaussian_Species::calculateFundamentalMeasures(bool needsTensor)
{
  if(use_discrete)
    FMT_Species.calculateFundamentalMeasures(needsTensor);
  
  long pos;
#pragma omp parallel for  private(pos)  schedule(static)				
  for(pos = 0; pos < density_.Ntot(); pos++)
    {
      FundamentalMeasures fm;
      dynamic_cast<GaussianDensity*>(&density_)->get_measures(pos,hsd_,fm);

      fmt_weighted_densities[EI()].add_fundamental_measure(pos,fm.eta);
      fmt_weighted_densities[SI()].add_fundamental_measure(pos,fm.s2);
      for(int j=0;j<3;j++)
	{
	  fmt_weighted_densities[VI(j)].add_fundamental_measure(pos,fm.v2[j]);
	  if(needsTensor)
	    for(int k=j;k<3;k++)
	      fmt_weighted_densities[TI(j,k)].add_fundamental_measure(pos,fm.T[j][k]);
	}
    }
}

  // The only reason for this is that the FMT function includes a "parity switch" which is not needed here.
void FMT_Gaussian_Species::set_fundamental_measure_derivatives(FundamentalMeasures &DPHI, long pos, bool needsTensor)
{
  double dPhi_dEta = DPHI.eta;
  double dPhi_dS = (DPHI.s0/(hsd_*hsd_)) + (DPHI.s1/hsd_) + DPHI.s2;
  double dPhi_dV[3] = {DPHI.v2[0] + DPHI.v1[0]/hsd_,
		       DPHI.v2[1] + DPHI.v1[1]/hsd_,
		       DPHI.v2[2] + DPHI.v1[2]/hsd_};

  fmt_weighted_densities[EI()].Set_dPhi(pos,dPhi_dEta);
  fmt_weighted_densities[SI()].Set_dPhi(pos,dPhi_dS);    

  for(int j=0;j<3;j++)
    {
      fmt_weighted_densities[VI(j)].Set_dPhi(pos, dPhi_dV[j]);	
      if(needsTensor)
	for(int k=j;k<3;k++)
	  fmt_weighted_densities[TI(j,k)].Set_dPhi(pos,(j == k ? 1 : 2)*DPHI.T[j][k]); // taking account that we only use half the entries
    }
}

void FMT_Gaussian_Species::Build_Force(bool needsTensor)
{
  // HERE
  // what to do about Species::Build_Force????

  
  GaussianDensity &density = *(static_cast<GaussianDensity*>(&density_));

  double dV = density.dV();    
  
  for(int ig=0; ig<density.number_of_gaussians(); ig++)
    {
      const Gaussian &g = density.get_gaussian(ig);
      
      for(long pos = 0; pos < density.Ntot(); pos++)
	{
	  // turn lattice position into spatial position
	  long ix,iy,iz;
	  density.cartesian(pos,ix,iy,iz);

	  double rx = density.getX(ix);
	  double ry = density.getY(iy);
	  double rz = density.getZ(iz);

	  // get derivatives
	  
	  FundamentalMeasures dfm[5];
	  density.get_dmeasures_for_gaussian(ig, rx, ry, rz, hsd_, dfm);

	  FundamentalMeasures fm1;
	  density.get_measures(rx, ry, rz, hsd_, fm1);

	  
	  double df = 0.0;
	  // here, the dPhi that is held by each weighted measure is actually the "collapsed" version:
	  // so, e.g. fmt_weighted_densities[1].dPhi is  dPhi_ds0/hsd*hsd + dPhi_ds1/hsd + dPhi_ds2

	  for(int j=0;j<5;j++)
	    {
	      df = fmt_weighted_densities[EI()].get_dPhi(pos);	  
	      dF_.IncrementBy(5*ig+j,df*dfm[j].eta*dV);

	      df = fmt_weighted_densities[SI()].get_dPhi(pos);	  	  
	      dF_.IncrementBy(5*ig+j,df*dfm[j].s2*dV);
	  
	      df = fmt_weighted_densities[VI(0)].get_dPhi(pos);	  	  	  
	      dF_.IncrementBy(5*ig+j,df*dfm[j].v2[0]*dV);

	      df = fmt_weighted_densities[VI(1)].get_dPhi(pos);	  	  	  
	      dF_.IncrementBy(5*ig+j,df*dfm[j].v2[1]*dV);

	      df = fmt_weighted_densities[VI(2)].get_dPhi(pos);	  	  	  
	      dF_.IncrementBy(5*ig+j,df*dfm[j].v2[2]*dV);

	      if(needsTensor)
		{
		  df = fmt_weighted_densities[TI(0,0)].get_dPhi(pos);	  	  	  
		  dF_.IncrementBy(5*ig+j,df*dfm[j].T[0][0]*dV);
		  
		  df = fmt_weighted_densities[TI(0,1)].get_dPhi(pos);	  	  	  
		  dF_.IncrementBy(5*ig+j,df*dfm[j].T[0][1]*dV);

		  df = fmt_weighted_densities[TI(0,2)].get_dPhi(pos);	  	  	  
		  dF_.IncrementBy(5*ig+j,df*dfm[j].T[0][2]*dV);
		  
		  df = fmt_weighted_densities[TI(1,1)].get_dPhi(pos);	  	  	  
		  dF_.IncrementBy(5*ig+j,df*dfm[j].T[1][1]*dV);
		  
		  df = fmt_weighted_densities[TI(1,2)].get_dPhi(pos);	  	  	  
		  dF_.IncrementBy(5*ig+j,df*dfm[j].T[1][2]*dV);

		  df = fmt_weighted_densities[TI(2,2)].get_dPhi(pos);	  	  	  
		  dF_.IncrementBy(5*ig+j,df*dfm[j].T[2][2]*dV);
		}
	    }
	}
    }
}


double FMT_Gaussian_Species::externalField(bool bCalcForces)
{
  double Fx =  -mu_*density_.getNumberAtoms();
  if(bCalcForces)
    {
      GaussianDensity &density = *(static_cast<GaussianDensity*>(&density_));
      for(int ig=0; ig<density.number_of_gaussians(); ig++)
	{
	  double dN_dx;
	  double dN_dalf;
	  density.get_dN(ig,dN_dx, dN_dalf);

	  dF_.IncrementBy(5*ig  , -mu_*dN_dx);
	  dF_.IncrementBy(5*ig+1, -mu_*dN_dalf);
	  //	  dF_.set(5*ig+1, 100000*dF_.get(5*ig+1)); // SCALE ALPHA FORCE
	  //	  dF_.set(5*ig+1, dF_.get(5*ig+1)); // SCALE ALPHA FORCE
	}
    }
  return Fx;  
}

void FMT_Gaussian_Species::beginForceCalculation()
{
  if(fixedMass_ > 0.0)
    throw std::runtime_error("Fixed mass not implemented in FMT_Gaussian_Species::beginForceCalculation");

  GaussianDensity &density = *(static_cast<GaussianDensity*>(&density_));
  if(density.use_discrete())
    density.fill_discrete_gaussian_field();
  
}

void FMT_Gaussian_Species::addToGaussianForce(vector<double> &dF)
{
  GaussianDensity &density = *(static_cast<GaussianDensity*>(&density_));  
  long pos = (density.use_discrete() ? density.Ntot() : 0);
  
  for(int i=0;i<dF.size();i++)
    dF_.IncrementBy(pos+i, dF[i]);

}
