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


void FMT_Gaussian_Species::set_density_from_alias(const DFT_Vec &x)
{
  GaussianDensity &density = *(static_cast<GaussianDensity*>(&density_));
  double L[] = {density.Lx(), density.Ly(), density.Lz()};
  
  for(int ig=0; ig<density.gaussians_.size(); ig++)  
    {
      Gaussian &g = density.gaussians_[ig];
      long pos = ig*5;
      g.set_parameters(x.get(pos), x.get(pos+1), x.get(pos+2), x.get(pos+3), x.get(pos+4), hsd_, L);
    }      
}

  
void FMT_Gaussian_Species::get_density_alias(DFT_Vec &x) const
{
  GaussianDensity &density = *(static_cast<GaussianDensity*>(&density_));
  double L[] = {density.Lx(), density.Ly(), density.Lz()};
  
  for(int ig=0; ig<density.gaussians_.size(); ig++)  
    {
      Gaussian &g = density.gaussians_[ig];
      double y,alf,Rx,Ry,Rz;

      g.get_parameters(y,alf,Rx,Ry,Rz,hsd_);

      long pos = ig*5;
      x.set(pos,   y);
      x.set(pos+1, alf);
      x.set(pos+2, Rx);
      x.set(pos+3, Ry);
      x.set(pos+4, Rz);
    }        
}

double FMT_Gaussian_Species::calculateFreeEnergyAndDerivatives_IdealGas_()
{
  GaussianDensity &density = *(static_cast<GaussianDensity*>(&density_));

  double F = 0.0;

  for(int ig=0; ig<density.gaussians_.size(); ig++)
    {
      Gaussian &g = density.gaussians_[ig];
      
      double x = g.prefactor();
      double alf = g.alf();
            
      F += x*log(x) + 1.5*x*log(alf/M_PI)-2.5*x;

      double dF_dX   = log(x) + 1.0 + 1.5*log(alf/M_PI)-2.5;
      double dF_dAlf = 1.5*x/alf;

      double dX_dY   = g.dprefactor_dx_;
      double dX_dAlf = g.dprefactor_dalf_;
      
      // flatten 
      int pos = 5*ig;
      dF_.IncrementBy(pos,   dF_dX*dX_dY);
      dF_.IncrementBy(pos+1, dF_dAlf + dF_dX*dX_dAlf);
    }
  return F;
}

void FMT_Gaussian_Species::calculateFundamentalMeasures(bool needsTensor)
{
  long pos;
#pragma omp parallel for  private(pos)  schedule(static)				
  for(pos = 0; pos < density_.Ntot(); pos++)
    {
      FundamentalMeasures fm;
      dynamic_cast<GaussianDensity*>(&density_)->get_measures(pos,hsd_,fm);

      fmt_weighted_densities[EI()].set_fundamental_measure(pos,fm.eta);
      fmt_weighted_densities[SI()].set_fundamental_measure(pos,fm.s2);
      for(int j=0;j<3;j++)
	{
	  fmt_weighted_densities[VI(j)].set_fundamental_measure(pos,fm.v2[j]);
	  if(needsTensor)
	    for(int k=j;k<3;k++)
	      fmt_weighted_densities[TI(j,k)].set_fundamental_measure(pos,fm.T[j][k]);
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
  GaussianDensity &density = *(static_cast<GaussianDensity*>(&density_));

  double dV = density.dV();    

  
  for(int ig=0; ig<density.gaussians_.size(); ig++)
    {
      Gaussian &g = density.gaussians_[ig];
      
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
      for(int ig=0; ig<density.gaussians_.size(); ig++)
	{
	  double dN_dx;
	  double dN_dalf;
	  density.get_dN(ig,dN_dx, dN_dalf);

	  dF_.IncrementBy(5*ig  , -mu_*dN_dx);
	  dF_.IncrementBy(5*ig+1, -mu_*dN_dalf);
	  dF_.set(5*ig+1, 100000*dF_.get(5*ig+1)); // SCALE ALPHA FORCE
	}
    }
  return Fx;  
}

double FMT_Gaussian_Species::FMF(double w0, double r0, vector<double> &x, vector<double> &w, DFT_Vec &dF) const
{
  double F = 0;

  GaussianDensity &density = *(static_cast<GaussianDensity*>(&density_));
  //  vector<double> dF(5*density.gaussians_.size());

  dF.zeros(5*density.gaussians_.size());
  F = density.FMF(w0,r0,x,w,dF);

  
  return F;
}
