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

#include "DFT.h"


template <class T>
DFT_VDW_Surfactant<T>::DFT_VDW_Surfactant(int Nx, int Ny, int Nz)
  : DFT_VDW<T>(Nx,Ny,Nz),  bFixedN_(false)
{}

template <class T>
void DFT_VDW_Surfactant<T>::addSurfactantPotential(Potential1 &pot, double kT)
{
  Species *species = ((VDW_Species*) DFT_VDW_Surfactant<T>::allSpecies_.front());
  
  int Nx = species->getLattice().Nx();
  int Ny = species->getLattice().Ny();
  int Nz = species->getLattice().Nz();

  double dx = species->getLattice().getDX();
  double dy = species->getLattice().getDY();
  double dz = species->getLattice().getDZ();

  surfactant_potential_.initialize(Nx,Ny,Nz);

  // Construct surfactant potential

  for(int ix=0;ix<Nx;ix++)
    for(int iy=0;iy<Ny;iy++)
      for(int iz=0;iz<Nz;iz++)
	{
	  long i = species->getLattice().pos(ix,iy,iz);

	  int dix = ix; if(dix > Nx/2) dix -= Nx; if(dix < -Nx/2) dix += Nx;
	  int diy = iy; if(diy > Ny/2) diy -= Ny; if(diy < -Ny/2) diy += Ny;
	  int diz = iz; if(diz > Nz/2) diz -= Nz; if(diz < -Nz/2) diz += Nz; 

	  double r =  sqrt(dix*dix*dx*dx+diy*diy*dy*dy+diz*diz*dz*dz);

	  surfactant_potential_.Real().set(i, pot.Watt(r));	  
	}
  surfactant_potential_.do_real_2_fourier();
}

template <class T>
double DFT_VDW_Surfactant<T>::calculateFreeEnergyAndDerivatives(bool onlyFex)
{
  double F = DFT_VDW<T>::calculateFreeEnergyAndDerivatives(onlyFex);

  // I will assume that species 0 is the **water** and species 1 is the surfactant.

  VDW_Species &water = *((VDW_Species*) DFT_VDW_Surfactant<T>::allSpecies_[0]);
  VDW_Species &surf = *((VDW_Species*) DFT_VDW_Surfactant<T>::allSpecies_[1]);
  
  const Density& surf_density = surf.getDensity();

  int Nx = surf_density.Nx();
  int Ny = surf_density.Ny();
  int Nz = surf_density.Nz();

  long Ntot = Nx*Ny*Nz;
  
  double dx = surf_density.getDX();
  double dy = surf_density.getDY();
  double dz = surf_density.getDZ();

  double dV = surf_density.dV();

  ////////////////////////////////////////////////////////////////////////////
  // Construct v2(r) for the water & do fft

  const DFT_Vec& vRealx = water.getV_Real(0);
  const DFT_Vec& vRealy = water.getV_Real(1);
  const DFT_Vec& vRealz = water.getV_Real(2);

  DFT_FFT v2(Nx, Ny, Nz);      
  for(long i=0;i<Ntot;i++)
    v2.Real().set(i, vRealx.get(i)*vRealx.get(i)+vRealy.get(i)*vRealy.get(i)+vRealz.get(i)*vRealz.get(i));
  v2.do_real_2_fourier();

  //////////////////////////////////////////////////////////////////////////
  // Surfactant force is just convoluton of the potential and v2
  // Lets do this the dumb way first:
  DFT_FFT dF(Nx, Ny, Nz);      

  dF.Four().Schur(surfactant_potential_.Four(),v2.Four(),false); //true);
  dF.Four().scaleBy(dF.Real().size());  // FFTW's convention
  dF.do_fourier_2_real();
  dF.Real().multBy(dV*dV); // this is d F/d rho_i and so the factors of dV survive

  surf.addToForce(dF.Real());
  
  /////////////////////////////////////////////////////////////////////////////////
  // Contribution to the free energy is this dotted with the surfactant density

  F += surf_density.getInteractionEnergy(dF.Real());


  ////////////////////////////////////////////////////////////////////////
  // Now we need the derivative wrt the water density

  for(int i=0;i<3;i++)
    { 
      dF.zeros();

      dF.Four().Schur(surf_density.getDK(), surfactant_potential_.Four(),false);
      dF.Four().scaleBy(dF.Real().size());  // FFTW's convention
      dF.do_fourier_2_real();  

      dF.Real().Schur(dF.Real(),water.getV_Real(i));

      dF.do_real_2_fourier();
      dF.Four().Schur(dF.Four(),water.getVweight_Four(i),false); // this includes the FFTW scale factor already
      dF.do_fourier_2_real();

      dF.Real().multBy(2*dV*dV);
      
      water.addToForce(dF.cReal());
    }

  return F;
}


template <class T>
void DFT_VDW_Surfactant<T>::coexistence(double& xliq, double &xgas) const
{
  DFT_VDW<T>::coexistence(xliq,xgas);
}

template <class T>
void DFT_VDW_Surfactant<T>::spinodal(double& xs1, double &xs2) const
{
  DFT_VDW<T>::spinodal(xs1,xs2);
}

template <class T>
double DFT_VDW_Surfactant<T>::findLiquidFromMu(double mu, double mu_coex, double xliq_coex) const
{
  return DFT_VDW<T>::findLiquidFromMu(mu, mu_coex, xliq_coex);
}


template class DFT_VDW_Surfactant<WhiteBearI>;
template class DFT_VDW_Surfactant<WhiteBearII>;
template class DFT_VDW_Surfactant<RSLT>;
template class DFT_VDW_Surfactant<RSLT2>;
