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
DFT_VDW_Surfactant<T>::DFT_VDW_Surfactant(Lattice &lattice, Potential1 &potential,  string& pointsFile, double kT, double Asurf, double rhosurf)
  : DFT_VDW<T>(lattice, potential, pointsFile, kT), Asurf_(Asurf), rho_surf_(rhosurf),
  surfactant_density_(lattice.Nx(), lattice.Ny(), lattice.Nz()),
  surfactant_potential_(lattice.Nx(), lattice.Ny(), lattice.Nz()), bFixedN_(false)
{

  int Nx = lattice.Nx();
  int Ny = lattice.Ny();
  int Nz = lattice.Nz();

  double dx = lattice.getDX();
  double dy = lattice.getDY();
  double dz = lattice.getDZ();
  
  // Construct surfactant potential

  cout << "Asurf = " << Asurf_ << " rho_surf_ = " << rho_surf_ << endl;
  
  double r20 = pow(2,1.0/3);

  for(int ix=0;ix<Nx;ix++)
    for(int iy=0;iy<Ny;iy++)
      for(int iz=0;iz<Nz;iz++)
	{
	  long i = lattice.pos(ix,iy,iz);

	  //	      double vpot = (i == 0 ? Asurf_*1/dV  : 0);

	  int dix = ix; if(dix > Nx/2) dix -= Nx; if(dix < -Nx/2) dix += Nx;
	  int diy = iy; if(diy > Ny/2) diy -= Ny; if(diy < -Ny/2) diy += Ny;
	  int diz = iz; if(diz > Nz/2) diz -= Nz; if(diz < -Nz/2) diz += Nz; 

	  double r2 =  dix*dix*dx*dx+diy*diy*dy*dy+diz*diz*dz*dz;

	  double vpot = 0;
	  if(r2 < r20) r2 = r20;
	  if(r2 < 6*6) vpot = Asurf_*(1/(r2*r2*r2*r2*r2*r2)-(1/(r2*r2*r2)));
	  surfactant_potential_.Real().set(i, vpot);	  
	}
  surfactant_potential_.do_real_2_fourier();

}


template <class T>
double DFT_VDW_Surfactant<T>::calculateFreeEnergyAndDerivatives(Density& density, double mu, DFT_Vec& dF, bool onlyFex)
{
  double F = DFT_VDW<T>::calculateFreeEnergyAndDerivatives(density, mu, dF,onlyFex);

  int Nx = density.Nx();
  int Ny = density.Ny();
  int Nz = density.Nz();

  double dx = density.getDX();
  double dy = density.getDY();
  double dz = density.getDZ();

  long Ntot = Nx*Ny*Nz;

  double dV = density.dV();

  const DFT_Vec& vReal0 = DFT_VDW<T>::dft_fmt_->getV_Real(0);
  const DFT_Vec& vReal1 = DFT_VDW<T>::dft_fmt_->getV_Real(1);
  const DFT_Vec& vReal2 = DFT_VDW<T>::dft_fmt_->getV_Real(2);

  cout << "Asurf = " << Asurf_ << " rho_surf_ = " << rho_surf_ << endl;
  
  // Construct v2
  double r20 = pow(2,1.0/3);
  DFT_FFT v2(Nx, Ny, Nz);      
  for(int ix=0;ix<Nx;ix++)
    for(int iy=0;iy<Ny;iy++)
      for(int iz=0;iz<Nz;iz++)
	{
	  long i = density.pos(ix,iy,iz);
	  v2.Real().set(i, vReal0.get(i)*vReal0.get(i)+vReal1.get(i)*vReal1.get(i)+vReal2.get(i)*vReal2.get(i));
	}

  v2.do_real_2_fourier();      
      
  // Construct surfactant density      
  surfactant_density_.Four().Schur(surfactant_potential_.Four(),v2.Four());
  surfactant_density_.do_fourier_2_real();
  for(long i=0;i<Ntot;i++)
    surfactant_density_.Real().set(i, exp(-surfactant_density_.Real().get(i)*dV/Ntot));

  double rawSum = surfactant_density_.Real().accu();
      
  if(bFixedN_)
    surfactant_density_.Real().multBy(Ntot/rawSum);

  surfactant_density_.Real().multBy(rho_surf_);      
  surfactant_density_.do_real_2_fourier();

  // Free energy is just the integral
  double Fs = -(surfactant_density_.Real().accu() - rho_surf_*Ntot);
  if(bFixedN_)
    Fs = -rho_surf_*Ntot*log(rawSum/Ntot);
      
  // Construct conv(surfactant_density, surfactant_potential_) 
  DFT_FFT Density_Conv_Potential(Nx,Ny,Nz);      
  Density_Conv_Potential.Four().Schur(surfactant_density_.Four(),surfactant_potential_.Four());   
  Density_Conv_Potential.do_fourier_2_real(); 
  Density_Conv_Potential.Real().multBy(dV/Ntot);
  
  // Construct conv(surfactant_density, surfactant_potential_)(i) v(i)
  // and convolute with vector weight function            
  DFT_FFT dFs(Nx,Ny,Nz);
  dFs.Four().zeros();            
  for(int J=0;J<3;J++)
    {
      DFT_FFT v_J(Nx, Ny, Nz);
      v_J.Real().Schur(Density_Conv_Potential.Real(), (J == 0 ? vReal0 : (J == 1 ? vReal1 : vReal2)));
      v_J.do_real_2_fourier();
      dFs.Four().incrementSchur(v_J.Four(),DFT_VDW<T>::dft_fmt_->getVweight_Four(J));
    }
  dFs.do_fourier_2_real();
  dFs.Real().multBy(2*dV); // N.B.: No factor of Ntot because it is already in wk() ...

  //Done! Add to acculated contributions from previous calculations
  F += Fs*dV; 
  dF.Increment_And_Scale(dFs.cReal(),1.0); 

  if(bFixedN_)
    cout << "Number of surfactant atoms in cell = " << surfactant_density_.Real().accu()*dV << endl;

  return F;
}



template <class T>
double DFT_VDW_Surfactant<T>::Mu(double x) const
{
  return DFT_VDW<T>::Mu(x);
}

template <class T>
double DFT_VDW_Surfactant<T>::Omega(double x) const
{
  return DFT_VDW<T>::Omega(x);
}

template <class T>
double DFT_VDW_Surfactant<T>::Fhelmholtz(double x) const
{
  return DFT_VDW<T>::Fhelmholtz(x);
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
