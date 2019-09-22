/* This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
 * To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
 *
 * Author: James F. Lutsko
 * www.lutsko.com
 */

#ifndef __LUTSKO__INTERACTION__
#define __LUTSKO__INTERACTION__

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <complex>
#include <complex.h>

#include "Species.h"

/**
  *  @brief This class encapsulates the interaction between two species (or one species with itself)
  */  

class Interaction
{
 public:
 Interaction(Species &s1, Species &s2, Potential1 &v, double kT) : s1_(s1), s2_(s2)
    {
      const Density &density = s1.getDensity();
      
      // The lattice
      long Nx = density.Nx();
      long Ny = density.Ny();
      long Nz = density.Nz();

      double dx = density.getDX();
      double dy = density.getDY();
      double dz = density.getDZ();
  
      // Set up the mean field potential
      w_att_.initialize(Nx,Ny,Nz);

      for(int nx = 0;nx<Nx;nx++)
	for(int ny = 0;ny<Ny;ny++)
	  for(int nz = 0;nz<Nz;nz++)
	    {
	      long pos = nz+Nz*(ny+Ny*nx);

	      double x = nx*dx;
	      double y = ny*dy;
	      double z = nz*dz;

	      if(nx > Nx/2) x -= Nx*dx;
	      if(ny > Ny/2) y -= Ny*dy;
	      if(nz > Nz/2) z -= Nz*dz;

	      double r2 = x*x+y*y+z*z;
	      w_att_.Real().addTo(pos,v.Watt(sqrt(r2))/kT);
	    }
      // Set the parameters in the VDW object  
      // Do this BEFORE the FFT which may corrupt the real-space part
      a_vdw_ = 0.5*w_att_.Real().accu()*dx*dy*dz;

      // Now save the FFT of the field  
      w_att_.do_real_2_fourier();     

    }

  double getInteractionEnergyAndForces()
  {
    const Density &density1 = s1_.getDensity();
    
    long Ntot = density1.Ntot();
    double dV = density1.dV();

    DFT_FFT v(density1.Nx(), density1.Ny(), density1.Nz());      
    double E = 0;
    
    if(s1_.getSequenceNumber() == s2_.getSequenceNumber())
      {
	v.Four().Schur(density1.getDK(),w_att_.Four());
	v.Four().multBy(dV*dV/Ntot);
	v.do_fourier_2_real(); 
	s1_.addToForce(v.Real());
	E = 0.5*density1.getInteractionEnergy(v.Real());	
      } else {
	v.Four().Schur(density1.getDK(),w_att_.Four());
	v.Four().multBy(0.5*dV*dV/Ntot);
	v.do_fourier_2_real(); 
	s2_.addToForce(v.Real());

	const Density &density2 = s2_.getDensity();

	v.Four().Schur(density2.getDK(),w_att_.Four());
	v.Four().multBy(0.5*dV*dV/Ntot);
	v.do_fourier_2_real(); 
	s1_.addToForce(v.Real());

	E = 0.5*density1.getInteractionEnergy(v.Real());	
    }      
    return E;
  }

  double Mu(const vector<double> &x, int species) const
  {
    double mu = 0.0;

    if(s1_.getSequenceNumber() == species)
      mu += a_vdw_*x[s2_.getSequenceNumber()];

    if(s2_.getSequenceNumber() == species)
      mu += a_vdw_*x[s1_.getSequenceNumber()];    

    return mu;
  }

  double Fhelmholtz(const vector<double> &x) const
  {
    return a_vdw_*x[s1_.getSequenceNumber()]*x[s2_.getSequenceNumber()];
  }  

  double getVDWParameter() const { return a_vdw_;}
  
 private:
  Species &s1_;
  Species &s2_;

  double a_vdw_;
  DFT_FFT w_att_;  

};



#endif // __LUTSKO__INTERACTION__
