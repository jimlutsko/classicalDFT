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

#include "Density.h"

void Density::initialize_from_file(const char *filename)
{
  readDensity(filename);
}

// We are assuming that the lattice spacing is the same ...
// Assume that the input density has lattice points ix=0,...,Nx1-1
// and that the new system has lattice points ix=0,...,Nx2-1
// and assume that Nx2 > Nx1
// I will demand that Nx2-Nx1 is even
// and the idea is that we just pad the extra terms
void Density::initialize_from_smaller_density(const Density &density)
{

  long Nx1 = density.Nx();
  long Ny1 = density.Ny();
  long Nz1 = density.Nz();
  
  long dNx = Nx_ - Nx1;
  long dNy = Ny_ - Ny1;
  long dNz = Nz_ - Nz1;

  long Mx = dNx/2;
  long My = dNy/2;
  long Mz = dNz/2;

  if(dNx != 2*Mx) cout << "Warning: difference in density lattices is not even in x direction" << endl;
  if(dNy != 2*My) cout << "Warning: difference in density lattices is not even in y direction" << endl;
  if(dNz != 2*Mz) cout << "Warning: difference in density lattices is not even in z direction" << endl;

  
  for(int ix=0;ix<Nx_;ix++)
      for(int iy=0;iy<Ny_;iy++)
	  for(int iz=0;iz<Nz_;iz++)
	    {
	      int jx = ix - Mx;
	      if(jx < 0) jx = 0;
	      else if(jx > Nx1-1) jx = Nx1-1;

	      int jy = iy - My;
	      if(jy < 0) jy = 0;
	      else if(jy > Ny1-1) jy = Ny1-1;

	      int jz = iz - Mz;
	      if(jz < 0) jz = 0;
	      else if(jz > Nz1-1) jz = Nz1-1; 

	      double d  = density.getDensity(jx,jy,jz);
	      
	      set_Density_Elem(ix,iy,iz,d);
	    }
}
