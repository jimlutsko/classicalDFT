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

#include "Density.h"
#include "visit_writer.h"

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

void Density::detectClusters(double threshold, vector< vector<long> > &clusters)
{
  /*
  // Copy the density
  DFT_Vec p(Density_.Real());

  double m = threshold;
  double i = -1;

  for(int j=0;j<p.size();j++)
    {
      if(p.get(j) > m) {m = p.get(j); i = j;}

      if(i < 0) break; // no more clusters

      vector<long> new_cluster;

      new_cluster.push_back(m);
      p.set(m,0.0);
      
      // check each neighbor of each staged point
      for(long j=0;j<stage.size();j++)
	{
	  long point = stage[j];
	  long ix,iy,iz;
	  cartesian(point,ix,iy,iz);

	  for(int jx=-1;jx<=1;jx++)
	    for(int jy=-1;jy<=1;jy++)
	      for(int jz=-1;jz<=1;jz++)
		{
		  if(jx == 0 && jy == 0 && jz == 0) continue;
		  long pos = get_PBC_Pos(ix+jx,iy+jy,iz+jz);
		  if(p.get(pos) > threshold) {new_cluster.push_back(pos); p.set(pos,0.0);}
		}
	}
      
      clusters.push_back(new_cluster);
    }
  */
  throw std::runtime_error("Density::detectClusters not implemented");

}

void Density::write_VTK_File(string &filename)
{
  // I don't understand why it has to be plus 1 ...
  int dims[] = {int(Nx_+1), int(Ny_+1), int(Nz_+1)};

  int nvars = 1;
  int vardims[] = {1};
  int centering[] = {0};
  const char *varnames[] = {"density"};

  unsigned Nmax = Ntot();
  
  float *density = new float[Nmax];
  float *vars[] = {(float *)density};

  /* Create zonal variable */
  unsigned pos = 0;
  for(int k = 0; k < Nz_; ++k)
    for(int j = 0; j < Ny_; ++j)
      for(int i = 0; i < Nx_; ++i)
	{
	  density[pos] = getDensity(i,j,k);
	  pos++;
	}
  /* Use visit_writer to write a regular mesh with data. */
  write_regular_mesh(filename.c_str(), 0, dims, nvars, vardims,
		     centering, varnames, vars);
  delete density;
}
