#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

#include <gsl/gsl_integration.h>

#include <boost/serialization/export.hpp>

using namespace std;

#include "Density.h"
#include "External_Field.h"
#include "visit_writer.h"

BOOST_CLASS_EXPORT(Density)

Density::Density(double dx, double L[])
  : Lattice(dx, L), Density_(){ Density_.initialize(Nx_,Ny_,Nz_); }  // Allows child classes to do their own initialization
Density::Density(double dx, double dy, double dz, double L[])
  : Lattice(dx, dy, dz, L), Density_() {Density_.initialize(Nx_,Ny_,Nz_); }  // Allows child classes to do their own initialization
Density::Density(): Lattice(), Density_(){}
Density::Density(const Density &dd)
  : Lattice(dd), Density_(dd.Density_){}


// We are assuming that the lattice spacing is the same ...
// Assume that the input density has lattice points ix=0,...,Nx1-1
// and that the new system has lattice points ix=0,...,Nx2-1
// and assume that Nx2 > Nx1
// I will demand that Nx2-Nx1 is even
// and the idea is that we just pad the extra terms
void Density::set_from_smaller_density(const Density &density)
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

	      double d  = density.get(jx,jy,jz);
	      
	      set(ix,iy,iz,d);
	    }
}


// Here we also assume equal lattice spacings and demand that Nx1-Nx2 is even.
void Density::crop_from_larger_density(const Density &density)
{
  int Nx1 = density.Nx();
  int Ny1 = density.Ny();
  int Nz1 = density.Nz();
  
  int dNx = Nx1 - Nx_;
  int dNy = Ny1 - Ny_;
  int dNz = Nz1 - Nz_;

  int Mx = dNx/2;
  int My = dNy/2;
  int Mz = dNz/2;

  if(dNx != 2*Mx) cout << "Warning: difference in density lattices is not even in x direction" << endl;
  if(dNy != 2*My) cout << "Warning: difference in density lattices is not even in y direction" << endl;
  if(dNz != 2*Mz) cout << "Warning: difference in density lattices is not even in z direction" << endl;
  
  for(int ix=0;ix<Nx_;ix++)
  for(int iy=0;iy<Ny_;iy++)
  for(int iz=0;iz<Nz_;iz++)
  {
    double d = density.get(ix+Mx,iy+My,iz+Mz);
    set(ix,iy,iz,d);
  }
}


// Here I assume the dimensions of the box (Lx, Ly, Lz) are the same.
// I also assume that the lengths Lx, Ly, Lz are multiples of both spacings.
// The density is a trilinear interpolation of the nearest points on the coarse lattice.
void Density::set_from_coarser_density(const Density &density)
{
  double dx1 = density.getDX();
  double dy1 = density.getDY();
  double dz1 = density.getDZ();
  
  double dV1 = dx1*dy1*dz1;
  double Natoms1 = density.get_number_of_atoms();
  
  int Nx1 = int( Lx()/dx1 + 0.5 );
  int Ny1 = int( Ly()/dy1 + 0.5 );
  int Nz1 = int( Lz()/dz1 + 0.5 );
  
  cout << endl;
  cout << "Fine   spacings: " << setw(12) << dx_ << setw(12) << dy_ << setw(12) << dz_ << endl;
  cout << "Coarse spacings: " << setw(12) << dx1 << setw(12) << dy1 << setw(12) << dz1 << endl;
  
  cout << endl;
  cout << "Box dimensions:  " << setw(12) << Lx() << setw(12) << Ly() << setw(12) << Lz() << endl;
  cout << "Coarse Nx1*dx1:  " << setw(12) << Nx1*dx1 << setw(12) << Ny1*dy1 << setw(12) << Nz1*dz1 << endl;
  
  if (fabs(Lx()-Nx1*dx1)>1e-8*dx1) throw runtime_error("Error: the box length is not a multiple of the lattice spacing in the x direction");
  if (fabs(Ly()-Ny1*dy1)>1e-8*dy1) throw runtime_error("Error: the box length is not a multiple of the lattice spacing in the y direction");
  if (fabs(Lz()-Nz1*dz1)>1e-8*dz1) throw runtime_error("Error: the box length is not a multiple of the lattice spacing in the z direction");
  
  for (int ix=0; ix<Nx_; ix++)
  for (int iy=0; iy<Ny_; iy++)
  for (int iz=0; iz<Nz_; iz++)
  {
    double x = getX(ix);
    double y = getY(iy);
    double z = getZ(iz);
    
    int ix1 = int( (x+Lx()/2+0.01*dx_)/dx1 ); // +0.01*dx_ makes sure it is rounded correctly
    int iy1 = int( (y+Ly()/2+0.01*dx_)/dy1 );
    int iz1 = int( (z+Lz()/2+0.01*dx_)/dz1 );
    
    double x1 = density.getX(ix1);
    double y1 = density.getY(iy1);
    double z1 = density.getZ(iz1);
    
    double lx = (x-x1)/dx1; while (lx>=1.0) lx-=1.0; while (lx<0.0) lx+=1.0;
    double ly = (y-y1)/dy1; while (ly>=1.0) ly-=1.0; while (ly<0.0) ly+=1.0;
    double lz = (z-z1)/dz1; while (lz>=1.0) lz-=1.0; while (lz<0.0) lz+=1.0;
    
    double d000 = density.get(ix1+0,iy1+0,iz1+0);
    double d100 = density.get(ix1+1,iy1+0,iz1+0);
    double d010 = density.get(ix1+0,iy1+1,iz1+0);
    double d001 = density.get(ix1+0,iy1+0,iz1+1);
    double d011 = density.get(ix1+0,iy1+1,iz1+1);
    double d101 = density.get(ix1+1,iy1+0,iz1+1);
    double d110 = density.get(ix1+1,iy1+1,iz1+0);
    double d111 = density.get(ix1+1,iy1+1,iz1+1);
    
    double d  = (1-lx)*(1-ly)*(1-lz)*d000;
           d +=    lx *(1-ly)*(1-lz)*d100;
           d += (1-lx)*   ly *(1-lz)*d010;
           d += (1-lx)*(1-ly)*   lz *d001;
           d += (1-lx)*   ly *   lz *d011;
           d +=    lx *(1-ly)*   lz *d101;
           d +=    lx *   ly *(1-lz)*d110;
           d +=    lx *   ly *   lz *d111;
    
    set(ix,iy,iz,d);
    
    // Check
    /*
    if (ix==Nx_-1 && iy==Ny_/2 && iz==Nz_/2+1)
    {
      cout << endl;
      cout << "Fine   ix, iy, iz: " << setw(12) << ix  << setw(12) << iy  << setw(12) << iz  << endl;
      cout << "Coarse ix, iy, iz: " << setw(12) << ix1 << setw(12) << iy1 << setw(12) << iz1 << endl;
      
      cout << endl;
      cout << "Fine   x, y, z: " << setw(12) << x  << setw(12) << y  << setw(12) << z  << endl;
      cout << "Coarse x, y, z: " << setw(12) << x1 << setw(12) << y1 << setw(12) << z1 << endl;
      
      cout << endl;
      cout << "lx, ly, lz: " << setw(12) << lx << setw(12) << ly << setw(12) << lz << endl;
      
      cout << endl;
      cout << "d000: " << setw(12) << d000 << endl;
      cout << "d100: " << setw(12) << d100 << endl;
      cout << "d010: " << setw(12) << d010 << endl;
      cout << "d001: " << setw(12) << d001 << endl;
      cout << "d011: " << setw(12) << d011 << endl;
      cout << "d101: " << setw(12) << d101 << endl;
      cout << "d110: " << setw(12) << d110 << endl;
      cout << "d111: " << setw(12) << d111 << endl;
      
      
      cout << endl;
      cout << "Trilinear interpolation d : " << setw(12) << d << endl;
    }
    */
  }
  
  cout << endl;
  cout << "Number of particles in coarse density: " << setw(12) << Natoms1 << endl;
  cout << "Number of particles in fine   density: " << setw(12) << get_number_of_atoms() << endl;
}


void Density::get_particles(double threshold, vector< vector<long> > &clusters)
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

double Density::getNumberAtoms() const //{ return Density_.cReal().accu()*dV();}
{
  Summation s;
  for(long i=0;i<Ntot(); i++)
    s.add(Density_.cReal().get(i));
  return s*dV();  
}

double Density::get_dot_with(DFT_Vec &v) const
{
  Summation s;
  for(long i=0;i<Ntot(); i++)
    s.add(Density_.cReal().get(i)*v.get(i));
  return s;
}

  // Center of mass
void Density::get_center_of_mass(double &rx, double &ry, double &rz) const
{
  rx = ry = rz = 0.0;
  double m = 0;
  for(int i=0;i<Nx_;i++)
    for(int j=0;j<Ny_;j++)
      for(int k=0;k<Nz_;k++)
	{
	  double d = get(i,j,k);
	  rx += d*i;
	  ry += d*j;
	  rz += d*k;
	  m += d;
	}
  rx /= m;
  ry /= m;
  rz /= m;
  return;
}

  // <R^2>
double Density::get_msd() const 
{
  double r2 = 0;
  double m = 0;
  for(int i=0;i<Nx_;i++)
    for(int j=0;j<Ny_;j++)
      for(int k=0;k<Nz_;k++)
	{
	  double x = getX(i);
	  double y = getY(j);
	  double z = getZ(k);

	  double d = get(i,j,k);
	  if(x*x+y*y+z*z < (0.5*Nx_*dx_)*(0.5*Nx_*dx_)) //+(0.5*Ny_*dy_)*(0.5*Ny_*dy_)+(0.5*Nz_*dz_)*(0.5*Nz_*dz_))
	    {
	      r2 += d*(x*x+y*y+z*z);
	      m += d;
	    }
	}
  return r2/m;
}


void Density::write_VTK_File(string filename) const
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
	  density[pos] = get(i,j,k);
	  pos++;
	}
  /* Use visit_writer to write a regular mesh with data. */
  write_regular_mesh(filename.c_str(), 0, dims, nvars, vardims,
		     centering, varnames, vars);
  delete density;
}

void Density::writeDensity(string filename) const
{
  ofstream of(filename.c_str(),ios::binary);
  Density_.cReal().save(of);
}

void Density::readDensity(const char *filename) // read from binary file: obsolete
{
  ifstream in(filename,ios::binary);
  if(! in.good())
    {
      stringstream s;
      s << "Could not open input file " << filename << " ... aborting ... " << endl;
      throw std::runtime_error(s.str());
    }
  // everything OK so read it. 
  Density_.Real().load(in);
}

double Density::get_ave_background_density() const
{
  double d = 0;
  double n = 0;
  

  for(int ix=0;ix<Nx_;ix++)
    for(int iy=0;iy<Ny_;iy++)
      {d += get(ix,iy,0); n++;}

  for(int ix=0;ix<Nx_;ix++)
    for(int iz=0;iz<Nz_;iz++)
      {d += get(ix,0,iz); n++;}

  for(int iy=0;iy<Ny_;iy++)
    for(int iz=0;iz<Nz_;iz++)    
      {d += get(0,iy,iz); n++;}  


  return d/n;
}

double Density::get_max_background_density() const
{
  double d = get(0,0,0);
  
  for(int ix=0;ix<Nx_;ix++)
    for(int iy=0;iy<Ny_;iy++)
      d = std::max(d,get(ix,iy,0));

  for(int ix=0;ix<Nx_;ix++)
    for(int iz=0;iz<Nz_;iz++)
      d = std::max(d,get(ix,0,iz));

  for(int iy=0;iy<Ny_;iy++)
    for(int iz=0;iz<Nz_;iz++)
      d = std::max(d,get(0,iy,iz));


  return d;
}

double Density::get_min_background_density() const
{
  double d = get(0,0,0);
  
  for(int ix=0;ix<Nx_;ix++)
    for(int iy=0;iy<Ny_;iy++)
      d = std::min(d,get(ix,iy,0));

  for(int ix=0;ix<Nx_;ix++)
    for(int iz=0;iz<Nz_;iz++)
      d = std::min(d,get(ix,0,iz));

  for(int iy=0;iy<Ny_;iy++)
    for(int iz=0;iz<Nz_;iz++)
      d = std::min(d,get(0,iy,iz));


  return d;
}

void Density::set_background_density(double val)
{
  for(int ix=0;ix<Nx_;ix++)
    for(int iy=0;iy<Ny_;iy++)
      set(ix,iy,0, val);

  for(int ix=0;ix<Nx_;ix++)
    for(int iz=0;iz<Nz_;iz++)
      set(ix,0,iz,val); 

  for(int iy=0;iy<Ny_;iy++)
    for(int iz=0;iz<Nz_;iz++)    
      set(0,iy,iz, val);
}



// No better place to stash this, for the moment. 
void Lattice::test_boundary_coding()
{
  long pos = -1;

  int ix,iy,iz;
  iz = 0;
  cout << Nx_ << " " << Ny_ << " " << Nz_ << endl;
  for(ix = 0;ix < Nx_; ix++)
    for(iy = 0; iy < Ny_; iy++)
      {
	pos++;

	long p = boundary_pos(ix,iy,iz);
	cout << ix << " " << iy << " " << iz << " " << pos << " " << p << endl;	
	if(p != pos) throw std::runtime_error("p1");

	int jx,jy,jz; jx = jy = jz = -1;
	boundary_cartesian(pos,jx,jy,jz);
	cout << ix << " " << iy << " " << iz << " " << pos << " " << jx << " " << jy << " " << jz << endl;
	if(ix != ix || jy != iy || jz != iz) throw std::runtime_error("p2");
      }

  iy = 0;
  for(ix = 0;ix < Nx_; ix++)
    for(iz = 1; iz < Nz_; iz++)
      {
	pos++;

	long p = boundary_pos(ix,iy,iz);
	cout << ix << " " << iy << " " << iz << " " << pos << " " << p << endl;
	if(p != pos) throw std::runtime_error("p3");

	int jx,jy,jz; jx = jy = jz = -1;
	boundary_cartesian(pos,jx,jy,jz);
	cout << ix << " " << iy << " " << iz << " " << pos << " " << jx << " " << jy << " " << jz << endl;	
	if(ix != ix || jy != iy || jz != iz) throw std::runtime_error("p4");
      }

  ix = 0;
  for(iy = 1;iy < Ny_; iy++)
    for(iz = 1; iz < Nz_; iz++)
      {
	pos++;

	long p = boundary_pos(ix,iy,iz);
	cout << ix << " " << iy << " " << iz << " " << pos << " " << p << endl;	
	if(p != pos) throw std::runtime_error("p5");

	int jx,jy,jz; jx = jy = jz = -1;
	boundary_cartesian(pos,jx,jy,jz);
	cout << ix << " " << iy << " " << iz << " " << pos << " " << jx << " " << jy << " " << jz << endl;	
	if(ix != ix || jy != iy || jz != iz) throw std::runtime_error("p6");
      }     
  cout << "pos = " << pos << endl;
  cout << "Test passed";
}


double Density::get_neighbor_values(long pos,double &xpx, double &xmx, double &xpy, double &xmy, double &xpz, double &xmz) const
{
  int ix, iy,iz;
  
  cartesian(pos,ix,iy,iz);

  xpx =  get(ix+1, iy,   iz);
  xmx =  get(ix-1, iy,   iz);
  xpy =  get(ix,   iy+1, iz);
  xmy =  get(ix,   iy-1, iz);
  xpz =  get(ix,   iy,   iz+1);
  xmz =  get(ix,   iy,   iz-1);
  return get(ix,   iy,   iz);
}

//template<typename Archive> void Density::serialize(Archive & ar, const unsigned int version)
/*
template<class Archive> void Density::serialize(Archive & ar, const unsigned int version)
{
  ar & boost::serialization::base_object<Lattice>(*this);
  ar & Density_;
  if(version == 0) // this is only for reading, just to get us to the right place in the input file if it is from the previous version ...
    {
      DFT_Vec vWall_;
      vWall_.zeros(Ntot_);
      ar & vWall_;
    }
}
*/


/*
template<typename Archive> void External_Field::serialize(Archive & ar, const unsigned int version)
{
  ar & boost::serialization::base_object<Lattice>(*this);
  ar & species_;
  //  ar & field_;
}
*/


template<class Archive>
void boost::serialization::save_construct_data(Archive & ar, const External_Field * t, const unsigned int file_version)
{
  ar << t->field_;
  ar << t->species_;

}

template<class Archive>
void boost::serialization::load_construct_data(Archive & ar, External_Field * t, const unsigned int file_version)
{
  ar >> t->field_;
  ar >> t->species_;
}
 
