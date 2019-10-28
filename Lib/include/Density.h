/* This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
 * To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
 *
 * Author: James F. Lutsko
 * www.lutsko.com
 */

#ifndef __LUTSKO__DENSITY_ARRAY__
#define __LUTSKO__DENSITY_ARRAY__

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <complex>
#include <complex.h>

static const double SMALL_VALUE = 1e-18;

#include "DFT_LinAlg.h"
#include "Lattice.h"
#include "Summation.h"



/**
  *  @brief Base class for the density: basically, a wrapper for the array holding the density at each lattice site
  */  

class Density : public Lattice
{
 public:
  /**
  *   @brief  Default  constructor for Density
  *  
  *   @param  dx is lattice spacing: assumed to be the same in all directions
  *   @param  L[] are the dimensions of the physical region
  *   @return nothing 
  */  
 Density(double dx, double L[])
   : Lattice(dx, L), Density_(), vWall_()
    {
      Density_.initialize(Nx_,Ny_,Nz_);  // Allows child classes to do their own initialization
      vWall_.zeros(Ntot_);
    }

  /**
  *   @brief  Copy constructor
  *  
  *   @param  dd Density object to copy
  *   @return nothing 
  */    
  Density(const Density &dd)
    : Lattice(dd), Density_(), vWall_()
    {
      Density_.initialize(Nx_,Ny_,Nz_);      
      vWall_ = dd.vWall_;
    }
  
  /**
  *   @brief  Destructor for Density
  *  
  *   @return nothing 
  */  
  ~Density(){}

  /**
  *   @brief Called to read initial density from file. This will fail if one attempts to read a density object of a different size. 
  *  
  *   @param  filename
  *   @return  none
  */  
  void initialize_from_file(const char *filename);

    /**
  *   @brief Called to create initial density from a smaller geometry.
  *  
  *   @param  filename
  *   @return  none
  */  
  void initialize_from_smaller_density(const Density &density);

  /**
  *   @brief Decendents of the Density object can implement this method to do graphical displays
  *  
  *   @param  title  the text information to display (e.g. as title) with image
  *   @param  file the file to write to (if doing this)
  *   @return  none
  */  
  virtual void doDisplay(string &title, string &file, int seq = 0) const = 0; //{}

  /**
  *   @brief  Total number of particles
  *  
  *   @param  none
  *   @return  Number of particles
  */  
  virtual double getNumberAtoms() const //{ return Density_.cReal().accu()*dV();}
  {

    double dV1 = dV();
    Summation s;
    for(long i=0;i<Ntot(); i++)
      s.add(Density_.cReal().get(i)*dV1);
    return s;

  }

  
  /**
  *   @brief  Center of mass
  *  
  *   @param  rx: cm in x direction
  *   @param  ry: cm in y direction
  *   @param  rz: cm in z direction
  *   @return  cm
  */  
  void getCenterOfMass(double &rx, double &ry, double &rz) const
  {
    rx = ry = rz = 0.0;
    double m = 0;
    for(int i=0;i<Nx_;i++)
      for(int j=0;j<Ny_;j++)
	for(int k=0;k<Nz_;k++)
	  {
	    double d = getDensity(i,j,k);
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


  /**
  *   @brief  Write the real-space array Density_ to a file in binary format.
  *  
  *   @param  filename: name of the file to write to.
  *   @return  
  */  
  void writeDensity(string &filename) const
  {
    ofstream of(filename.c_str(),ios::binary);
    Density_.cReal().save(of);
  }

  /**
  *   @brief  Read the real-space array Density_ from a file in binary format.
  *  
  *   @param  filename: name of the file to read from.
  *   @return  
  */  
  void readDensity(const char *filename)
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

  /**
  *   @brief  set density at a given lattice position
  *  
  *   @param  i: lattice position in x direction
  *   @param  j: lattice position in y direction
  *   @param  k: lattice position in z direction
  *   @param  val: density of the given point
  *   @return none
  */  
  void set_Density_Elem(int i, int j, int k, double val)   { set_Density_Elem(pos(i,j,k),val);}
    //  { Density_.Real().set(pos(i,j,k),val);}

  /**
  *   @brief  set density at a given lattice position
  *  
  *   @param  i: multi-index lattice positionlattice position in x direction
  *   @param  val: density of the given point
  *   @return none
  */    
  void set_Density_Elem(unsigned i, double val)  { Density_.Real().set(i,val);}

  /**
  *   @brief  set density at all positions based on amplitude: d = SMALL_VALUE + amplitude*amplitude: used during minimization to assure a positive definite density.
  *  
  *   @param  x: amplitude
  *   @return none
  */  
  void set_density_from_amplitude(const DFT_Vec &x) { Density_.Real().setFromAmplitude(SMALL_VALUE,x);}

  /**
  *   @brief  set density by copying given vector
  *  
  *   @param  x: density to copy
  *   @return none
  */  
  void set(const DFT_Vec &x) { Density_.Real().set(x);}

  /**
  *   @brief  scale the density by a scalar
  *  
  *   @param  a: all densities are multiplied by a.
  *   @return none
  */  
  void scale_to(double a) { Density_.Real().MultBy(a);} 

  /**
  *   @brief  Get density at point i
  *  
  *   @param  i: index of element to be returned
  *   @return Density_.Real[i]
  */  
  virtual double getDensity(long i) const { return Density_.cReal().get(i);} 

  /**
  *   @brief  Get density at coordinates ix,iy,iz
  *  
  *   @param  ix: index of point in x-direction
  *   @param  iy: index of point in y-direction
  *   @param  iz: index of point in z-direction
  *   @return density of given point
  */  
  virtual double getDensity(int ix, int iy, int iz) const
  { 
    putIntoBox(ix,iy,iz);
    return Density_.cReal().get(pos(ix,iy,iz));
  }

  /**
  *   @brief  Translate coordinates ix,iy,iz into an index taking account of the periodic boundaries
  *  
  *   @param  ix: index of point in x-direction
  *   @param  iy: index of point in y-direction
  *   @param  iz: index of point in z-direction
  *   @return index
  */  
  virtual long get_PBC_Pos(int ix, int iy, int iz) const
  { 
    putIntoBox(ix,iy,iz);
    return pos(ix,iy,iz);
  }

  /**
   *   @brief  Read-only accessor for entire real-space density array
   *  
   *   @return Density_.Real()
   */  
  const DFT_Vec& getDensity() const { return Density_.cReal();}

  /**
   *   @brief  Directly access the data
   *  
   *   @return double* array
   */  
  double* getData() { return Density_.Real().memptr();}

  /**
   *   @brief  Read-only accessor for array holding fft of density;
   *  
   *   @return Density_.Four()
   */  
  const DFT_Vec_Complex & getDK() const { return Density_.cFour();} 

  /**
  *   @brief  Density fft at wavevector i
  *  
  *   @return Density_.Four()[i]
  */  
  const complex<double> DensityK(long i) const { return Density_.cFour().get(i);}

  /**
  *   @brief  Accessor for array holding the wall potential
  *  
  *   @return vWall_
  */  
  const DFT_Vec& getExternalField() const { return vWall_;}

  /**
  *   @brief  Perform fft of density (from real to fourier space)
  *  
  *   @return none
  */  
  void doFFT() {Density_.do_real_2_fourier();} 

  /**
   *  @brief  Returns true for lattice points for which the external potential is infinite and the density is therefore zero:
   *  
   *   @return boolean true/false
   */  
  virtual bool IsInUnphysicalRegion(int i, int j, int k) const {return false;}

  /**
  *   @brief  Get the value of the external (wall) field at point (ijk)
  *  
  *   @param  ijk are the positions
  *   @return  
  */  
  virtual double getField(int i, int j, int k) const { return vWall_.get(pos(i,j,k));}

  /**
   *  @brief  Calculates the wall contribution to the energy (but without factors of dV)
   *  
   *   @return returns sum_i Density_.Real()[i] vwall_[i]
   */  
  double getExternalFieldEnergy() const { return vWall_.dotWith(Density_.cReal());} 

  /**
   *  @brief  Calculates the mean field contribution to the energy (but without factors of dV)
   *
   * @param vdr is meant to be the integral of the interaction potential, v, and the density (in real space)
   *  
   *   @return returns sum_i Density_.Real()[i] vdr[i]
   */  
  double getInteractionEnergy(DFT_Vec &vdr) const
  {
    Summation s;
    for(long i=0;i<Ntot(); i++)
      s.add(Density_.cReal().get(i)*vdr.get(i));
    return s;
    
    //    return Density_.cReal().dotWith(vdr);
  }

  /**
  *   @brief  Get the value of the derivative of the external field in a given direction. This is needed for DDFT calculations.
  *  
  *   @param  ix is the index of the desired point in the x-direction.
  *   @param  iy is the index of the desired point in the y-direction.
  *   @param  iz is the index of the desired point in the z-direction.
  *   @param  direction can be 0,1 or 2 for the direction of the derivative.
  *   @return Because of the needs of the DDFT code, this actually returns the derivative evaluated at iz-0.5.   
  */  
  virtual double getFieldDeriv(int ix, int iy, int iz, int direction) const { return 0.0;}

  /**
  *   @brief  Detects clusters: e.g. particles
  *  
  *   @param  threshold is min density - all points with density below this are treated as noise
  *   @param  vector of clusters: each cluster is a vector of indices corresponding to density poitns
  *   @return none
  */  
  void detectClusters(double threshold, vector< vector<long> > &clusters);

  
 protected:

  DFT_FFT Density_;  // The arrays for the real and fourier components of the density

  DFT_Vec vWall_;            ///< Array holding wall potential at each lattice point: i = pos(ix,iy,iz)

 private:

};

#endif // __LUTSKO__DENSITY_ARRAY__
