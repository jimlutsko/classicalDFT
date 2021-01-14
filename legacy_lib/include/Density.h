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

#include "Lattice.h"

// There is some sort of a conflict between Boost headers and <complex> so this needs to come last in order to compile.

#include <complex>

#ifdef USE_OMP
#include <omp.h>
#endif

static const double SMALL_VALUE = 1e-18;

#include "DFT_LinAlg.h"
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
  Density() : Lattice(), Density_(), vWall_(){}
  Density(const Density &dd)
    : Lattice(dd), Density_(), vWall_()
  {
    Density_.initialize(Nx_,Ny_,Nz_);      
    vWall_ = dd.vWall_;
  }
  ~Density(){}

  virtual size_t number_of_parameters() const { return Ntot_;}
  
  virtual void initialize_from_smaller_density(const Density &density);
  virtual void initialize_from_file(const char *filename);

  void set(int i, int j, int k, double val) { set(pos(i,j,k),val);}
  void set(unsigned i, double val)  { Density_.Real().set(i,val);}
  void set(const DFT_Vec &x) { Density_.Real().set(x);}

  virtual double get(long pos) const { return Density_.cReal().get(pos);} 
  virtual double get(int ix, int iy, int iz) const { putIntoBox(ix,iy,iz);return Density_.cReal().get(pos(ix,iy,iz));}

  void scale_to(double a) { Density_.Real().MultBy(a);} 

  //const DFT_Vec& getDensity() const { return Density_.cReal();}  

  
  /**
  *   @brief Decendents of the Density object can implement this method to do graphical displays
  *  
  *   @param  title  the text information to display (e.g. as title) with image
  *   @param  file the file to write to (if doing this)
  *   @return  none
  */  
  virtual void doDisplay(string &title, string &file, int seq = 0) const {}
  
  /**
  *   @brief  Total number of particles
  *  
  *   @param  none
  *   @return  Number of particles
  */  
  virtual double getNumberAtoms() const { return Density_.cReal().accu()*dV();}
  virtual double min() const { return Density_.cReal().min();}
  virtual double max() const { return Density_.cReal().max();}
  
  
  /**
  *   @brief  Center of mass
  *  
  *   @param  rx: cm in x direction
  *   @param  ry: cm in y direction
  *   @param  rz: cm in z direction
  *   @return  cm
  */  
  void getCenterOfMass(double &rx, double &ry, double &rz) const;

    /**
  *   @brief  <R^2>
  *  
  *   @return R2 (for computing Lindemann parameter of uniform solid)
  */  
  double getRsquared() const ;


  /**
   *   @brief  Read-only accessor for array holding fft of density;
   *  
   *   @return Density_.Four()
   */  
  const DFT_Vec_Complex & getDK() const { return Density_.cFour();} 

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

  /**
  *   @brief  Write the density into a vtk-readable file
  *  
  *   @param  filename: name of file (should end in .vtk).
  *   @return none
  */  
  void write_VTK_File(string &filename);


  //  DFT_FFT& getFullVector() { return Density_;}

  friend ostream &operator<<(ostream &of, const Density &d) {of << static_cast<const Lattice &>(d) << d.Density_ << d.vWall_; return of;}  
  friend istream &operator>>(istream  &in, Density &d )     {in >> static_cast<Lattice &>(d) >> d.Density_ >> d.vWall_; return in;}

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Lattice>(*this);
    ar & Density_;
    ar & vWall_;
  }

  
  /**
  *   @brief  Write the real-space array Density_ to a file in binary format. THIS IS A LEGACY FUNCTION THAT SHOULD BE REMOVED AT SOME POINT.
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
  *   @brief  Read the real-space array Density_ from a file in binary format.  THIS IS A LEGACY FUNCTION THAT SHOULD BE REMOVED AT SOME POINT.
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
  
protected:

  DFT_FFT Density_;  // The arrays for the real and fourier components of the density

  DFT_Vec vWall_;            ///< Array holding wall potential at each lattice point: i = pos(ix,iy,iz)

};

#endif // __LUTSKO__DENSITY_ARRAY__
