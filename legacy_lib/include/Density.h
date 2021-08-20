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


static const double SMALL_VALUE = 1e-18;

#include "DFT_LinAlg.h"
#include "Summation.h"



/**
  *  @brief Base class for the density: basically, a wrapper for the array holding the density at each lattice site
  */  

class Density : public Lattice
{
 public:

  //Constructors/Destructors
  Density(double dx, double L[]);
  Density(double dx, double dy, double dz, double L[]);
  Density();
  Density(const Density &dd);
  
  ~Density(){}

  // Initialize
  void initialize_from_file(const char *filename);
  void initialize_from_smaller_density(const Density &density);

  //set/get density 
  void set(int i, int j, int k, double val)   { set(pos(i,j,k),val);}
  void set(unsigned pos, double val)  { Density_.Real().set(pos,val);}  
  void set(const DFT_Vec &x) { Density_.Real().set(x);}
  void set_background_density(double val);

  virtual double  get(long pos) const { return Density_.cReal().get(pos);}   
  virtual double  get(int ix, int iy, int iz) const {putIntoBox(ix,iy,iz); return Density_.cReal().get(pos(ix,iy,iz));}
  complex<double> get_fourier_value(long pos) const {return Density_.cFour().get(pos);}
  
  //Decendents of the Density object can implement this method to do graphical displays
  virtual void doDisplay(string &title, string &file, int seq = 0) const {}

  //Some simple properties
  double getNumberAtoms() const;
  double get_number_of_atoms() const { return Density_.cReal().accu()*dV();}  
  double get_min() const { return Density_.cReal().min();}
  double get_max() const { return Density_.cReal().max();}
  void   get_center_of_mass(double &rx, double &ry, double &rz) const;
  double get_msd() const;
  void   get_particles(double threshold, vector< vector<long> > &clusters); // detect particles
  double get_ave_background_density() const;
  double get_max_background_density() const;
  double get_min_background_density() const;
  double get_dot_with(DFT_Vec &v) const;
  
  // do stuff to the density
  void   scale_to(double a) { Density_.Real().MultBy(a);}
  void   shift(DFT_Vec &v)  { Density_.Real().IncrementBy(v);} 
  void   doFFT()            {Density_.do_real_2_fourier();}

  // access arrays  
  const double*          get_density_pointer() { return Density_.Real().memptr();}

  const DFT_Vec&         get_density_real()    const { return Density_.cReal();}
  const DFT_Vec_Complex& get_density_fourier() const { return Density_.cFour();} 
  const DFT_Vec&         get_external_field()  const { return vWall_;}

  // external field (wall_)
  virtual double get_field(int i, int j, int k) const { return vWall_.get(pos(i,j,k));}
  virtual double get_field_deriv(int ix, int iy, int iz, int direction) const { return 0.0;}  
  double         get_field_dot_density() const { return vWall_.dotWith(Density_.cReal());}   
  
  // Read/Write
  void write_VTK_File(string &filename);
  void writeDensity(string &filename) const; // write in binary format: obsolete
  void readDensity(const char *filename);  // read from binary file: obsolete

  friend ostream &operator<<(ostream &of, const Density &d) {of << static_cast<const Lattice &>(d) << d.Density_ << d.vWall_; return of;}  
  friend istream &operator>>(istream  &in, Density &d )     {in >> static_cast<Lattice &>(d) >> d.Density_ >> d.vWall_; return in;}
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Lattice>(*this);
    ar & Density_;
    ar & vWall_;
  }
protected:

  DFT_FFT Density_;  // The arrays for the real and fourier components of the density
  DFT_Vec vWall_;            ///< Array holding wall potential at each lattice point: i = pos(ix,iy,iz)

};

#endif // __LUTSKO__DENSITY_ARRAY__
