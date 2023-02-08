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


class String_Slave;
class String_Slave2;


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

  Density& operator=(const Density &a){copy(a);return *this;}

  void copy(const Density &ref)
  {
    Lattice::copy((Lattice&) ref);
    Density_ = ref.Density_;
  }
  
  // Initialize
  
  //set/get density 
  void set_from_smaller_density(const Density &density);

  void set(const char *filename){ readDensity(filename);} // read from file
  void set(const DFT_Vec &x) { Density_.set(x);}
  void set(double density) { Density_.set(density);} 

  void set(int i, int j, int k, double val)   { set(pos(i,j,k),val);}
  void set(unsigned pos, double val)  { Density_.set(pos,val);}  

  void set_background_density(double val);
  
  virtual double  get(long pos) const { return Density_.cReal().get(pos);}   
  virtual double  get(int ix, int iy, int iz) const {putIntoBox(ix,iy,iz); return Density_.cReal().get(pos(ix,iy,iz));}
  complex<double> get_fourier_value(long pos) const {return Density_.cFour().get(pos);}
  
  //Decendents of the Density object can implement this method to do graphical displays
  virtual void doDisplay(string &title, string &file, int seq = 0, void *param = NULL) const {}
  virtual void turn_on_display(){}
  
  //Some simple properties
  double getNumberAtoms()      const;
  double get_number_of_atoms() const { return getNumberAtoms();}
  double get_mass()            const { return getNumberAtoms();}

  double get_min()  const { return Density_.cReal().min();}
  double get_max()  const { return Density_.cReal().max();}
  void   get_center_of_mass(double &rx, double &ry, double &rz) const;
  double get_msd() const;
  void   get_particles(double threshold, vector< vector<long> > &clusters); // detect particles
  double get_ave_background_density() const;
  double get_max_background_density() const;
  double get_min_background_density() const;
  double get_dot_with(DFT_Vec &v) const;
  // get neighboring values: the return value is the density at pos.
  double get_neighbor_values(long pos,double &xpx, double &xmx, double &xpy, double &xmy, double &xpz, double &xmz) const;  
  
  // do stuff to the density

  void   operator*=(double a)   { Density_.MultBy(a);}
  void   operator+=(DFT_Vec &v) { Density_.IncrementBy(v);} 
  void   shift(DFT_Vec &direction, double scale) { Density_.Real().IncrementBy_Scaled_Vector(direction, scale);}
  
  void   doFFT() {Density_.do_real_2_fourier();}
  
  // access arrays  
  const double*          get_density_pointer() { return Density_.Real().memptr();}
  const DFT_Vec&         get_density_real()    const { return Density_.cReal();}
  const DFT_Vec_Complex& get_density_fourier() const { return Density_.cFour();} 
  
  // Read/Write
  void write_VTK_File(string filename) const;

  void writeDensity(string filename) const; // write in binary format: obsolete
  void readDensity(const char *filename);  // read from binary file: obsolete

  friend ostream &operator<<(ostream &of, const Density &d) {of << static_cast<const Lattice &>(d) << d.Density_; return of;}  
  friend istream &operator>>(istream  &in, Density &d )     {in >> static_cast<Lattice &>(d) >> d.Density_; return in;}

  template<typename Archive> void serialize(Archive & ar, const unsigned int version)
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

  friend String_Slave;
  friend String_Slave2;
  
protected:
  DFT_FFT Density_;  // The arrays for the real and fourier components of the density
};

class External_Field
{
public:
  External_Field() {}
  ~External_Field() {}

  int get_species() const {return species_;}
  const DFT_Vec& get_field() const { return field_;}

  template<class Archive> void serialize(Archive &ar, const unsigned int version)
  {
    ar & field_;
    ar & species_;
  }  
protected:
  DFT_Vec field_;
  int     species_ = 1;
};


BOOST_CLASS_VERSION(Density, 1)
BOOST_CLASS_VERSION(External_Field, 1)

#endif // __LUTSKO__DENSITY_ARRAY__
