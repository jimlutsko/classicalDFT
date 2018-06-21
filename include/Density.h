#ifndef __LUTSKO__DENSITY_ARRAY__
#define __LUTSKO__DENSITY_ARRAY__

//#include "spherical.h"

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <complex>
#include <complex.h>

#include <mgl2/mgl.h>
//#include <mgl2/fltk.h>

static const double SMALL_VALUE = 1e-18;

#include "DFT_LinAlg.h"

/**
  *  @brief This class encapsulates all information relating to the lattice: the number of points, the lattice spaceing and the length of the cell in each direction. Allowance is made for different spacings in each direction but current implementations force all spacings to be equal.
  */  

class Lattice
{
 public:
  /**
  *   @brief  Default  constructor for Lattice
  *  
  *   @param  dx is lattice spacing: assumed to be the same in all directions
  *   @param  L[] are the dimensions of the physical region
  *   @return nothing 
  */  

  Lattice(double dx, double L[])
    : Nx_(0), Ny_(0), Nz_(0), dx_(dx), dy_(dx), dz_(dx)
    {
      L_[0] =L[0];
      L_[1] =L[1];
      L_[2] =L[2];

      Nx_ = long((L[0]+1e-8*dx_)/dx_);
      Ny_ = long((L[1]+1e-8*dy_)/dy_);
      Nz_ = long((L[2]+1e-8*dz_)/dz_);

      L_[0] = (Nx_-0)*dx_;
      L_[1] = (Ny_-0)*dy_;
      L_[2] = (Nz_-0)*dz_;


      if(fabs(L[0]-Nx_*dx_) > 1e-5*dx_)
	throw std::runtime_error("Box length incommensurate with lattice in x direction");

      if(fabs(L[1]-Ny_*dy_) > 1e-5*dy_)
	throw std::runtime_error("Box length incommensurate with lattice in y direction");

      if(fabs(L[2]-Nz_*dz_) > 1e-5*dz_)
	throw std::runtime_error("Box length incommensurate with lattice in z direction");

      Ntot_ = Nx_*Ny_*Nz_;
      Nout_ = Nx_*Ny_*((Nz_/2)+1);
    }
  /**
  *   @brief  Copy constructor for Lattice
  *   @param  ref: the Lattice to be copied
  *   @return nothing 
  */  

  Lattice(const Lattice& ref)
    {
      Nx_ = ref.Nx_;
      Ny_ = ref.Ny_;
      Nz_ = ref.Nz_;

      Ntot_ = ref.Ntot_;
      Nout_ = ref.Nout_;

      dx_ = ref.dx_;
      dy_ = ref.dy_;
      dz_ = ref.dz_;

      L_[0] = ref.L_[0];
      L_[1] = ref.L_[1];
      L_[2] = ref.L_[2];
    }

  /**
  *   @brief  Translate a (Cartesian) x-index into a position in the cubic box.
  *  
  *   @param  x-index
  *   @return position x
  */  
  double getX(int i) const { return dx_*(i-(Nx_-1)/2);}
  /**
  *   @brief  Translate a (Cartesian) y-index into a position in the cubic box.
  *  
  *   @param  y-index
  *   @return position y
  */  
  double getY(int i) const { return dy_*(i-(Ny_-1)/2);}
  /**
  *   @brief  Translate a (Cartesian) z-index into a position in the cubic box.
  *  
  *   @param  z-index
  *   @return position z
  */  
  double getZ(int i) const { return dz_*(i-(Nz_-1)/2);}

  /**
  *   @brief  Accessor for spaceing in x-direction
  *  
  *   @return dx_
  */  
  double getDX() const {return dx_;}
  /**
  *   @brief  Accessor for spaceing in y-direction
  *  
  *   @return dy_
  */  
  double getDY() const {return dy_;}
  /**
  *   @brief  Accessor for spaceing in z-direction
  *  
  *   @return dz_
  */  
  double getDZ() const {return dz_;}

  /**
   *  @brief  Translated (ix,iy,iz) cartesian indices into single index
   *  
   *   @return i(ix,iy,iz)
   */  
  inline long pos(long i, long j, long k) const { return k + Nz_*(j +  Ny_*i);}

  /**
   *  @brief  Accessor for (total) number of lattice points in x direction
   *  
   *   @return returns Nx_
   */  
  long Nx() const { return Nx_;}
  /**
   *  @brief  Accessor for (total) number of lattice points in y direction
   *  
   *   @return returns Ny_
   */    
  long Ny() const { return Ny_;}
  /**
   *  @brief  Accessor for (total) number of lattice points in z direction
   *  
   *   @return returns Nz_
   */  
  long Nz() const { return Nz_;}

    /**
   *  @brief  Accessor for box length in x direction
   *  
   *   @return returns L_[0]
   */  
  double Lx() const { return L_[0];}
  /**
   *  @brief  Accessor for nox length in y direction
   *  
   *   @return returns L_[1]
   */    
  double Ly() const { return L_[1];}
  /**
   *  @brief  Accessor for box length in z direction
   *  
   *   @return returns L_[2]
   */  
  double Lz() const { return L_[2];}

  /**
   *  @brief  Accessor for (total) number of lattice points redundant with @size()
   *  
   *   @return returns Ntot_
   */  
  long Ntot() const { return Ntot_;}

  /**
   *  @brief  Accessor for (total) number of lattice points : redundant with @Ntot()
   *  
   *   @return returns Ntot_
   */  
  long size() const { return Ntot_;}

  /**
   *  @brief  Accessor for (total) number of (complex) wave vectors returned by fftw3
   *  
   *   @return returns Nx_*Ny_*(int(Nz_/2)+1)
   */  
  long Nout() const { return Nout_;}

  /**
   *  @brief  Accessor for volume of lattice cell
   *  
   *   @return returns dx_*dy_*dz_
   */  
  double dV() const { return dx_*dy_*dz_;}

  /**
   *  @brief  Accessor for total simulation volume
   *  
   *   @return returns L[0]*L[1]*L[2]
   */  
  double getVolume() const { return L_[0]*L_[1]*L_[2];}

  /**
  *   @brief  Apply periodic boundary conditions to the lattice coordinates
  *  
  *   @param  ix- index of point in x-direction
  *   @param  iy- index of point in y-direction
  *   @param  iz- index of point in z-direction
  *   @return none
  */  
  void putIntoBox(int &ix, int &iy, int &iz) const
  {
    if(ix < 0)  ix += Nx_; 
    if(ix >= Nx_) ix -= Nx_;
    
    if(iy < 0)  iy += Ny_; 
    if(iy >= Ny_) iy -= Ny_;
    
    if(iz < 0)  iz += Nz_; 
    if(iz >= Nz_) iz -= Nz_; 
  }

 protected:
  long Nx_; ///< Number of lattice points in x direction
  long Ny_; ///< Number of lattice points in y direction
  long Nz_; ///< Number of lattice points in z direction

  long Ntot_;  ///< Total number of lattice points
  long Nout_;  ///< Total number of lattice points in fourier space (a la fftw)

  double dx_;   ///< Lattice spacing in x direction
  double dy_;   ///< Lattice spacing in y direction
  double dz_;   ///< Lattice spacing in z direction

  double L_[3];    ///< Dimensions of physical box in hard-sphere units
};
 
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
      Density_.initialize(Nx_,Ny_,Nz_);
      
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
  *   @brief  Destructorfor Density
  *  
  *   @return nothing 
  */  
  ~Density(){}

  virtual void initialize(double density, double d2) {}
  
  /**
  *   @brief Called to read initial density from file. This will fail if one attempts to read a density object of a different size. 
  *  
  *   @param  filename
  *   @return  none
  */  
  void initialize_from_file(const char *filename);

    /**
  *   @brief Called to creqte initial density from a smaller geometry.
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
  virtual void doDisplay(string &title, string &file){}

  /**
  *   @brief  Total number of particles
  *  
  *   @param  none
  *   @return  Number of particles
  */  
  virtual double getNumberAtoms() const { return Density_.cReal().accu()*dV();}

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
  void set_Density_Elem(int i, int j, int k, double val)  { Density_.Real().set(pos(i,j,k),val);}
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
  void scale_to(double a) { Density_.Real().multBy(a);} 

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
  void doFFT() {Density_.do_real_2_fourier();} //fftw_execute(p_);}

  /**
   *   @brief  Allows derived classes to initialize an array for displaying two-d graphics since only they know the size and shape of array needed.
   *  
   *   @param  a: the array to be intitialized
   *   @return none
   */  
  virtual void initialize_2D_data( mglData& a) const {throw std::runtime_error("Cannot be here: Density::initialize_2D_data");}

  /**
   *   @brief  Allows derived classes to populate an array for displaying two-d graphics
   *  
   *   @param  a: the array to be filled
   *   @return none
   */  
  virtual void fill_2D_data(mglData& a) const {throw std::runtime_error("Cannot be here: Density::fill_2D_data");}


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
  double getInteractionEnergy(DFT_Vec &vdr) const { return Density_.cReal().dotWith(vdr);}

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

 protected:

  DFT_FFT Density_;  // The arrays for the real and fourier components of the density

  DFT_Vec vWall_;            ///< Array holding wall potential at each lattice point: i = pos(ix,iy,iz)
};

#endif // __LUTSKO__DENSITY_ARRAY__
