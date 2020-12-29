#ifndef __LUTSKO__LATTICE__
#define __LUTSKO__LATTICE__

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

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
	{
	  cout << "L[0] = " << L[0] << " L_[0] = " << L_[0] << " Nx*dx = " << Nx_*dx_ << endl;
	  throw std::runtime_error("Box length incommensurate with lattice in x direction");
	}
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
  *   @brief  Emptry constructor
  */  

  Lattice()
    {
      Nx_ = Ny_ = Nz_ = 0;

      Ntot_ = Nout_ = 0;

      dx_ = dy_ = dz_ = 0;

      L_[0] =  L_[1] = L_[2] = 0;
    }
  

  /**
  *   @brief  Translate a (Cartesian) index into a position in the box.
  */  
  double getX(int i) const { return dx_*(i-(Nx_-1)/2);}
  double getY(int i) const { return dy_*(i-(Ny_-1)/2);}
  double getZ(int i) const { return dz_*(i-(Nz_-1)/2);}


  // Translate a position into the nearest coordinate
  int getIX(double x) const { return round(x/dx_)+(Nx_-1)/2;}
  int getIY(double y) const { return round(y/dy_)+(Ny_-1)/2;}
  int getIZ(double z) const { return round(z/dz_)+(Nz_-1)/2;}
  
  // spacings
  double getDX() const {return dx_;}
  double getDY() const {return dy_;}
  double getDZ() const {return dz_;}

  /**
   *  @brief  Translated (ix,iy,iz) cartesian indices into single index
   *  
   *   @return i(ix,iy,iz)
   */  
  inline long pos(long i, long j, long k) const
  {
    if(i < 0 || i >= Nx_) throw std::runtime_error("ix out of bounds");
    if(j < 0 || j >= Ny_) throw std::runtime_error("iy out of bounds");
    if(k < 0 || k >= Nz_) throw std::runtime_error("iz out of bounds");
    return k + Nz_*(j +  Ny_*i);
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
   *  @brief  Translated single index into cartesian indices
   *  
   *   @param pos - the input index
   *   @param i - output x index
   *   @param j - output y index
   *   @param k - output z index
   */  
  inline void cartesian(long pos, long &i, long &j, long &k) const { k = pos%(Nz_); pos = (pos-k)/Nz_; j = pos%Ny_; i = pos/Ny_;}

  long Nx() const { return Nx_;}
  long Ny() const { return Ny_;}
  long Nz() const { return Nz_;}

  double Lx() const { return L_[0];}
  double Ly() const { return L_[1];}
  double Lz() const { return L_[2];}

  long Ntot() const { return Ntot_;}
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
    while(ix < 0)  ix += Nx_; 
    while(ix >= Nx_) ix -= Nx_;
    
    while(iy < 0)  iy += Ny_; 
    while(iy >= Ny_) iy -= Ny_;
    
    while(iz < 0)  iz += Nz_; 
    while(iz >= Nz_) iz -= Nz_; 
  }

  friend ostream &operator<<(ostream &of, const Lattice &l) 
  {    
    of.write((char*) &l.Nx_, sizeof(long));
    of.write((char*) &l.Ny_, sizeof(long));
    of.write((char*) &l.Nz_, sizeof(long));

    of.write((char*) &l.Ntot_, sizeof(long));
    of.write((char*) &l.Nout_, sizeof(long));    

    of.write((char*) &l.dx_, sizeof(double));
    of.write((char*) &l.dy_, sizeof(double));
    of.write((char*) &l.dz_, sizeof(double));

    of.write((char*) &l.L_, 3*sizeof(double));
    
    return of;
  }

  friend istream &operator>>(istream  &in, Lattice &l )     
  {
    in.read((char*) &l.Nx_, sizeof(long));
    in.read((char*) &l.Ny_, sizeof(long));
    in.read((char*) &l.Nz_, sizeof(long));
    
    in.read((char*) &l.Ntot_, sizeof(long));
    in.read((char*) &l.Nout_, sizeof(long));    
    
    in.read((char*) &l.dx_, sizeof(double));
    in.read((char*) &l.dy_, sizeof(double));
    in.read((char*) &l.dz_, sizeof(double));
    
    in.read((char*) &l.L_, 3*sizeof(double));
    
    return in;
  }    
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & Nx_;
    ar & Ny_;
    ar & Nz_;

    ar & Ntot_;
    ar & Nout_;

    ar & dx_;
    ar & dy_;
    ar & dz_;

    ar & L_;
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

#endif
