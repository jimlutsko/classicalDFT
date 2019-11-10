#ifndef __LUTSKO__LATTICE__
#define __LUTSKO__LATTICE__


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
   *  @brief  Translated single index into cartesian indices
   *  
   *   @param pos - the input index
   *   @param i - output x index
   *   @param j - output y index
   *   @param k - output z index
   */  
  inline void cartesian(long pos, long &i, long &j, long &k) const { k = pos%(Nz_); pos = (pos-k)/Nz_; j = pos%Ny_; i = pos/Ny_;}

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

#endif
