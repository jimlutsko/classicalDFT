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
  { init(L);}

  Lattice(double dx, double dy, double dz, double L[])
    : Nx_(0), Ny_(0), Nz_(0), dx_(dx), dy_(dy), dz_(dz)
  { init(L);}  
  
  Lattice()
  {
    Nx_ = Ny_ = Nz_ = 0;
    Ntot_ = Nout_ = 0;
    dx_ = dy_ = dz_ = 0;
    L_[0] =  L_[1] = L_[2] = 0;
  }

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


void init(double L[])
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
  
  // accessors
  double getDX() const {return dx_;}
  double getDY() const {return dy_;}
  double getDZ() const {return dz_;}
  long Nx() const { return Nx_;}
  long Ny() const { return Ny_;}
  long Nz() const { return Nz_;}
  double Lx() const { return L_[0];}
  double Ly() const { return L_[1];}
  double Lz() const { return L_[2];}
  long Ntot() const { return Ntot_;}

  long size() const { return Ntot_;}
  long Nout() const { return Nout_;}
  double dV() const { return dx_*dy_*dz_;}
  double getVolume() const { return L_[0]*L_[1]*L_[2];}
  
  // Translate index to cartesian coordinates
  double getX(int i) const { return dx_*(i-(Nx_-1)/2);}
  double getY(int i) const { return dy_*(i-(Ny_-1)/2);}
  double getZ(int i) const { return dz_*(i-(Nz_-1)/2);}


  // Cartesian indices into vector index
  inline long pos(int v[]) const { return pos(v[0],v[1],v[2]);}
  inline long pos(int i, int j, int k) const
  {
    if(i < 0 || i >= Nx_) throw std::runtime_error("ix out of bounds");
    if(j < 0 || j >= Ny_) throw std::runtime_error("iy out of bounds");
    if(k < 0 || k >= Nz_) throw std::runtime_error("iz out of bounds");
    return k + Nz_*(j +  Ny_*i);
  }
  // this is for complex arrays generated by FFTW3
  inline long kpos(int i, int j, int k) const
  {
    if(i < 0 || i >= Nx_) throw std::runtime_error("ix out of bounds");
    if(j < 0 || j >= Ny_) throw std::runtime_error("iy out of bounds");
    if(k < 0 || k >= 1+(Nz_/2)) throw std::runtime_error("iz out of bounds");
    return k + (1 + (Nz_/2))*(j +  Ny_*i);
  }  
  
  // Vector index to Cartesian indices
  inline void cartesian(long pos, int v[]) const
  {
    int ix,iy,iz; cartesian(pos,ix,iy,iz);
    v[0] = ix; v[1] = iy; v[2] = iz;
  }
  
  inline void cartesian(long pos, int &i, int &j, int &k) const { k = pos%(Nz_); pos = (pos-k)/Nz_; j = pos%Ny_; i = pos/Ny_;}  

  // Get vector index taking account of PBC
  virtual long get_PBC_Pos(int v[]) const {return get_PBC_Pos(v[0],v[1],v[2]);}  
  virtual long get_PBC_Pos(int ix, int iy, int iz) const
  { 
    putIntoBox(ix,iy,iz);
    return pos(ix,iy,iz);
  }

  // Apply PBC to cartesian indices
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

  ////////// Characterizing boundary points:
  // We map points on the boundary (ix or iy or iz = 0)
  // to an array of NxNy+Nx(Nz-1)+(Ny-1)(Nz-1) points


  long get_Nboundary() const { return Nx_*Ny_+Nx_*Nz_+Ny_*Nz_-Nx_-Ny_-Nz_+1;}
  // cartesian coordinates to position
  long boundary_pos(int v[]) const { return boundary_pos(v[0],v[1],v[2]);}

  long boundary_pos(int ix, int iy, int iz) const    
  {
    long pos = -1;

    while(ix < 0) ix += Nx_;
    while(iy < 0) iy += Ny_;
    while(iz < 0) iz += Nz_;

    while(ix > Nx_-1) ix -= Nx_;
    while(iy > Ny_-1) iy -= Ny_;
    while(iz > Nz_-1) iz -= Nz_;

    if(iz == 0)
      pos = ix*Ny_+iy;
    else {
      pos = Nx_*Ny_;
      if(iy == 0)
	pos += ix*(Nz_-1)+iz-1;
      else if(ix == 0) {
	pos += Nx_*(Nz_-1) + (iy-1)*(Nz_-1)+(iz-1);
      } else {
	throw std::runtime_error("Non-boundary point given to Lattice::boundary_pos");
      }
    }
    return pos;
  }


  void boundary_cartesian(long pos,int v[]) const
  {
    int ix,iy,iz;
    boundary_cartesian(pos,ix,iy,iz);
    v[0] = ix; v[1] = iy; v[2] = iz;
  }
  
  void boundary_cartesian(long pos,int &ix, int &iy, int &iz) const
  {
    ix = iy = iz = 0;
    
    if(pos < Nx_*Ny_)
      {
	ix = int(pos/Ny_);
	iy = pos-Ny_*ix;
      } else {
      pos -= Nx_*Ny_;
      if(pos < Nx_*(Nz_-1))
	{
	  ix = int(pos/(Nz_-1));
	  iz = 1+int(pos-ix*(Nz_-1));
	} else {
	pos -= Nx_*(Nz_-1);
	if(pos < (Ny_-1)*(Nz_-1))
	  {
	    iy = 1+int(pos/(Nz_-1));
	    iz = 1+int(pos - (iy-1)*(Nz_-1));
	  } else {
	  throw std::runtime_error("Non-boundary point given to Lattice::boundary_cartesian");
	}
      }      
    }
  }

  // a function for testing the boundary stuff
  void test_boundary_coding();

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
