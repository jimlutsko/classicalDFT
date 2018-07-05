 #ifndef __LUTSKO_DROPLET__
#define __LUTSKO_DROPLET__

#include "Density.h"
 
/**
  *  @brief Class that specializes Density to a slit-pore geometry with a fixed spherical particle
  */  

class Droplet : public Density
{
 public:
  /**
  *   @brief  Default  constructor for Droplet 
  *  
  *   @param  dx is lattice spacing: assumed to be the same in all directions
  *   @param  L[] are the dimensions of the physical region
  *   @param  PointsPerHardSphere is number of points per hard-sphere diameter 
  *   @param  R is the radius of the fixed particle
  *   @param  zpos is the position of the fixed particle: it is always centered in x,y directions
  *   @return nothing 
  */  
 Droplet(double dx, double L[], int PointsPerHardSphere, double R, double zpos) 
   : Density(dx,L), R_(R), zPos_(zpos)
  {
  }

  /**
  *   @brief  Generates an intial guess at the density
  *  
  *   @param  density is the bulk density
  *   @param  d2 is unused
  *   @return  
  */  
  virtual void initialize(double density, double d2, Density *inputdensity = NULL)
  {
    Density::initialize(density,d2);

    for(int i=0;i<Nx_;i++)
      for(int j=0;j<Ny_;j++)
	for(int k=0;k<Nz_; k++)
	  {
	    if(inputdensity)
	      if(i == 0 || j == 0 || k == 0)
		set_Density_Elem(i,j,k,inputdensity->getDensity(i,j,k));
	    
	    double x = getX(i);
	    double y = getY(j);
	    double z = getZ(k);
	    
	    double r2 = x*x+y*y+z*z;

	    double dd = (r2 < R_*R_ ? density : d2);

	    set_Density_Elem(i,j,k,dd);
	  }
  }

  void initializeTo(Droplet& first, Droplet& last, double ratio)
  {
    Density::initialize(0.0,0.0);
    for(long i=0;i<Ntot();i++)
      set_Density_Elem(i, ratio*first.getDensity(i) + (1-ratio)*last.getDensity(i));    
  }

  void initializeToSqueezedImage(Droplet& last, double ratio, double NN)
  {
    cout << "ratio = " << ratio << endl;
    initialize(100,100);

    float squeeze = 1/ratio;

    int ix0 = int((Nx_-1)/2);
    int iy0 = int((Ny_-1)/2);
    int iz0 = int((Nz_-1)/2);

    double Mass = 0.0;
    double V = 0.0;
    double dmax = 0;
    for(int ix=0;ix<Nx_;ix++)
      {
	int IX = ix0+squeeze*(ix-ix0);
	if(IX < 0) IX = 0;
	if(IX > Nx_-1) IX = Nx_-1;
          for(int iy=0;iy<Ny_;iy++)
	    {
	      int IY = iy0+squeeze*(iy-iy0);
	      if(IY < 0) IY = 0;
	      if(IY > Ny_-1) IY = Ny_-1;
	        for(int iz=0;iz<Nz_;iz++)
		  {
		    int IZ = iz0+squeeze*(iz-iz0);
		    if(IZ < 0) IZ = 0;
		    if(IZ > Nz_-1) IZ = Nz_-1;
		    
		    double x = getX(IX);
		    double y = getY(IY);
		    double z = getZ(IZ);

		    set_Density_Elem(ix,iy,iz,last.getDensity(IX,IY,IZ));
		    double d = getDensity(ix,iy,iz);
		    Mass += d*dV();
		    V += dV();
		    if(d > dmax) dmax = d;
		  }
	    }
      }
    double deltaN = NN-Mass;
    double deltaV = Ntot()*dV()-V;
    double Nmax = dmax*Ntot()*dV();
    double alpha = (NN-Mass)/(Nmax-Mass);
    
    for(long i=0;i<Ntot();i++)
      {
	double d = getDensity(i);
	set_Density_Elem(i, d+(dmax-d)*alpha);
      }
      //      if(getDensity(i) > 90.0)
      //	set_Density_Elem(i, deltaN/deltaV);    
  }
  
    
  void initializeToCNT(double NN, double centralDensity, double ratio, double Rmax)
  {
    cout << "ratio = " << ratio << " NN = " << NN << endl;

    double Radius = ratio*Rmax;
    double Radius2 = Radius*Radius;
    
    Density::initialize(0.0,0.0);

    long Ncentral = 0;
    for(int ix=0;ix<Nx_;ix++)
          for(int iy=0;iy<Ny_;iy++)
	        for(int iz=0;iz<Nz_;iz++)
		  {
		    double x = getX(ix);
		    double y = getY(iy);
		    double z = getZ(iz);
		    double r2 = x*x+y*y+z*z;
		    if(r2 < Radius2)
		      {
			set_Density_Elem(ix,iy,iz,centralDensity);
			Ncentral++;
		      }
		  }

    double backgroundDensity = (NN - Ncentral*centralDensity*dV())/((Ntot_-Ncentral)*dV());

    for(int ix=0;ix<Nx_;ix++)
          for(int iy=0;iy<Ny_;iy++)
	        for(int iz=0;iz<Nz_;iz++)
		  {
		    double x = getX(ix);
		    double y = getX(iy);
		    double z = getX(iz);
		    double r2 = x*x+y*y+z*z;		    
		    if(r2 >= Radius2)
		      set_Density_Elem(ix,iy,iz,backgroundDensity);
		  }
  }


  // Used to generate 2D graphics: this holds the data to be displayed
    virtual void initialize_2D_data( mglData& a) const {a.Create(Nx_,Ny_);}

  // Used to generate 2D graphics: this is called to load the data to be displayed. 
  virtual void fill_2D_data(mglData& a) const 
  {  
    for(int i=0;i<Nx_;i++)
      for(int j=0;j<Ny_;j++)
	a.a[i+Nx_*j] = log(getDensity(i,j, int((Nz_-1)/2))); 	
	//	a.a[i+Nx_*j] = getDensity(i,j, int(0.85*Nz_));
  }

  void translate(ofstream &of)
  {
    for(int i = 0; i < Nx_; i++)
      for(int j = 0; j < Ny_; j++)
	for(int k = 0; k < Nz_; k++)
	  of << i*dx_ << "," << j*dy_ << "," << k*dz_ << "," << getDensity(i,j,k) << endl;
  }

 protected:
  double R_;        ///< Radius of the spherical cap droplet
  double zPos_;     ///< Position of the center of the sphere defining the droplet
};



#endif // __LUTSKO_DROPLET__
