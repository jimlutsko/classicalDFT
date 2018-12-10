#ifndef __LUTSKO_DROPLET__
#define __LUTSKO_DROPLET__

#include "Density.h"
#include "Display.h"

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
  *   @param  sigV is the length scale of the wall potential
  *   @param  epsV is the energy scale of the wall potential
  *   @return nothing 
  */  
 Droplet(double dx, double L[], int PointsPerHardSphere, double R, double zpos, double sigV, double epsV) 
   : Density(dx,L), R_(R), zPos_(zpos), sigV_(sigV), epsV_(epsV)
  {
    display_ = new Display(Nx_, Ny_+5);
  }

  /**
  *   @brief  Generates an intial guess at the density
  *  
  *   @param  density is the bulk density
  *   @param  d2 is unused
  *   @return  
  */  
  virtual void initialize(double density, double d2)
  {
    Density::initialize(density,d2);

    double zwall = (L_[2]/2)-1;

    cout << "zwall = " << zwall << endl;

    for(int i=0;i<Nx_;i++)
      for(int j=0;j<Ny_;j++)
	for(int k=0;k<Nz_; k++)
	  {
	    double x = getX(i);
	    double y = getY(j);
	    double z = getZ(k); //-zwall;
	    
	    double r2 = x*x+y*y+z*z;

	    double dd = (r2 < R_*R_ ? density : d2);

	    double d = exp(-getField(i,j,k)*dd);
	   
	    set_Density_Elem(i,j,k,max((d > 0.1 ? dd : 0),SMALL_VALUE));
	    if(fabs(dd-max((d > 0.1 ? dd : 0),SMALL_VALUE))) cout << "What???" << endl;
	    //	    set_Density_Elem(i,j,k,dd,SMALL_VALUE));
	  }
  }


  void doDisplay(string &title, string &file)
  {


    Density::doDisplay(title, file);
    
    for(int i=0;i<Nx_;i++)
      for(int j=0;j<Ny_;j++)
	display_->setData(i+Nx_*j, getDensity(i,j, int((Nz_-1)/2)));

    display_->doDisplay(title, file);
  }
  /*
  
  // Used to generate 2D graphics: this holds the data to be displayed
    virtual void initialize_2D_data( mglData& a) const {a.Create(Nx_,Ny_+5);}

  // Used to generate 2D graphics: this is called to load the data to be displayed. 
  virtual void fill_2D_data(mglData& a) const 
  {  
    for(int i=0;i<Nx_;i++)
      for(int j=0;j<Ny_;j++)
	a.a[i+Nx_*j] = getDensity(i,j, int((Nz_-1)/2)); //min(getDensity(i,Ny_/2,j-1),1.0);
    //	a.a[i+Nx_*j] = getDensity(i,Ny_/2,j-1); //min(getDensity(i,Ny_/2,j-1),1.0);

    /
    for(int i=0;i<Nx_;i++)
      for(int j=Nz_;j<Nz_+5;j++)
	a.a[i+Nx_*j] = 10.0;
    /
  }
*/
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
  double sigV_;     ///< length scale of wall potential
  double epsV_;     ///< energy scale of wall potential

  Display *display_;
};



#endif // __LUTSKO_DROPLET__