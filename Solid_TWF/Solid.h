#ifndef __LUTSKO_PERIOIDC__
#define __LUTSKO_PERIOIDC__

#include "Density.h"
#include "Display.h"
 
/**
  *  @brief Class that specializes DensityArray to a purely perdiodic box
  */  

class Solid : public Density
{
 public:
  /**
  *   @brief  Default  constructor for Solid 
  *  
  *   @param  dx is lattice spacing: assumed to be the same in all directions
  *   @param  L[] are the dimensions of the physical region
  *   @param  PointsPerHardSphere is number of points per hard-sphere diameter 
  *   @return nothing 
  */  
 Solid(double dx, double L[], double a, int ncopy) 
   : Density(dx,L), a_(a), ncells_(ncopy), display(NULL) { vWall_.zeros(Ntot_); }

  /**
  *   @brief  This is always zero
  *  
  *   @param  ijk are the positions
  *   @return  
  */  
  virtual double getField(int i, int j, int k) const { return 0.0;}
  
  /**
  *   @brief  Generates an intial guess at the density
  *  
  *   @param  f(z) gives the density as a function of position z
  *   @return  
  */  
  virtual void initialize(double alpha, double prefac, double Ntarget)
  {
    display = new Display(Nx_,Ny_);

    double atoms[4][3];

    atoms[0][0] = 0.0;
    atoms[0][1] = 0.0;
    atoms[0][2] = 0.0;

    atoms[1][0] = a_/2;
    atoms[1][1] = a_/2;
    atoms[1][2] = 0.0;

    atoms[2][0] = a_/2;
    atoms[2][1] = 0.0;
    atoms[2][2] = a_/2;

    atoms[3][0] = 0.0;
    atoms[3][1] = a_/2;
    atoms[3][2] = a_/2;


    cout << endl;
    cout << "Density Initialization ..." << endl;
    cout << "alpha = " << alpha << " prefac = " << prefac << endl;

    long Nweak = 0; // number of cells with low density
    double densityLimit = 0.01;
    
    for(int i=0;i<Nx_;i++)
      for(int j=0;j<Ny_;j++)
	for(int k=0;k<Nz_; k++)
	  {
	    double x = getX(i);
	    double y = getY(j);
	    double z = getZ(k);

	    double dsum = 0.0;

	    for(int icell=0;icell < ncells_; icell++)
	      for(int jcell=0;jcell < ncells_; jcell++)
		for(int kcell=0;kcell < ncells_; kcell++)
		  for(int l=0;l<4;l++)	     
		    {
		      double dx = fabs(x-atoms[l][0]-icell*a_); if(dx > L_[0]/2) dx -= L_[0];
		      double dy = fabs(y-atoms[l][1]-jcell*a_); if(dy > L_[1]/2) dy -= L_[1];
		      double dz = fabs(z-atoms[l][2]-kcell*a_); if(dz > L_[2]/2) dz -= L_[2];

		      double r2 = dx*dx+dy*dy+dz*dz;
		      dsum += prefac*pow(alpha/M_PI,1.5)*exp(-alpha*r2);
		    }
	    dsum = max(dsum, SMALL_VALUE);
	    set_Density_Elem(i,j,k,dsum);
	    if(dsum < densityLimit) Nweak++;
	  }


      // Adjust to have correct mass
  double Ncurrent = getNumberAtoms();

  double density_to_add = (Ntarget-Ncurrent)/(Nweak*dV());

  cout << "Ntarget = " << Ntarget << " Ncurrent = " << Ncurrent << " density adjustment = " << density_to_add << endl;


  for(int i=0;i<Nx_;i++)
    for(int j=0;j<Ny_;j++)
      for(int k=0;k<Nz_; k++)
	if(getDensity(i,j,k) < densityLimit) set_Density_Elem(i,j,k,getDensity(i,j,k)+1e-2); //density_to_add);


    
    cout << endl;
  }


  virtual void doDisplay(string &title, string &file) const
  {
    for(int i=0;i<Nx_;i++)
      for(int j=0;j<Ny_;j++)
	  display->setData(i+Nx_*j, getDensity(i,j,Nz_/2));

    display->doDisplay(title,file);
  }

  
  virtual void initialize_2D_data( mglData& a) const {a.Create(Nx_,Ny_);}
  virtual void fill_2D_data(mglData& a) const 
  {  
    for(int i=0;i<Nx_;i++)
      for(int j=0;j<Ny_;j++)
	a.a[i+Nx_*j] = min(getDensity(i,j,Nz_/2),1.0);
  }

 private:
  double a_;
  int ncells_;
  Display *display;
  
};



#endif // __LUTSKO_PERIOIDC__
