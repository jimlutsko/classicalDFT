#ifndef __LUTSKO_PERIOIDC__
#define __LUTSKO_PERIOIDC__

#include "Density.h"
 
/**
  *  @brief Class that specializes Density to a purely perdiodic box
  */  

class Periodic : public Density
{
 public:
  /**
  *   @brief  Default  constructor for Periodic 
  *  
  *   @param  dx is lattice spacing: assumed to be the same in all directions
  *   @param  L[] are the dimensions of the physical region
  *   @param  PointsPerHardSphere is number of points per hard-sphere diameter 
  *   @return nothing 
  */  
 Periodic(double dx, double L[], Grace *g) 
   : Density(dx,L), g_(g), dft_(NULL) { vWall_.zeros(Ntot_); }


  void setDFT(DFT_VDW_Surfactant<RSLT> *d){dft_ = d;}
  
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
  virtual void initialize(double x0, double x1)
  {

    cout << "x0 = " << x0 << " or " << exp(-getField(10,0,0))*x0 << endl;

    for(int i=0;i<Nx_;i++)
      for(int j=0;j<Ny_;j++)
	for(int k=0;k<Nz_; k++)
	  {
	    double z = getZ(k);
	    if(z < 0)
	      set_Density_Elem(i,j,k,exp(-getField(i,j,k))*x0);
	    else set_Density_Elem(i,j,k,exp(-getField(i,j,k))*x1);
	  }
  }


  virtual void initialize_2D_data( mglData& a) const {a.Create(Nx_,Ny_);}
  virtual void fill_2D_data(mglData& a) const 
  {  
    for(int i=0;i<Nx_;i++)
      for(int j=0;j<Ny_;j++)
	a.a[i+Nx_*j] = min(getDensity(i,j,Nz_/2),1.0);
  }

  virtual void doDisplay(string &title, string &file) const
  {
    if(g_ == NULL) return;

    g_->deleteDataSet(0);
    if(dft_) g_->deleteDataSet(1);


    for(int i= 0; i < Nz(); i++)
      {
	double x = getZ(i);
	g_->addPoint(x,getDensity(Nx()/2,Ny()/2,i),0);
	if(dft_)
	  {
	    double y = dft_->getSurfactant(pos(Nx()/2,Ny()/2,i));
	    g_->addPoint(x,y,1);

	  }
      }
    
    g_->setTitle(title.c_str());
    g_->setColor(2,1);
    g_->redraw();
  }

 private:
  Grace *g_;
  DFT_VDW_Surfactant<RSLT> *dft_;
};



#endif // __LUTSKO_PERIOIDC__
