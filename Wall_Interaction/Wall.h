#ifndef __LUTSKO_SLIT_PORE2__
#define __LUTSKO_SLIT_PORE2__

#include "Density.h"

#ifdef USE_MGL
#include "Display.h" 
#endif

#ifdef USE_GRACE
#include "Grace.h"
#endif

/**
  *  @brief Class that specializes Density to a slit-pore geometry with a fixed spherical particle
  */  

class Wall : public Density
{
 public:
  /**
  *   @brief  Default  constructor for Wall 
  *  
  *   @param  dx is lattice spacing: assumed to be the same in all directions
  *   @param  L[] are the dimensions of the physical region
  *   @param  sigV is the length scale of the wall potential
  *   @param  epsV is the energy scale of the wall potential
  *   @return nothing 
  */  
Wall(double dx, double L[], double sigV, double epsV, double hsd, void *grace) 
: Density(dx,L),  sigV_(sigV), epsV_(epsV), hsd_(hsd)
  {
#ifdef USE_GRACE
grace_ = (Grace*) grace;
#endif

    for(int ix=0;ix<Nx_;ix++)
      for(int iy=0;iy<Ny_;iy++)
	for(int iz=0;iz<Nz_;iz++)
	  {
	    double V = 0;
	    double z = getZ(iz);
	    double zwall = (L_[2]/2)-1;
    
	    if(z > -zwall)  
	      {
		double y = (z+zwall)/sigV_;       
		V += 0.4*pow(y,-10)-pow(y,-4)-(sqrt(2.0)/3)*pow(y+0.61/sqrt(2),-3);		
	      }
	    if(z <  zwall) 
	      {
		double y = (zwall-z)/sigV_;	
		V += 0.4*pow(y,-10)-pow(y,-4)-(sqrt(2.0)/3)*pow(y+0.61/sqrt(2),-3);
	      }	    
	    if(z >= zwall-hsd_/2 || z <= -zwall+hsd_/2) V = 100/epsV_;

	    double x = getX(ix); 
	    double y = getY(iy); 
	    
	    vWall_.set(pos(ix,iy,iz),2*M_PI*epsV_*V);
	  }
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

    for(int i=0;i<Nx_;i++)
      for(int j=0;j<Ny_;j++)
	for(int k=0;k<Nz_; k++)
	  {
	      double z = getZ(k);
	      double d = exp(-getField(i,j,k)*density);
	      set_Density_Elem(i,j,k,(d > 0.1 ? (i < Nx_/2 ? density : d2) : 0));
	  }
  }

  virtual void doDisplay(string &title, string &file, int seq) const
  {
#ifdef USE_GRACE
    if(grace_ == NULL) return;

    grace_->deleteDataSet(seq);
    //    grace_->deleteDataSet(seq+2);

    double emax = 0;
    
    for(int i= 0; i < Nz(); i++)
      {
	double x = getZ(i);
	//	grace_->addPoint(x,(M_PI*hsd_*hsd_*hsd_/6)*getDensity(Nx()/2,Ny()/2,i),seq+2);
	//	grace_->addPoint(x,(M_PI*hsd_*hsd_*hsd_/6)*getDensity(3*Nx()/4,3*Ny()/4,i),seq);
	grace_->addPoint(x,getDensity(3*Nx()/4,3*Ny()/4,i),seq);
	emax = max(emax,(M_PI*hsd_*hsd_*hsd_/6)*getDensity(Nx()/2,Ny()/2,i));
	/*	if(dft_)
	  {
	    double y = dft_->getSurfactant(pos(Nx()/2,Ny()/2,i));
	    grace_->addPoint(x,y,1);

	  }
	*/
      }
      grace_->setTitle(title.c_str());
    //    grace_->setColor(2,1);
    grace_->redraw();

#endif    
    
  }  

 protected:
  double sigV_;     ///< length scale of wall potential
  double epsV_;     ///< energy scale of wall potential
  int sequence_;
  double hsd_;

#ifdef USE_GRACE
      Grace *grace_;
#endif
};



#endif // __LUTSKO_SLIT_PORE__
