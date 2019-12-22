#ifndef __LUTSKO_SOLIDFCC__
#define __LUTSKO_SOLIDFCC__

#include <iostream>
#include <vector>
#include <cmath>

#include "Density.h"

#ifdef USE_MGL
#include "Display.h" 
#endif

#ifdef USE_GRACE
#include "Grace.h"
#endif

/**
  *  @brief Class that specializes Density to a FCC cell
  */  

class SolidFCC : public Density
{
 public:
  /**
   *   @brief  Default  constructor for a FCC cell in a periodic box 
   *  
   *   @param  dx is lattice spacing: assumed to be the same in all directions
   *   @param  L[] are the dimensions of the physical region (FCC cube side)
   *   @return nothing 
   */  
 SolidFCC(double dx, double L[])
   : Density(dx,L)
   {
#ifdef USE_GRACE
      grace_ = new Grace();
#endif
#ifdef USE_MGL
    // This is an object that writes a png snapshot whenever doDisplay gets called.
    display_ = new Display(Nx_,Ny_);
#endif
      
   }

   ~SolidFCC()
   {
#ifdef USE_GRACE
     if(grace_  != NULL) {grace_->close(); delete grace_;}
#endif
#ifdef USE_MGL     
     if(display_ != NULL) delete display_;
#endif
   }
   
  /**
   *   @brief  Generates an initial guess at the density
   *
   *   @detailed This guess is the sum of multiple gaussian profiles located at each FCC atom position for one cell
   *  
   *   @param  cell dimensions
   *   @param  standard deviation of gaussian profile
   *   @param  number of atoms in the FCC unit cell (taking vacancies into account)
   *   @return  
   */  
  virtual void initialize(double alpha, int ncells, double prefac, double dmin)
  {
    double a_latt = L_[0]/ncells;
    
    double atoms[4][3];

    atoms[0][0] = 0.0;
    atoms[0][1] = 0.0;
    atoms[0][2] = 0.0;

    atoms[1][0] = a_latt/2;
    atoms[1][1] = a_latt/2;
    atoms[1][2] = 0.0;

    atoms[2][0] = a_latt/2;
    atoms[2][1] = 0.0;
    atoms[2][2] = a_latt/2;

    atoms[3][0] = 0.0;
    atoms[3][1] = a_latt/2;
    atoms[3][2] = a_latt/2;

    for(int i=0;i<Nx_;i++)
      for(int j=0;j<Ny_;j++)
	for(int k=0;k<Nz_; k++)
	  {
	    double x = getX(i);
	    double y = getY(j);
	    double z = getZ(k);

	    double dsum = 0;

	    for(int icell=0;icell < ncells; icell++)
	      for(int jcell=0;jcell < ncells; jcell++)
		for(int kcell=0;kcell < ncells; kcell++)
		  for(int l=0;l<4;l++)	     
		    {
		      double dx = fabs(x-atoms[l][0]-icell*a_latt); if(dx > L_[0]/2) dx -= L_[0];
		      double dy = fabs(y-atoms[l][1]-jcell*a_latt); if(dy > L_[1]/2) dy -= L_[1];
		      double dz = fabs(z-atoms[l][2]-kcell*a_latt); if(dz > L_[2]/2) dz -= L_[2];

		      double r2 = dx*dx+dy*dy+dz*dz;
		      dsum += prefac*pow(alpha/M_PI,1.5)*exp(-alpha*r2);
		    }
	    if(dsum < dmin) dsum = dmin;
	    set_Density_Elem(i,j,k,dsum);	    
	  }
    //    string title("dummy");
    //    string file("dummy.dat");
    //    doDisplay(title,file,0);
  }

  virtual void initializeUniform(double density)
  {
    for(int i=0;i<Nx_;i++)
      for(int j=0;j<Ny_;j++)
	for(int k=0;k<Nz_; k++)
	  set_Density_Elem(i,j,k,density);
  }  

  // This gets called after every update and for each species: seq is the species number. 
  
  virtual void doDisplay(string &title, string &file, int seq) const
  {
#ifdef USE_GRACE    
    // Write to a grace window
    if(grace_ == NULL) return;

    grace_->deleteDataSet(seq);
    grace_->deleteDataSet(seq+1);

    int jx = 0;
    int jy = 0;
    int jz = 0;
    double dmax = 0;

    for(int ix = 0; ix < Nx(); ix++)
      for(int iy = 0; iy < Ny(); iy++)
	for(int iz = 0; iz < Nz(); iz++)
	    if(getDensity(ix,iy,iz) > dmax)
	      {
		dmax = getDensity(ix,iy,iz);
		jx = ix;
		jy = iy;
		jz = iz;
	      }
    bool pause = false;
    double px = 0;
    for(int i= 0; i <= Nz(); i++)
      {
	int iz = i-jz;
	while(iz > Nz()/2) iz -= Nz();
	while(iz < -Nz()/2) iz += Nz();

	
        double x = getZ(iz);
	grace_->addPoint(x,getDensity(jx,jy,i),seq);
	grace_->addPoint(x,std::max(1e-20, fabs(species_->getDF().get(get_PBC_Pos(jx,jy,i)))/dV()),seq+1);

	//	if(fabs(species_->getDF().get(get_PBC_Pos(jx,jy,i)))/dV() < 1e-20) {pause = true; px = x; cout << jx << " " << jy << " " << i << " " << getDensity(jx,jy,i) << " " << species_->getDF().get(get_PBC_Pos(jx,jy,i)) << endl;}
      }

    int iz = jz-Nz()/2;
    while(iz >= Nz()) iz -= Nz();
    while(iz < 0) iz += Nz();
    //    cout << "density = " <<  getDensity(jx,jy,iz) << " dF = " << species_->getDF().get(pos(jx,jy,iz)) << "df_max = " << species_->getDF().get(pos(jx,jy,iz)) << " jx = " << jx <<  " jy = " << jy << " iz = " << iz << " pos = " << pos(jx,jy,iz) << endl;    
    grace_->setTitle(title.c_str());
    grace_->redraw();
    if(pause){cout << px << endl;  grace_->pause();}
#endif
#ifdef USE_MGL
	// find max density on plane
	
	double max_density = 0;
	
	for(int i=0;i<Nx_;i++)
	  for(int j=0;j<Ny_;j++)
	    {
	      double d = getDensity(i,j,0);
	      if (d>max_density) max_density = d;
	    }
	
    // make a png snapshot --> Adjusted to max density
    
    for(int i=0;i<Nx_;i++)
      for(int j=0;j<Ny_;j++)
        display_->setData(i+Nx_*j, min(getDensity(i,j,jz)/max_density,1.0)); // TODO original line doesn't work when 1.0 is replaced by a bigger number in min(,) 

    stringstream title1; title1 << "Species " << seq << " " << title;
    stringstream file1;  file1 << seq << "_" << file;

    string t1 = title1.str();
    string f1 = file1.str();
    
    display_->doDisplay(t1, f1);
#endif        
  }  

  void setSpecies(Species* s) { species_ = s;}
  
 protected:
  int sequence_;

  Species *species_;
  
#ifdef USE_MGL
  Display *display_ = NULL;
#endif  
  

#ifdef USE_GRACE  
  Grace *grace_ = NULL;
#endif  

};



#endif // SLIT_PORE2__
