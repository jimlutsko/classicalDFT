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
 Droplet(double dx, double L[], int PointsPerHardSphere, string &infile, double alpha, double baseline)
   : Density(dx,L), infile_(infile), alpha_(alpha), baseline_(baseline)
  {}

  /**
  *   @brief  Generates an intial guess at the density
  *  
  *   @param  density is the bulk density
  *   @param  d2 is unused
  *   @return  
  */  
  virtual void initialize()
  {
    Density::initialize(1,1);

    ifstream in(infile_.c_str());
    if(!in.good())
      throw std::runtime_error("Bad input file");

    string line;
    getline(in,line); // get empty line

    vector< vector<double> > pos;
    
    double xav,yav,zav;
    xav = yav = zav = 0.0;
    
    double x,y,z;
    while(in >> x >> y >> z)
      {
	vector<double> d;
	d.push_back(x);
	d.push_back(y);
	d.push_back(z);
	pos.push_back(d);
	xav += x;
	yav += y;
	zav += z;
      }
    
    int N = pos.size();
    xav /= N;
    yav /= N;
    zav /= N;
    
    cout << "av pos: (" << xav << "," << yav << "," << zav << ")" << endl;

    // make into integers
    xav = int(xav);
    yav = int(yav);
    zav = int(zav);

    cout << "av pos: (" << xav << "," << yav << "," << zav << ")" << endl;
    
    for(int i=0;i<pos.size();i++)    
      {pos[i][0] -= xav; pos[i][1] -= yav; pos[i][2] -= zav;}

    double norm = pow(alpha_/M_PI,1.5);
    
    for(int i=0;i<Nx_;i++)
      for(int j=0;j<Ny_;j++)
	for(int k=0;k<Nz_; k++)
	  {
	    double x = getX(i);
	    double y = getY(j);
	    double z = getZ(k);

	    double density = 0;
	    
	    for(int l=0; l<pos.size(); l++)
	      {
	        double r2 = (x-pos[l][0])*(x-pos[l][0])+(y-pos[l][1])*(y-pos[l][1])+(z-pos[l][2])*(z-pos[l][2]);
		density += exp(-alpha_*r2)*norm;
	      }
	    if(density < 1e-10) density += baseline_;
	    set_Density_Elem(i,j,k,max(density,SMALL_VALUE));		
	  }
    
    cout << setprecision(12) << "Natoms = " << getNumberAtoms() << endl;


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
  string infile_;
  double alpha_;
  double baseline_;
  
  Display *display_;
};



#endif // __LUTSKO_DROPLET__
