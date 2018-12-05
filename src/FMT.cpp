#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

using namespace std;

#ifdef USE_OMP
#include <omp.h>
#endif

#ifdef USE_MPI
#include <mpi.h>  
#endif

#include "FMT.h"
#include "Enskog.h"


ostringstream Eta_Too_Large_Exception::cnvt;

double sigma2 = 1e-4;


FMT::FMT(Lattice &lattice,  double hsd, string& pointsFile, double hsd1)
  : dx_(lattice.getDX()), dy_(lattice.getDY()), dz_(lattice.getDZ()), hsd0_(hsd), hsd1_(hsd1), etaMax_(1e30)
{
  long Nx = lattice.Nx();
  long Ny = lattice.Ny();
  long Nz = lattice.Nz();

  
  cout << endl;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << "/////  Calculating FMT weighting functions" << endl;
  cout << "/////  and setting up FFT's " << endl;

  // Initialize the array of weighted densities
  initializeWeightedDensities(d0_,hsd0_,Nx,Ny,Nz, pointsFile);

  if(hsd1_ > 0)
    initializeWeightedDensities(d1_,hsd1_,Nx,Ny,Nz,pointsFile);

  long Nout = lattice.Nout();
  long Ntot = lattice.Ntot();

  dPhi_.initialize(Nx,Ny,Nz);

  cout << "/////  Finished.  " << endl;
  cout << "///////////////////////////////////////////////////////////" << endl;

}

void FMT::initializeWeightedDensities(vector<FMT_Weighted_Density> &dd, double hsd, long Nx, long Ny, long Nz, string& pointsFile)
{
  dd.resize(11);
  for(FMT_Weighted_Density &d: dd)
    d.initialize(Nx, Ny, Nz);

  generateWeights(dd,hsd,pointsFile,Nx,Ny,Nz);

  for(FMT_Weighted_Density &d: dd)
    d.transformWeights();


  
  /*
  for(int i=0;i<Nx;i++)
    for(int j=0;j<Ny;j++)
      for(int k=0;k<(Nz/2)+1;k++)
	{
	  double kx = 2*M_PI*(2*i > Nx ? i-Nx : i)*hsd/(Nx*dx_);
	  double ky = 2*M_PI*(2*j > Ny ? j-Ny : j)*hsd/(Ny*dy_);
	  double kz = 2*M_PI*(2*k > Nz ? k-Nz : k)*hsd/(Nz*dz_);
	  double kk = sqrt(kx*kx+ky*ky+kz*kz);

	  if(i != 0 || j != 0 || k != 0)
	    {
	      kx /= kk;
	      ky /= kk;
	      kz /= kk;
	    }

	  long pos = k+((Nz/2)+1)*(j+Ny*i);
	  double fac = 1.0/(Nx*Ny*Nz);

	  double j0 = gsl_sf_bessel_j0(kk/2);
	  double j1 = gsl_sf_bessel_j1(kk/2);
	  double j2 = gsl_sf_bessel_j2(kk/2);

	  dd[0].setWk(pos,fac*(M_PI*hsd*hsd*hsd/6)*(j0+j2),0.0);   
	  dd[1].setWk(pos,fac*(M_PI*hsd*hsd)*j0,0.0);
	  //v
	  dd[2].setWk(pos,0.0,fac*(M_PI*hsd*hsd)*kx*j1);
	  dd[3].setWk(pos,0.0,fac*(M_PI*hsd*hsd)*ky*j1);
	  dd[4].setWk(pos,0.0,fac*(M_PI*hsd*hsd)*kz*j1);
	  //T
	  dd[5].setWk(pos,fac*(M_PI*hsd*hsd/3)*(j0+j2-3*kx*kx*j2),0.0);
	  dd[6].setWk(pos,fac*(M_PI*hsd*hsd/3)*(     -3*kx*ky*j2),0.0);
	  dd[7].setWk(pos,fac*(M_PI*hsd*hsd/3)*(     -3*kx*kz*j2),0.0);
	  dd[8].setWk(pos,fac*(M_PI*hsd*hsd/3)*(j0+j2-3*ky*ky*j2),0.0);
	  dd[9].setWk(pos,fac*(M_PI*hsd*hsd/3)*(     -3*ky*kz*j2),0.0);

	  dd[10].setWk(pos,fac*(M_PI*hsd*hsd/3)*(j0+j2-3*kz*kz*j2),0.0);
	}
  */
  cout << "Lx = " << Nx*dx_ << endl;
  cout << "Ly = " << Ny*dy_ << endl;
  cout << "Lz = " << Nz*dz_ << endl;
  
}


// The idea here is as follows. Consider the s weighted density which is just an integral over a sphere of radius hsd/2.
// In principle, we will need the volume integral
//    s(r_i) = int w(|r_i - r'|) rho(r') delta(r'-(hsd/2)) dr'
// The first step is to approximate as some discrete sum on the spherical shell:
//    s(r_i) ~ sum_{k} w(|r_i - r'_{k}*(hsd/2)|) rho(r'_{k}*(hsd/2)) ww_{k}
// where I am assuming that the points r'_{k} are given for a unit sphere.
// Next, we evaluate rho using trilinear interpolation. This will involve a sum over the 8 corners of the cubic cell
// containing the integration point (i.e. its nearest neighbors on the lattice)
// with corresponding weights (which I call u):
//    s(r_i) ~ sum_{k} w(|r_i - r'_{k}*(hsd/2)|) ww_{k} sum_{l \in nearest neighbors of r'_{k}*(hsd/2)} rho(r_{l}) u_{l}(r'_{k}*(hsd/2))
// or
//     s(r_i) ~ sum_{j} A_{ij} rho(r_{k})
// where
//       A_{ij} = sum_{k} w(|r_i - r'_{k}*(hsd/2)|) ww_{k} u_{l}(r'_{k}*(hsd/2))  
// with the understanding that the weight u_{l} is zero if l is not a nearest neighbor of the point r'_{k}*(hsd/2).
// Finally, A is translationally invariant so we really just need it for i=0. It is the arrays A_{0,j} (i.e. the "weights") which are computed here. 

void FMT::generateWeights(vector<FMT_Weighted_Density> &densities, double hsd, string& pointsFile, long Nx, long Ny, long Nz)
{
  // Generate integration points for spherical surface

  cout << endl;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << "/////  Generating integration points on sphere " << endl;

  vector < vector<double> > points; 

#ifndef USE_MPI  
  // Read points from file : C++ way
  
  ifstream in(pointsFile.c_str());
  if(!in.good())
    throw std::runtime_error("input file cannot be opened");
#else
  // Read points from file : C code for use with MPI functions
    MPI_Status    status;
    MPI_File      fh;

    int error = MPI_File_open(MPI_COMM_WORLD, pointsFile.c_str(),
                  MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if(error != MPI_SUCCESS) 
      throw std::runtime_error("Could not open points file");

  // get num bytes
  MPI_Offset numbytes;
  error = MPI_File_get_size(fh, &numbytes);
  if(error != MPI_SUCCESS) 
    throw std::runtime_error("Could not get points file size");

  // allocate buffer
  char *buffer = (char*)calloc(numbytes, sizeof(char));	

  if(buffer == NULL)
    throw std::runtime_error("could not allocate buffer for points file");

  // copy
  MPI_File_seek(fh, 0, MPI_SEEK_SET);
  error = MPI_File_read(fh, buffer, numbytes, MPI_BYTE, &status);
  if(error != MPI_SUCCESS) 
    throw std::runtime_error("Could not read points file");

  MPI_File_close(&fh);
  
  string sbuffer(buffer,numbytes/sizeof(char));

  delete buffer;
  istringstream in(sbuffer);
#endif
  
  for(string buf; getline(in,buf);)
    {
      vector<double> line;
      
      stringstream st(buf);
      double d;
      while(st >> d) {line.push_back(d);}
      line.push_back(4*M_PI);

      double r2 = line[0]*line[0]+line[1]*line[1]+line[2]*line[2];
      double r = sqrt(r2);
      line[0] /= r;
      line[1] /= r;
      line[2] /= r;

      points.push_back(line);
    }
  // adjust the weights
  for(vector<double> &point : points)
    point[3] /= points.size();

  cout << "/////  Finished: there are " << points.size() << " points " << endl;
  cout << "///////////////////////////////////////////////////////////" << endl;


  // Add up the weights for each point.
  // We throw in all permutations (iperm loop) of the integration points and all reflections (is loop)
  // to establish as much symmetry as possible. Probably unnecessary.    
  cout << endl;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << "/////  Generating integration points for sphere volume " << endl;

  // This is for the radial integral
  int Nr = 128*4;
  double r = 0.5*hsd; // hard sphere radius
  // Get some Gauss-Lagrange integration points and weights for the radial integral
  gsl_integration_glfixed_table *tr = gsl_integration_glfixed_table_alloc(Nr);

  long count = 0;
  for(int pos=0;pos < points.size(); pos++)
    for(int iperm = 0; iperm < 6; iperm++)      // add permutations
      {
	int ii = 0;
	int jj = 1;
	int kk = 2;
	if(iperm == 0) {ii = 0; jj = 1; kk = 2;}
	if(iperm == 1) {ii = 0; jj = 2; kk = 1;}
	if(iperm == 2) {ii = 1; jj = 0; kk = 2;}
	if(iperm == 3) {ii = 1; jj = 2; kk = 0;}
	if(iperm == 4) {ii = 2; jj = 0; kk = 1;}
	if(iperm == 5) {ii = 2; jj = 1; kk = 0;}

	double x0 = r*points[pos][ii];
	double y0 = r*points[pos][jj];
	double z0 = r*points[pos][kk];

	for(int is=0;is<8;is++) // add in reflections too
	  {
	    double x = x0;
	    double y = y0;
	    double z = z0;
	    // is == 0: do nothing
	    if(is == 1) {x = -x;}
	    if(is == 2) {y = -y;}
	    if(is == 3) {z = -z;}
	    if(is == 4) {x = -x; y = -y;}
	    if(is == 5) {x = -x; z = -z;}
	    if(is == 6) {y = -y; z = -z;}
	    if(is == 7) {x = -x; y = -y; z = -z;}


	    // scale the original weight by the number of redundancies (6 x 7 or 6 x 8 depending on if we keep the last reflection)
	    double ww = points[pos][3]/48.0;

	    double v[3]; // used for the vector and tensor densities
	    v[0] = x/r;
	    v[1] = y/r;
	    v[2] = z/r;

	    // Find "smallest" corner of the cube that contains this point. 
	    int ix0 = int(x/dx_);
	    int iy0 = int(y/dy_);
	    int iz0 = int(z/dz_);

	    if(x < 0) ix0 -=1;
	    if(y < 0) iy0 -=1;
	    if(z < 0) iz0 -=1;

	    double ddx = (x/dx_)-ix0;
	    double ddy = (y/dy_)-iy0;
	    double ddz = (z/dz_)-iz0;

	    // Do the trilinear interpolation
	    //  and add to the weights
	    for(int io=0;io<2;io++)
	      for(int jo=0;jo<2;jo++)
		for(int ko=0;ko<2;ko++)		  
		  {
		    double cc = (io == 0 ? 1-ddx : ddx)*(jo == 0 ? 1-ddy : ddy)*(ko == 0 ? 1-ddz : ddz); // trilinear interpolation

		    int nx = ix0+io;
		    int ny = iy0+jo;
		    int nz = iz0+ko;

		    if(nx < 0) nx += Nx;
		    if(ny < 0) ny += Ny;
		    if(nz < 0) nz += Nz;

		    long pos = nz+Nz*(ny+Ny*nx);
	      
		    S(densities).addToWeight(pos,cc*ww*r*r);
	      
		    for(int jv=0;jv<3;jv++)
		      {
			V(jv, densities).addToWeight(pos, -cc*ww*r*r*v[jv]);

			for(int kT=0;kT<=jv;kT++)
			  T(jv,kT, densities).addToWeight(pos,cc*ww*r*r*v[jv]*v[kT]);	
		      }
		    count++;
		  }

	    // Now add in the radial integral as well ...
	    // This just adds another sum 
	    for(int k=0;k<Nr;k++)
	      {
		double r;
		double wr;
	      
		gsl_integration_glfixed_point(0, 0.5*hsd, k, &r, &wr, tr);
	      
		double x = r*points[pos][ii];
		double y = r*points[pos][jj];
		double z = r*points[pos][kk];

		// is = 0: do nothing
		if(is == 1) {x = -x;}
		if(is == 2) {y = -y;}
		if(is == 3) {z = -z;}
		if(is == 4) {x = -x; y = -y;}
		if(is == 5) {x = -x; z = -z;}
		if(is == 6) {y = -y; z = -z;}
		if(is == 7) {x = -x; y = -y; z = -z;}

		double ww = points[pos][3]/48.0;
	      
		// Find the cube that contains this point. 
		int ix0 = int(x/dx_);
		int iy0 = int(y/dy_);
		int iz0 = int(z/dz_);

		if(x < 0) ix0 -=1;
		if(y < 0) iy0 -=1;
		if(z < 0) iz0 -=1;

		double ddx = (x/dx_)-ix0;
		double ddy = (y/dy_)-iy0;
		double ddz = (z/dz_)-iz0;

		for(int io=0;io<2;io++)
		  for(int jo=0;jo<2;jo++)
		    for(int ko=0;ko<2;ko++)		  
		      {

			int nx = ix0+io;
			int ny = iy0+jo;
			int nz = iz0+ko;

			if(nx < 0) nx += Nx;
			if(ny < 0) ny += Ny;
			if(nz < 0) nz += Nz;

			long pos = nz+Nz*(ny+Ny*nx);

			double cc = (io == 0 ? 1-ddx : ddx)*(jo == 0 ? 1-ddy : ddy)*(ko == 0 ? 1-ddz : ddz); // trilinear interpolation
			Eta(densities).addToWeight(pos,cc*ww*wr*r*r);      
		      }
	      }
	  }
      }
  // This is some overly-cautious adjusting to assure that some exact relations are satsified. (Numerical noise may have violated them slightly).
  for(long pos = 0;pos <Nx*Ny*Nz; pos++) 
    {
      double sum = T(0,0, densities).getWeight(pos)+T(1,1, densities).getWeight(pos)+T(2,2, densities).getWeight(pos);
      if(fabs(sum) > 1e-16)
	{
	  T(0,0, densities).setWeight(pos, T(0,0, densities).getWeight(pos)*S(densities).getWeight(pos)/sum);
	  T(1,1, densities).setWeight(pos, T(1,1, densities).getWeight(pos)*S(densities).getWeight(pos)/sum);
	  T(2,2, densities).setWeight(pos, T(2,2, densities).getWeight(pos)*S(densities).getWeight(pos)/sum);
	}
    }  
  cout << "/////  Finished.  " << endl;
  cout << "///////////////////////////////////////////////////////////" << endl;
}
/*
// This is the old, explicit interpolation. It is now directly incorporated. 
void FMT::interpolate(Point &p, Point &p0, Point &p1, double coeffs[8])
{
  double dx = (p.X() - p0.X())/(p1.X()-p0.X());
  double dy = (p.Y() - p0.Y())/(p1.Y()-p0.Y());
  double dz = (p.Z() - p0.Z())/(p1.Z()-p0.Z());

  coeffs[0] = (1-dx)*(1-dy)*(1-dz); // 000
  coeffs[1] =    dx *(1-dy)*(1-dz); // 100
  coeffs[2] = (1-dx)*   dy *(1-dz); // 010
  coeffs[3] = (1-dx)*(1-dy)*   dz ; // 001
  coeffs[4] =    dx *   dy *(1-dz); // 110
  coeffs[5] =    dx *(1-dy)*   dz ; // 101
  coeffs[6] = (1-dx)*   dy *   dz ; // 011
  coeffs[7] =    dx *   dy *   dz ; // 111
}
*/






double FMT::dPHI(long i)
{
  // the weighted densities for the lattice site under consideration (i).
  // The basic densities are the scalare eta and s, the vector v and the tensor T.
  // Note that s0,s1 etc only differ in constant prefactors if there is only one species present.
  double eta = 0.0;
  double s0  = 0.0;
  double s1  = 0.0;
  double s2  = 0.0;

  double v1[3] = {0.0,0.0,0.0};
  double v2[3] = {0.0,0.0,0.0};
  double T0[3][3] = {{0.0,0.0,0.0},
		 {0.0,0.0,0.0},
		 {0.0,0.0,0.0}};

  // Collect the contributions to the various weighted densities at lattice position i
  // hsd1_ has to do with mixtures (this has not been fully implementd. hsd1_ > 0 is a flag for whether or not there is a second species).
  addWeightedDensityContributions(i, d0_,hsd0_,eta,s0,s1,s2,v1,v2,T0);
  if(hsd1_ > 0)
    addWeightedDensityContributions(i, d1_,hsd1_,eta,s0,s1,s2,v1,v2,T0);

  // touch up the normalization of the tensor density to account for any numerical errors
  double ss = T0[0][0]+T0[1][1]+T0[2][2];

  T0[0][0] +=  (s2-ss)/3;
  T0[1][1] +=  (s2-ss)/3;
  T0[2][2] +=  (s2-ss)/3;

  // Calculate some useful quantities: v dot T, T dot T etc
  double vT[3] = { T0[0][0]*v2[0] +T0[0][1]*v2[1] +T0[0][2]*v2[2],
		   T0[1][0]*v2[0] +T0[1][1]*v2[1] +T0[1][2]*v2[2],
		   T0[2][0]*v2[0] +T0[2][1]*v2[1] +T0[2][2]*v2[2]};
  double TT[3][3] =  { {T0[0][0]*T0[0][0] +T0[0][1]*T0[1][0] +T0[0][2]*T0[2][0],
		     T0[1][0]*T0[0][0] +T0[1][1]*T0[1][0] +T0[1][2]*T0[2][0],
		     T0[2][0]*T0[0][0] +T0[2][1]*T0[1][0] +T0[2][2]*T0[2][0]},

		    {T0[0][0]*T0[0][1] +T0[0][1]*T0[1][1] +T0[0][2]*T0[2][1],
		     T0[1][0]*T0[0][1] +T0[1][1]*T0[1][1] +T0[1][2]*T0[2][1],
		     T0[2][0]*T0[0][1] +T0[2][1]*T0[1][1] +T0[2][2]*T0[2][1]},

		    {T0[0][0]*T0[0][2] +T0[0][1]*T0[1][2] +T0[0][2]*T0[2][2],
		     T0[1][0]*T0[0][2] +T0[1][1]*T0[1][2] +T0[1][2]*T0[2][2],
		     T0[2][0]*T0[0][2] +T0[2][1]*T0[1][2] +T0[2][2]*T0[2][2]}};

  double v1_v2  = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]; //dot(v1,v2);
  double v2_v2  = v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]; //dot(v1,v2);
  double vTv    = v2[0]*vT[0]+v2[1]*vT[1]+v2[2]*vT[2]; //dot(v1,v2);
  double T2     = TT[0][0]+TT[1][1]+TT[2][2]; //trace(TT);
  double T3     = TT[0][0]*T0[0][0]+TT[0][1]*T0[1][0]+TT[0][2]*T0[2][0]
    +TT[1][0]*T0[0][1]+TT[1][1]*T0[1][1]+TT[1][2]*T0[2][1]
    +TT[2][0]*T0[0][2]+TT[2][1]*T0[1][2]+TT[2][2]*T0[2][2];

  double fT = 0;
  double fT1 = 0;

  // These are the eta-dependendent cofactors that lie at the heart of FMT
  double f1 = f1_(eta);
  double f2 = f2_(eta);
  double f3 = f3_(eta);

  // Now, construct the local value of phi.
  double phi = 0;
  phi -= (1/M_PI)*s0*f1; //log(1-eta);
  phi += (1/(2*M_PI))*(s1*s2-v1_v2)*f2;
  phi += Phi3(s2,v2_v2,vTv,T2,T3)*f3; 


  if(etaMax_ > 1.0)
    if(eta > 0.5 && 1-eta < 0.0)
      {
	throw Eta_Too_Large_Exception();
      }

  // Also add in the contributions to the derivative of phi (at lattice site i) wrt the various weighted densities
  // (part of the chain-rule evaluation of dPhi/drho(j) = dPhi/deta(i) * deta(i)/drho(j) + ...)  
  add_dPhi_d_WeightedDensity(i, d0_,hsd0_,eta,s0,s1,s2,v1,v2,T0,v1_v2,v2_v2,vTv,T2,T3,f2,f3, vT, TT);
  if(hsd1_ > 0) add_dPhi_d_WeightedDensity(i, d1_,hsd1_,eta,s0,s1,s2,v1,v2,T0,v1_v2,v2_v2,vTv,T2,T3,f2,f3, vT, TT);

  return phi;
}

void FMT::add_dPhi_d_WeightedDensity(int i, vector<FMT_Weighted_Density> &dd, double hsd, double eta, double s0, double s1, double s2, double v1[], double v2[], double T0[3][3], double v1_v2, double v2_v2, double vTv, double T2, double T3, double f2, double f3, double vT[], double TT[3][3])
{
  double f1 = f1_(eta);
  
  double f1p = f1p_(eta);
  double f2p = f2p_(eta);
  double f3p = f3p_(eta);

  double dPhi_dEta = 0;
  dPhi_dEta -= (1/M_PI)*s0*f1p; ///(1-eta);
  dPhi_dEta += (1/(2*M_PI))*(s1*s2-v1_v2)*f2p;
  dPhi_dEta += Phi3(s2,v2_v2,vTv,T2,T3)*f3p; 

  double dPhi_dS2 = 0;
  dPhi_dS2 += -(1/M_PI)*f1/(hsd*hsd);
  dPhi_dS2 += (1/(2*M_PI))*((s2/hsd)+s1)*f2;
  dPhi_dS2 += dPhi3_dS2(s2,v2_v2,vTv,T2,T3)*f3; //(1.0/(8*M_PI))*(-v2_v2+T2)*f3; 

  double dPhi_dV0[3];

  for(int k=0;k<3;k++)
    {
      dPhi_dV0[k] = 0.0;
      dPhi_dV0[k] += (1/(2*M_PI))*(-v1[k]-v2[k]/hsd)*f2;
      dPhi_dV0[k] += dPhi3_dV2(k, s2, v2_v2, v2, vT)*f3; //(1.0/(8*M_PI))*(2*vT[k]-2*s2*v2[k])*f3;
    }

  double dPhi_dT[3][3];

  for(int k=0;k<3;k++)
    for(int j=0;j<3;j++)
      dPhi_dT[j][k] = dPhi3_dT(j,k,s2,v2,T0, TT)*f3; // (1.0/(8*M_PI))*(v2[j]*v2[k]-3*TT(j,k)+2*s2*T0(j,k))*f3; 


  Eta(dd).Set_dPhi(i,dPhi_dEta);
  S(dd).Set_dPhi(i,dPhi_dS2);

  for(int j=0;j<3;j++)
    {
      V(j,dd).Set_dPhi(i,dPhi_dV0[j]);
      for(int k=j;k<3;k++)	   
	T(j,k,dd).Set_dPhi(i,(j ==k ? 1 : 2)*dPhi_dT[j][k]); // taking account that we only use half the entries
    }

}

// This collects the contributions to the weighted density for each species.
void FMT::addWeightedDensityContributions(int i, vector<FMT_Weighted_Density> &dd, double hsd,double &eta,double &s0, double &s1, double &s2,double v1[], double v2[], double T0[3][3])
{
  eta += Eta(dd).r(i);

  s0 += S(dd).r(i)/(hsd*hsd) ;
  s1 += S(dd).r(i)/hsd;
  s2 += S(dd).r(i);

  for(int j=0;j<3;j++)
    {
      v1[j] += V(j,dd).r(i)/hsd;
      v2[j] += V(j,dd).r(i);
      for(int k=0;k<3;k++)
	T0[j][k] += T(j,k,dd).r(i);
    }
}

double FMT::calculateFreeEnergy(Density& density)
{
  long   Ntot = density.Ntot();
  double dV   = density.dV();

  // Compute the FFT of density 
  density.doFFT();

  // reference to Fourier-space array of density
  const DFT_Vec_Complex &rho_k = density.getDK();

  // This does the convolution of the density and the weight for each weighted density after which it converts back to real space 
  // ( so this computes the weighted densities n(r) = int w(r-r')rho(r')dr'). The results are all stored in parts of FMT_Weighted_Density
  for(FMT_Weighted_Density &d: d0_)
    d.convolute(rho_k);

  // Now compute the free energy. Here we loop over all lattice sites and compute Phi(r_i) for each one. This presupposes that we did the convolution above. 
  
  double F = 0;
  int chunk = Ntot/20;
  long i;

  // There is a problem throwing exceptions from an OMP loop - I think that only the affected thread stops and the others continue.
  // So we eat the exception, remember it, and rethrow it when all loops are finished.
  bool hadCatch = false;

  #pragma omp parallel for		\
    shared( chunk    )				\
    private(i)					\
    schedule(static,chunk)			\
    reduction(+:F)
  for(i=0;i<Ntot;i++)
    {
      try {
	F += dPHI(i);
      } catch( Eta_Too_Large_Exception &e) {
	hadCatch = true;
      }
    }

  // rethrow exception if it occurred: this messiness is do to the parallel evaluation. 
  if(hadCatch) 
    throw   Eta_Too_Large_Exception();
  
  return F;
}

//HERE: this cannot be right for a mixture ... TBD

// Calculate dF[i] = dPhi/drho(i)
//                 = sum_j dV * dPhi(j)/drho(i)
//                 = sum_j dV * sum_alpha dPhi(j)/dn_{alpha}(j)  dn_{alpha}(j)/drho(i)
//                 = sum_j dV * sum_alpha dPhi(j)/dn_{alpha}(j)  w_{alpha}(j,i)
// This is done with convolutions: FMT_Weighted_Density is an array with the index (alpha)
// and holds the weights, w_{alpha}(j,i) and their FFT's AND dPhi(j)/dn_{alpha}(j).
// It therefore FFT's both of these and adds them to dPhi_.Four.
// Once this is done of all alpha, dPhi_ FFT's back to real space and the result is put into dF (with a factor of dV thrown in).

void FMT::calculateFreeEnergyDerivatives(vector<FMT_Weighted_Density> &dd, double dV, DFT_Vec& dF)
{
  dPhi_.Four().zeros();
  
  for(FMT_Weighted_Density &d: dd)
    d.add_to_dPhi(dPhi_.Four());

  dPhi_.do_fourier_2_real();

  dF.set(dPhi_.cReal());

  dF.multBy(dV);
}



static double NSurfactant = -1;

double FMT::calculateFreeEnergyAndDerivatives_fourier_space1(Density& density, DFT_Vec& dF0) 
{
  double F = 0;

  try {
    F = calculateFreeEnergy(density);
  } catch( Eta_Too_Large_Exception &e) {
    throw e;
  }

  double dV   = density.dV();
  calculateFreeEnergyDerivatives(d0_, dV, dF0);

  return F*dV;
};

