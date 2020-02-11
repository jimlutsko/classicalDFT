#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>

#include <gsl/gsl_integration.h>

using namespace std;

#ifdef USE_MPI
#include <mpi.h>  
#endif

#include "Species.h"
#include "myColor.h"

int Species::SequenceNumber_ = 0;

FMT_Species::FMT_Species(Density& density, double hsd, double mu, int seq): Species(density,mu,seq), hsd_(hsd), d_(11)
{
  long Nx = density_.Nx();
  long Ny = density_.Ny();
  long Nz = density_.Nz();

  for(FMT_Weighted_Density &d: d_)
    d.initialize(Nx, Ny, Nz);
}

/*

void FMT_Species::reset(string& pointsFile)
{
  long Nx = density_.Nx();
  long Ny = density_.Ny();
  long Nz = density_.Nz();

  for(FMT_Weighted_Density &d: d_)
    d.initialize(Nx, Ny, Nz);

  generateWeights(pointsFile);

  for(FMT_Weighted_Density &d: d_)
    d.transformWeights();  
}
*/

FMT_Species_Numeric::FMT_Species_Numeric(Density& density, double hsd, string &pointsFile, double mu, int seq):
  FMT_Species(density,hsd,mu,seq)
{
  long Nx = density_.Nx();
  long Ny = density_.Ny();
  long Nz = density_.Nz();

  stringstream ss1;
  ss1 << "weights_"
      << Nx << "_"
      << Ny << "_"
      << Nz << "_"
      << hsd_ << "_"
      << density_.getDX() << "_"
      << ".dat";
  
  bool readWeights = true;
  
  ifstream in(ss1.str().c_str(), ios::binary);
  if(!in.good())
    {readWeights = false; cout << "\n" <<  "Could not open file with weights ... " << endl;}
  else {
    string buf;
    getline(in,buf);

    stringstream ss2(buf);
    int nx, ny, nz;
    double d, dx;
    ss2 >> nx >> ny >> nz >> d >> dx;

    if(nx != Nx)
      {readWeights = false; cout << "\n" <<  "Mismatch in Nx: expected " << Nx << " but read " << nx <<  endl;}
    if(ny != Ny)
      {readWeights = false; cout << "\n" <<  "Mismatch in Ny: expected " << Ny << " but read " << ny <<  endl;}
    if(nz != Nz)
      {readWeights = false; cout << "\n" <<  "Mismatch in Nz: expected " << Nz << " but read " << nz <<  endl;}
    if(fabs(hsd_-d) > 1e-8*(hsd_+d))
      {readWeights = false; cout << "\n" <<  "Mismatch in hsd: expected " << hsd_ << " but read " << d << endl;}
    if(fabs(density_.getDX()-dx) > 1e-8*(density_.getDX()+dx))
      {readWeights = false; cout << "\n" <<  "Mismatch in Dx: expected " << density_.getDX() << " but read " << dx << endl;}

    getline(in,buf);
    stringstream ss3(buf);
    if(ss3.str().compare(pointsFile) != 0)
      {readWeights = false; cout << "\n" <<  "Mismatch in points file: expected " << pointsFile << " but read " << ss3.str() <<  endl;}
  }
  if(readWeights)
    {
      for(auto& s: d_)
	s.load(in);
    } else {
    generateWeights(pointsFile);
  }
  /*
  cout << "Sample: " << endl;
  int counter = 0;
  for(long i=0;i<d_[EI()].Real().size() && counter < 20; i++)
    {
      double s = d_[SI()].getWeight(i);
      if(fabs(s) > 1e-10)
	{
	  cout << i << " " << d_[EI()].getWeight(i) << " " << d_[SI()].getWeight(i) << " " << d_[VI(0)].getWeight(i) << " " << d_[VI(1)].getWeight(i) << " " << d_[VI(2)].getWeight(i) << endl;
	  counter++;
	}
    }
  exit(0);
  */
  for(FMT_Weighted_Density &d: d_)
    d.transformWeights();
}


// The idea here is as follows. Consider the s weighted density which is just an integral over a sphere of radius hsd/2.
// In principle, we will need the volume integral
//    s(r_i) = int w(|r_i - r'|) rho(r') delta(r'-(hsd/2)) dr'
// The first step is to approximate as some discrete sum on the spherical shell:
//    s(r_i) ~ sum_{k} w(|r_i - r'_{k}*(hsd/2)|) rho(r'_{k}*(hsd/2)) ww_{k}
// where I am assuming that the points r'_{k} are given for a unit sphere.
// Next, we evaluate rho using trilinear interpolation. This will involve a sum over the 8 corners of the cubic cell
// containing the integration point (i.e. its nearest neighbors on the density_)
// with corresponding weights (which I call u):
//    s(r_i) ~ sum_{k} w(|r_i - r'_{k}*(hsd/2)|) ww_{k} sum_{l \in nearest neighbors of r'_{k}*(hsd/2)} rho(r_{l}) u_{l}(r'_{k}*(hsd/2))
// or
//     s(r_i) ~ sum_{j} A_{ij} rho(r_{k})
// where
//       A_{ij} = sum_{k} w(|r_i - r'_{k}*(hsd/2)|) ww_{k} u_{l}(r'_{k}*(hsd/2))  
// with the understanding that the weight u_{l} is zero if l is not a nearest neighbor of the point r'_{k}*(hsd/2).
// Finally, A is translationally invariant so we really just need it for i=0. It is the arrays A_{0,j} (i.e. the "weights") which are computed here. 

void FMT_Species_Numeric::generateWeights(string &pointsFile)
{
  long Nx = density_.Nx();
  long Ny = density_.Ny();
  long Nz = density_.Nz();

  double dx = density_.getDX();
  double dy = density_.getDY();
  double dz = density_.getDZ();

  
  // Generate integration points for spherical surface

  cout << endl;
  cout << myColor::GREEN;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << "/////  Generating integration points on sphere " << endl;

  vector < vector<double> > points; 

#ifndef USE_MPI  
  // Read points from file : C++ way
  
  ifstream in(pointsFile.c_str());
  if(!in.good())
    {
      stringstream ss;
      ss << "File " << pointsFile << " could not be opened";      
      throw std::runtime_error(ss.str().c_str());
    }
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
      line.push_back(4*M_PI); // the weight

      double r2 = line[0]*line[0]+line[1]*line[1]+line[2]*line[2];
      double r = sqrt(r2);
      line[0] /= r;
      line[1] /= r;
      line[2] /= r;

      points.push_back(line); // points includes the weights as the fourth entry
    }
  // adjust the weights
  for(vector<double> &point : points)
    point[3] /= points.size();

  cout << "/////  Finished: there are " << points.size() << " points " << endl;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << myColor::RESET << endl;


  // Add up the weights for each point.
  // We throw in all permutations (iperm loop) of the integration points and all reflections (is loop)
  // to establish as much symmetry as possible. Probably unnecessary.    
  cout << endl;
  cout << myColor::GREEN;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << "/////  Generating integration points for sphere volume " << endl;
  cout << myColor::RESET << endl;

  // This is for the radial integral
  int Nr = 128*4;
  double r = 0.5*hsd_; // hard sphere radius
  // Get some Gauss-Lagrange integration points and weights for the radial integral
  gsl_integration_glfixed_table *tr = gsl_integration_glfixed_table_alloc(Nr);

  // Decide whether or not to do the permuations:
  int iperm_max = 6; // 1 or 6
  int is_max = 8; // 1 or 8
  double inorm = iperm_max*is_max;
  
  long count = 0;
  for(int pos=0;pos < points.size(); pos++)
    {
      
      if(pos%1000 == 0) {if(pos > 0) cout << '\r'; cout << "\t" << int(double(pos)*100.0/points.size()) << "% finished: " << pos << " out of " << points.size(); cout.flush();}
      for(int iperm = 0; iperm < iperm_max; iperm++)      // add permutations
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

	  for(int is=0;is<is_max;is++) // add in reflections too
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
	      double ww = points[pos][3]/inorm;

	      double v[3]; // used for the vector and tensor densities
	      v[0] = x/r;
	      v[1] = y/r;
	      v[2] = z/r;

	      // Find "smallest" corner of the cube that contains this point. 
	      int ix0 = int(x/dx);
	      int iy0 = int(y/dy);
	      int iz0 = int(z/dz);

	      if(x < 0) ix0 -=1;
	      if(y < 0) iy0 -=1;
	      if(z < 0) iz0 -=1;

	      double ddx = (x/dx)-ix0;
	      double ddy = (y/dy)-iy0;
	      double ddz = (z/dz)-iz0;

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

		      while(nx < 0) nx += Nx;
		      while(ny < 0) ny += Ny;
		      while(nz < 0) nz += Nz;

		      while(nx >= Nx) nx -= Nx;
		      while(ny >= Ny) ny -= Ny;
		      while(nz >= Nz) nz -= Nz;

		      long pos = nz+Nz*(ny+Ny*nx);
	      
		      d_[SI()].addToWeight(pos,cc*ww*r*r);
	      
		      for(int jv=0;jv<3;jv++)
			{
			  d_[VI(jv)].addToWeight(pos, -cc*ww*r*r*v[jv]);

			  for(int kT=0;kT<=jv;kT++)
			    d_[TI(jv,kT)].addToWeight(pos,cc*ww*r*r*v[jv]*v[kT]);	
			}
		      count++;
		    }

	      // Now add in the radial integral as well ...
	      // This just adds another sum 
	      for(int k=0;k<Nr;k++)
		{
		  double r;
		  double wr;
	      
		  gsl_integration_glfixed_point(0, 0.5*hsd_, k, &r, &wr, tr);
	      
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

		  double ww = points[pos][3]/inorm;
	      
		  // Find the cube that contains this point. 
		  int ix0 = int(x/dx);
		  int iy0 = int(y/dy);
		  int iz0 = int(z/dz);

		  if(x < 0) ix0 -=1;
		  if(y < 0) iy0 -=1;
		  if(z < 0) iz0 -=1;

		  double ddx = (x/dx)-ix0;
		  double ddy = (y/dy)-iy0;
		  double ddz = (z/dz)-iz0;

		  for(int io=0;io<2;io++)
		    for(int jo=0;jo<2;jo++)
		      for(int ko=0;ko<2;ko++)		  
			{

			  int nx = ix0+io;
			  int ny = iy0+jo;
			  int nz = iz0+ko;

			  while(nx < 0) nx += Nx;
			  while(ny < 0) ny += Ny;
			  while(nz < 0) nz += Nz;
			
			  while(nx >= Nx) nx -= Nx;
			  while(ny >= Ny) ny -= Ny;
			  while(nz >= Nz) nz -= Nz;			

			  long pos = nz+Nz*(ny+Ny*nx);

			  double cc = (io == 0 ? 1-ddx : ddx)*(jo == 0 ? 1-ddy : ddy)*(ko == 0 ? 1-ddz : ddz); // trilinear interpolation
			  d_[EI()].addToWeight(pos,cc*ww*wr*r*r);      
			}
		}
	    }
	}
    }
  cout << endl;
  // This is some overly-cautious adjusting to assure that some exact relations are satsified. (Numerical noise may have violated them slightly).
  for(long pos = 0;pos <Nx*Ny*Nz; pos++) 
    {
      double sum = d_[TI(0,0)].getWeight(pos)+d_[TI(1,1)].getWeight(pos)+d_[TI(2,2)].getWeight(pos);
      if(fabs(sum) > 1e-16)
	{
	  d_[TI(0,0)].setWeight(pos, d_[TI(0,0)].getWeight(pos)*d_[SI()].getWeight(pos)/sum);
	  d_[TI(1,1)].setWeight(pos, d_[TI(1,1)].getWeight(pos)*d_[SI()].getWeight(pos)/sum);
	  d_[TI(2,2)].setWeight(pos, d_[TI(2,2)].getWeight(pos)*d_[SI()].getWeight(pos)/sum);
	}
    }
  cout << myColor::GREEN;
  cout << "/////  Finished.  " << endl;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << myColor::RESET << endl;

  size_t Npos = density_.Ntot();
  
  for(int i=0;i<11;i++)
    {
      size_t nonzero = 0;
      for(size_t pos = 0; pos < Npos; pos++)
	if(fabs(d_[i].getWeight(pos)) > 1e-10)
	  nonzero++;

      cout << "Weight " << i << " has " << nonzero << " non-zero elements out of a total of " << Npos << " or a fraction of " << double(nonzero)/Npos << endl;
    }

  
  /// Dump the weights
  stringstream ss;
  ss << "weights_"
      << Nx << "_"
      << Ny << "_"
      << Nz << "_"
      << hsd_ << "_"
      << density_.getDX() << "_"
      << ".dat";
  
  ofstream of(ss.str().c_str(), ios::binary);

  of.flags (std::ios::scientific);
  of.precision (std::numeric_limits<double>::digits10 + 1);

  of << Nx << " " << Ny << " " << Nz << " " << hsd_ << " " << density_.getDX() << endl;
  of << pointsFile << endl;
  
  for(auto& s: d_)
    s.dump(of);  
}
///////////////////////////////////////////////////////
// Beginning of FMT_Species_Analytic section

FMT_Species_Analytic::FMT_Species_Analytic(Density& density, double hsd, double mu, int seq):
  FMT_Species(density,hsd,mu,seq)
{
  generateWeights();

  for(FMT_Weighted_Density &d: d_)
    d.transformWeights();
}

double J_eta_helper(double R, double Sx, double Sy, double Sz, double c, double b, double a)
{
  double j1 = 0.0;
  j1 += a*a*(3*a*a*a*Sz+20*b*b*b*Sy+30*a*a*b*Sy+30*b*b*c*Sx-40*a*b*Sy*Sz-60*b*c*Sx*Sy)/240.0;
  j1 += R*R*(3*a*Sz*R*R+12*a*b*Sy*Sz-2*a*a*a*Sz-12*a*a*b*Sy)/48.0;

  j1 += -Sx*Sy*(3*(R*R-a*a)*(R*R-a*a)+4*a*Sz*(3*R*R-a*a))*(asin(b/sqrt(R*R-a*a))-M_PI/4)/24.0;
  j1 += Sx*Sy*Sz*R*R*R*asin(a*b/(sqrt(R*R-a*a)*sqrt(R*R-b*b)))/6.0;
  j1 += -Sx*(4*R*R*R*R+8*a*a*b*b-25*Sy*a*a*b+20*a*b*Sy*Sz+8*a*a*a*a-10*Sz*a*a*a-16*R*R*a*a+25*Sz*R*R*a)*sqrt(R*R-a*a-b*b)/120.0;
  return j1;
}

double J_eta(double R, double Sx, double Sy, double Sz, int wx, int wy, int wz)
{
  int va[2] = {0,wz};
  int vb[2] = {0,wy};
  int vc[2] = {0,wx};

  double j = 0.0;

  for(int ia = 0; ia < 2; ia++)
    for(int ib = 0; ib < 2; ib++)
      for(int ic = 0; ic < 2; ic++)
	{
	  double a = Sz - va[ia];
	  double b = Sy - vb[ib];
	  double c = Sx - vc[ic];

	  int sgn = 1;
	  if(ia == 1) sgn *= -1;
	  if(ib == 1) sgn *= -1;
	  if(ic == 1) sgn *= -1;
			  
	  if(R*R < a*a+b*b+c*c) continue;

	  double j0 = 0.0;
	  
	  j0 += (R*R-a*a-b*b-c*c)*(R*R-a*a-b*b-c*c)*(R*R-a*a-b*b-c*c)/48.0;
	  j0 += a*b*c*Sx*Sy*Sz;
	  j0 -= (M_PI/6)*Sx*Sy*Sz*R*R*R;

	  j0 += J_eta_helper(R,Sx,Sy,Sz,c,b,a);
	  j0 += J_eta_helper(R,Sx,Sz,Sy,c,a,b);
	  j0 += J_eta_helper(R,Sy,Sx,Sz,b,c,a);
	  j0 += J_eta_helper(R,Sy,Sz,Sx,b,a,c);
	  j0 += J_eta_helper(R,Sz,Sx,Sy,a,c,b);
	  j0 += J_eta_helper(R,Sz,Sy,Sx,a,b,c);

	  j += sgn * j0;
	}
  
  return j;
}

double J_Vz_helper(double R, double Sx, double Sy, double Sz, double c, double b, double a)
{
  double j1 = 0.0;
  
  j1 += a*(15*a*a*a*Sz+3*a*a*a*a+40*b*b*b*Sy+80*a*a*b*Sy+60*b*b*c*Sx-120*a*b*Sy*Sz-120*b*c*Sx*Sy)/240.0;
  j1 += b*b*(60*a*a*Sz+20*a*a*a+30*b*b*Sz+30*b*b*a-40*b*Sy*Sz-40*a*b*Sy-60*c*Sx*Sz)/240.0;
  j1 += c*c*(30*b*b*Sz+30*b*b*a-60*b*Sy*Sz-60*a*b*Sy)/240.0;

  j1 += R*R*(3*(4*b*Sy+R*R)*(a+Sz)-6*a*a*Sz-2*a*a*a-24*a*b*Sy)/48.0;
  j1 += -R*R*b*(b-Sy)*(a+Sz)/4;
  
  j1 += -Sx*Sy*(3*(R*R-a*a)*(R*R-a*a)+4*a*Sz*(3*R*R-a*a))*(a*b/((R*R-a*a)*sqrt(R*R-a*a-b*b)))/24.0;
  j1 += -Sx*Sy*(2*a*a*a+3*Sz*(R*R-a*a))*(asin(b/sqrt(R*R-a*a))-M_PI/4)/6.0;
  j1 += -Sx*(3*(R*R-b*b)*(R*R-b*b)+4*b*Sy*(3*R*R-b*b))*(asin(a/sqrt(R*R-b*b))-M_PI/4)/24.0;
  j1 += -Sx*Sz*(3*(R*R-b*b)*(R*R-b*b)+4*b*Sy*(3*R*R-b*b))*(1.0/sqrt(R*R-a*a-b*b))/24.0;
  j1 += -Sy*(3*(R*R-c*c)*(R*R-c*c)+4*c*Sx*(3*R*R-c*c))*(asin(b/sqrt(R*R-c*c))-M_PI/4)/24.0;
  
  j1 += -Sx*(16*a*b*b-30*Sy*a*b+20*b*Sy*Sz+22*a*a*a-30*Sz*a*a-7*R*R*a+25*Sz*R*R)*sqrt(R*R-a*a-b*b)/120.0;
  j1 += -Sx*(-9*a*b*b-25*Sz*b*b+20*b*Sy*Sz+20*a*b*Sy)*sqrt(R*R-a*a-b*b)/120.0;
  j1 += -(4*R*R*R*R+8*c*c*b*b-25*Sy*c*c*b+20*c*b*Sy*Sx+8*c*c*c*c-10*Sx*c*c*c-16*R*R*c*c+25*Sx*R*R*c)*sqrt(R*R-c*c-b*b)/120.0;  
  j1 += Sx*(4*R*R*R*R+8*a*a*b*b-25*Sy*a*a*b+20*a*b*Sy*Sz+8*a*a*a*a-10*Sz*a*a*a-16*R*R*a*a+25*Sz*R*R*a)*(a/sqrt(R*R-a*a-b*b))/120.0;
  j1 += Sx*(4*R*R*R*R+8*a*a*b*b-25*Sz*a*b*b+20*a*b*Sy*Sz+8*b*b*b*b-10*Sy*b*b*b-16*R*R*b*b+25*Sy*R*R*b)*(a/sqrt(R*R-a*a-b*b))/120.0;
  
  return j1;
}

double J_Vz(double R, double Sx, double Sy, double Sz, int wx, int wy, int wz)
{
  int va[2] = {0,wz};
  int vb[2] = {0,wy};
  int vc[2] = {0,wx};

  double j = 0.0;

  for(int ia = 0; ia < 2; ia++)
    for(int ib = 0; ib < 2; ib++)
      for(int ic = 0; ic < 2; ic++)
	{
	  double a = Sz - va[ia];
	  double b = Sy - vb[ib];
	  double c = Sx - vc[ic];

	  int sgn = 1;
	  if(ia == 1) sgn *= -1;
	  if(ib == 1) sgn *= -1;
	  if(ic == 1) sgn *= -1;
			  
	  if(R*R < a*a+b*b+c*c) continue;

	  double j0 = 0.0;
	  
	  j0 += -a*(R*R-a*a-b*b-c*c)*(R*R-a*a-b*b-c*c)/8.0;
	  j0 += b*c*Sx*Sy*Sz+a*b*c*Sx*Sy;
	  j0 += -(M_PI/6)*Sx*Sy*R*R*R;

	  j0 += Sx*Sy*R*R*R*(asin(a*b/(sqrt(R*R-a*a)*sqrt(R*R-b*b)))
			     +asin(a*c/(sqrt(R*R-a*a)*sqrt(R*R-c*c)))
			     +asin(b*c/(sqrt(R*R-b*b)*sqrt(R*R-c*c))))/3.0;
	  j0 += Sx*Sy*Sz*R*R*R*R*((b/((R*R-a*a)*sqrt(R*R-a*a-b*b)))
				  +(c/((R*R-a*a)*sqrt(R*R-a*a-c*c))))/3.0;
	  
	  j0 += J_Vz_helper(R,Sx,Sy,Sz,c,b,a);
	  j0 += J_Vz_helper(R,Sy,Sx,Sz,b,c,a);

	  j += sgn * j0;
	}
  
  return j;
}

double J_s_helper(double R, double Sx, double Sy, double Sz, double c, double b, double a)
{
  double j1 = 0.0;
  j1 += R*(6*a*Sz*R*R+12*a*b*Sy*Sz-2*a*a*a*Sz-12*a*a*b*Sy)/24.0;
  j1 += -Sx*Sy*R*(R*R-a*a+2*a*Sz)*(asin(b/sqrt(R*R-a*a))-M_PI/4)/2.0;
  j1 += -Sx*Sy*(3*(R*R-a*a)*(R*R-a*a)+4*a*Sz*(3*R*R-a*a))*(-R*b/((R*R-a*a)*sqrt(R*R-a*a-b*b)))/24.0;
  j1 += Sx*Sy*Sz*R*R*asin(a*b/(sqrt(R*R-a*a)*sqrt(R*R-b*b)))/2.0;
  j1 += Sx*Sy*Sz*R*R*R*(a*b/((R*R-a*a)*(R*R-b*b)))*((a*a+b*b-2*R*R)/sqrt(R*R-a*a-b*b))/6.0;
  j1 += -Sx*(16*R*R*R-32*R*a*a+50*Sz*R*a)*sqrt(R*R-a*a-b*b)/120.0;
  j1 += -Sx*(4*R*R*R*R+8*a*a*b*b-25*Sy*a*a*b+20*a*b*Sy*Sz+8*a*a*a*a-10*Sz*a*a*a-16*R*R*a*a+25*Sz*R*R*a)*(R/sqrt(R*R-a*a-b*b))/120.0;
  return j1;
}

double J_s(double R, double Sx, double Sy, double Sz, int wx, int wy, int wz)
{
  int va[2] = {0,wz};
  int vb[2] = {0,wy};
  int vc[2] = {0,wx};

  double j = 0.0;

  for(int ia = 0; ia < 2; ia++)
    for(int ib = 0; ib < 2; ib++)
      for(int ic = 0; ic < 2; ic++)
	{
	  double a = Sz - va[ia];
	  double b = Sy - vb[ib];
	  double c = Sx - vc[ic];

	  int sgn = 1;
	  if(ia == 1) sgn *= -1;
	  if(ib == 1) sgn *= -1;
	  if(ic == 1) sgn *= -1;
			  
	  if(R*R < a*a+b*b+c*c) continue;

	  double j0 = 0.0;
	  
	  j0 += R*(R*R-a*a-b*b-c*c)*(R*R-a*a-b*b-c*c)/8.0;
	  j0 -= (M_PI/2)*Sx*Sy*Sz*R*R;

	  j0 += J_s_helper(R,Sx,Sy,Sz,c,b,a);
	  j0 += J_s_helper(R,Sx,Sz,Sy,c,a,b);
	  j0 += J_s_helper(R,Sy,Sx,Sz,b,c,a);
	  j0 += J_s_helper(R,Sy,Sz,Sx,b,a,c);
	  j0 += J_s_helper(R,Sz,Sx,Sy,a,c,b);
	  j0 += J_s_helper(R,Sz,Sy,Sx,a,b,c);

	  j += sgn * j0;
	}
  
  return j;
}

void FMT_Species_Analytic::generateWeights()
{
  double dx = density_.getDX();
  double dy = density_.getDY();
  double dz = density_.getDZ();

  double dV = dx*dy*dz;
  double dS = dx*dy;

  double hsr = hsd_/2; // the hard-sphere radius
  
  // This saves having to a code a useless special case
  if(hsd_ < dx)
    throw std::runtime_error("hsd is less than the lattice spacing ... aborting");
  
  int Sx_max = 2+int(hsr/dx);
  int Sy_max = 2+int(hsr/dy);
  int Sz_max = 2+int(hsr/dz);

  long pmax = (Sx_max+1)*(Sy_max+1)*(Sz_max+1);
  
  int v[2] = {-1,1};

  cout << endl;
  cout << myColor::GREEN;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << "/////  Generating weights using analytic formulae" << endl;
  cout << myColor::RESET << endl;

  long counter = 0;
  
  for(int Sx = 0; Sx <= Sx_max; Sx++)
    for(int Sy = 0; Sy <= Sy_max; Sy++)
      for(int Sz = 0; Sz <= Sz_max; Sz++)
	{
	  counter++;
	  if(counter%1000 == 0) {if(counter > 0) cout << '\r'; cout << "\t" << int(double(counter)*100.0/pmax) << "% finished: " << counter << " out of " << pmax; cout.flush();}

	  double R2_min = (Sx-1)*(Sx-1)+(Sy-1)*(Sy-1)+(Sz-1)*(Sz-1);
	  double R2_max = (Sx+1)*(Sx+1)+(Sy+1)*(Sy+1)+(Sz+1)*(Sz+1);
	  
	  double w_eta = 0.0;
	  double w_s   = 0.0;
	  double w_v[3] = {0.0,0.0,0.0};

	  // The weight for point Sx,Sy,Sz has contributions from all adjoining cells
	  // The furthest corner is Sx+1,Sy+1,Sz+1 and if this less than hsd*hsd, the volume weights are 1, and the surface weights are zero
	  // else, the nearest corner is Sx-1,Sy-1,Sz-1 and if this is less than hsd*hsd, then the boundary is between these limits and we must compute
	  // else, all hsd boundary is less than the nearest corner and all weights are zero.

	  double R = hsr/dx;
	  
	  if(R*R > R2_max) {w_eta = dV;}
	  else if(R*R > R2_min)
	    for(int ix:v)
	      for(int iy:v)
		for(int iz:v)
		  {
		    w_eta  += dV*J_eta(R,Sx+ix,Sy+iy,Sz+iz,ix,iy,iz);
		    w_s    += dS*  J_s(R,Sx+ix,Sy+iy,Sz+iz,ix,iy,iz);
		    w_v[2] += dS* J_Vz(R,Sx+ix,Sy+iy,Sz+iz,ix,iy,iz);
		    w_v[1] += dS* J_Vz(R,Sx+ix,Sz+iz,Sy+iy,iy,ix,iz);
		    w_v[0] += dS* J_Vz(R,Sz+iz,Sy+iy,Sx+ix,iz,iy,ix);
		  }
	  
	  // Add in for all symmetries
	  for(int ix = 0; ix < (Sx == 0 ? 1 : 2); ix++)
	    for(int iy = 0; iy < (Sy == 0 ? 1 : 2); iy++)
	      for(int iz = 0; iz < (Sz == 0 ? 1 : 2); iz++)
		{		  
		  long pos = density_.get_PBC_Pos((1-2*ix)*Sx,(1-2*iy)*Sy,(1-2*iz)*Sz);	  
		  d_[EI()].addToWeight(pos,w_eta);
		  d_[SI()].addToWeight(pos,w_s);
		  for(int iv = 0;iv < 3;iv++)
		    d_[VI(iv)].addToWeight(pos,(iv == 0 ? (1-2*ix) : (iv == 1 ? (1-2*iy) : (1-2*iz)))*w_v[iv]);
		}
	}
  cout << endl;
  cout << myColor::GREEN;
  cout << "/////  Finished.  " << endl;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << myColor::RESET << endl;

}





/*
VDW_Species::VDW_Species(Density& density, double hsd, string &pointsFile, Potential1& potential, double kT)
  : FMT_Species(density,hsd,pointsFile), potential_(potential), a_vdw_(0)
{
  // The lattice
  long Nx = density.Nx();
  long Ny = density.Ny();
  long Nz = density.Nz();

  double dx = density.getDX();
  double dy = density.getDY();
  double dz = density.getDZ();
  
  // Set up the mean field potential
  w_att_.initialize(Nx,Ny,Nz);

  for(int nx = 0;nx<Nx;nx++)
    for(int ny = 0;ny<Ny;ny++)
      for(int nz = 0;nz<Nz;nz++)
	{
	  long pos = nz+Nz*(ny+Ny*nx);
	  
	  double x = nx*dx;
	  double y = ny*dy;
	  double z = nz*dz;

	  if(nx > Nx/2) x -= Nx*dx;
	  if(ny > Ny/2) y -= Ny*dy;
	  if(nz > Nz/2) z -= Nz*dz;

	  double r2 = x*x+y*y+z*z;
	  w_att_.Real().IncrementBy(pos,potential.Watt(sqrt(r2))/kT); 
 	}

  // Set the parameters in the VDW object  
  // Do this BEFORE the FFT which may corrupt the real-space part
  a_vdw_ = 0.5*w_att_.Real().accu()*dx*dy*dz;

  // Now save the FFT of the field  
  w_att_.do_real_2_fourier();
}
*/
