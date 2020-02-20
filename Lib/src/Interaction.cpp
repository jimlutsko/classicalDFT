/* This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
 * To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
 *
 * Author: James F. Lutsko
 * www.lutsko.com
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

#include <gsl/gsl_integration.h>

using namespace std;

#ifdef USE_OMP
#include <omp.h>
#endif

#include "Interaction.h"

void Interaction_Base::initialize()
{
  if(initialized_) return; // no need to do twice

  const Density &density = s1_.getDensity();
  long Nx = density.Nx();
  long Ny = density.Ny();
  long Nz = density.Nz();

  w_att_.initialize(Nx,Ny,Nz);      
  a_vdw_ = 0.0;
  
  // First, try to read the weights from a file
  stringstream ss1;
  
  if(!readWeights(ss1))
    generateWeights(v_, ss1, log_);
    
  // Introduce the temperature
  w_att_.Real().MultBy(1.0/kT_);
  a_vdw_ /= kT_;
  
  // Now generate the FFT of the field  
  w_att_.do_real_2_fourier();

  initialized_ = true;
}    

bool Interaction_Base::readWeights(stringstream &ss1)
{
  const Density &density = s1_.getDensity();
  long Nx = density.Nx();
  long Ny = density.Ny();
  long Nz = density.Nz();
        
  ss1 << "weights_"
      << s1_.getSequenceNumber() << "_"
      << s2_.getSequenceNumber() << "_"
      << Nx << "_"
      << Ny << "_"
      << Nz << "_"
      << v_.getIdentifier() << "_";
  ss1 << ".dat";    
  
  bool readWeights = true;
    
  ifstream in(ss1.str().c_str(), ios::binary);
  if(!in.good())
    {
      readWeights = false;
      cout << myColor::GREEN << endl;
      cout << "///////////////////////////////////////////////////////////" << endl;
      cout << myColor::GREEN  << "\n" <<  "Could not open file with potential kernal: it will be generated" << myColor::RESET << endl;
    } else {
    string buf;
    getline(in,buf);

    stringstream ss2(buf);
    int nx, ny, nz;
    double dx;
    ss2 >> nx >> ny >> nz >> dx;

    if(nx != Nx)
      {readWeights = false; cout << "\n" <<  "Mismatch in Nx: expected " << Nx << " but read " << nx <<  endl;}
    if(ny != Ny)
      {readWeights = false; cout << "\n" <<  "Mismatch in Ny: expected " << Ny << " but read " << ny <<  endl;}
    if(nz != Nz)
      {readWeights = false; cout << "\n" <<  "Mismatch in Nz: expected " << Nz << " but read " << nz <<  endl;}
    if(fabs(density.getDX()-dx) > 1e-8*(density.getDX()+dx))
      {readWeights = false; cout << "\n" <<  "Mismatch in Dx: generating weights: expected " << density.getDX() << " but read " << dx << endl;}      

    getline(in,buf);
    stringstream ss4(buf);      
    string identifier;
    ss4 >> identifier;
    string pot = v_.getIdentifier();
    if(identifier.compare(pot) != 0)
      {readWeights = false; cout << "\n" <<  "Mismatch in potential: expected " << pot << " but read " << identifier <<  endl;}
      
    readWeights = checkWeightsFile(in);
  }
  if(readWeights)
    {
      w_att_.Real().load(in);
      a_vdw_ = w_att_.Real().accu();
    } else {
    ofstream of(ss1.str().c_str(), ios::binary);

      of.flags (std::ios::scientific);
      of.precision (std::numeric_limits<double>::digits10 + 1);

      of << Nx << " " << Ny << " " << Nz << " " << density.getDX() << endl;
      of << v_.getIdentifier() << endl;
      of.close();
  }
  return readWeights;
}

  // Note that the matrix w already contains a factor of dV
double Interaction_Base::getInteractionEnergyAndForces()
{
  if(!initialized_)
    initialize();

  const Density &density1 = s1_.getDensity();
    
  long Ntot = density1.Ntot();
  double dV = density1.dV();

  DFT_FFT v(density1.Nx(), density1.Ny(), density1.Nz());
  v.zeros();
  double E = 0;
    
  if(s1_.getSequenceNumber() == s2_.getSequenceNumber())
    {
      v.Four().Schur(density1.getDK(),w_att_.Four(),false);
      v.do_fourier_2_real();
      v.Real().MultBy(dV/Ntot);
      s1_.addToForce(v.Real());
      E = 0.5*density1.getInteractionEnergy(v.Real());
    } else {
    v.Four().Schur(density1.getDK(),w_att_.Four());
    v.Four().MultBy(0.5*dV/Ntot);
    v.do_fourier_2_real(); 
    s2_.addToForce(v.Real());

    const Density &density2 = s2_.getDensity();

    v.Four().Schur(density2.getDK(),w_att_.Four());
    v.Four().MultBy(0.5*dV*dV/Ntot);
    v.do_fourier_2_real(); 
    s1_.addToForce(v.Real());

    E = 0.5*density1.getInteractionEnergy(v.Real());	
  }
  return E;
}

double Interaction_Base::checkCalc(int jx, int jy, int jz)
{
  const Density &density1 = s1_.getDensity();
    
  long Ntot = density1.Ntot();
  double dV = density1.dV();
    
  int Nx = density1.Nx();
  int Ny = density1.Ny();
  int Nz = density1.Nz();

  long sj = density1.pos(jx,jy,jz);
  double dd = 0.0;    
  for(int ix=0;ix<Nx; ix++)
    for(int iy=0;iy<Ny; iy++)
      for(int iz=0;iz<Nz; iz++)
	{
	  int kx = jx-ix;  while(kx < 0) kx += Nx; while (kx > Nx) kx -= Nx;
	  int ky = jy-iy;  while(ky < 0) ky += Ny; while (ky > Ny) ky -= Ny;
	  int kz = jz-iz;  while(kz < 0) kz += Nz; while (kz > Nz) kz -= Nz;
		    
	  long si = density1.pos(ix,iy,iz);
	  long sk = density1.pos(kx,ky,kz);
		    
	  dd += w_att_.Real().get(sk)*density1.getDensity(si);
	}
  return dd*dV*dV;
}



bool Interaction::checkWeightsFile(ifstream &in)
{
  string buf;
  getline(in,buf);
  stringstream ss3(buf); // probably could just use buf ...

  bool readWeights = (ss3.str().compare(pointsFile_) == 0);
    
  if(!readWeights)
    cout << "\n" <<  "Mismatch in points file: expected " << pointsFile_ << " but read " << ss3.str() <<  endl;
    
  return readWeights;
}
  
void Interaction::generateWeights(Potential1 &v, stringstream &ss, Log& log)
{    
  const Density &density = s1_.getDensity();
      
  // The lattice
  long Nx = density.Nx();
  long Ny = density.Ny();
  long Nz = density.Nz();

  double dx = density.getDX();
  double dy = density.getDY();
  double dz = density.getDZ();
  
  // Set up the mean field potential
  // We need to take into account the whole contribution of the potential out to its cutoff of Rc.
  // This may mean going beyond nearest neighbors in certain conditions.
  // We also compute the vdw parameter at the same time.
      
  double Rc = v.getRcut();

  // Generate integration points for spherical surface

  cout << endl;
  cout << myColor::GREEN;

  vector < vector<double> > points; 

#ifndef USE_MPI  
  // Read points from file : C++ way
  
  ifstream in(pointsFile_.c_str());
  if(!in.good())
    {
      stringstream ss;
      ss << "File " << pointsFile_ << " could not be opened";
      throw std::runtime_error(ss.str().c_str());
    }
#else
  // Read points from file : C code for use with MPI functions
  MPI_Status    status;
  MPI_File      fh;

  int error = MPI_File_open(MPI_COMM_WORLD, pointsFile_.c_str(),
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

  cout << myColor::RESET << endl;


  // Add up the weights for each point.
  // We throw in all permutations (iperm loop) of the integration points and all reflections (is loop)
  // to establish as much symmetry as possible. Probably unnecessary.    
  // Get some Gauss-Lagrange integration points and weights for the radial integral
  int Nr = 128*4;
  gsl_integration_glfixed_table *tr = gsl_integration_glfixed_table_alloc(Nr);

  // Decide whether or not to do the permuations:
  int iperm_max = 6; // 1 or 6
  int is_max = 8; // 1 or 8
  double inorm = iperm_max*is_max;
    
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

	  for(int is=0;is<is_max;is++) // add in reflections too
	    {
	      for(int k=0;k<Nr;k++)
		{
		  double r;
		  double wr;
	      
		  gsl_integration_glfixed_point(0, Rc, k, &r, &wr, tr);
		  
		  double watt = v.Watt(r); 
		  
		  double x = r*points[pos][0];
		  double y = r*points[pos][1];
		  double z = r*points[pos][2];
		  
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

			  double w = cc*ww*wr*r*r*watt;     
			  a_vdw_ += w;
			  w_att_.Real().IncrementBy(pos,w);
			}
		}
	    }
	}
    }
  cout << myColor::GREEN;
  cout << "/////  Finished.  " << endl;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << myColor::RESET << endl;
   
  /// Dump the weights
  ofstream of(ss.str().c_str(), ios::binary | ios::app);  
  of << pointsFile_ << endl;
  w_att_.Real().save(of);
}

void Interaction_Full::generateWeights(Potential1 &v, stringstream &ss, Log& log)
{    
  const Density &density = s1_.getDensity();
      
  // The lattice
  long Nx = density.Nx();
  long Ny = density.Ny();
  long Nz = density.Nz();

  double dx = density.getDX();
  double dy = density.getDY();
  double dz = density.getDZ();
  
  // Set up the mean field potential
  // We need to take into account the whole contribution of the potential out to its cutoff of Rc.
  // This may mean going beyond nearest neighbors in certain conditions.
  // We also compute the vdw parameter at the same time.
      
  double Rc = v.getRcut();

  int Nx_lim = 1+int(Rc/dx);
  int Ny_lim = 1+int(Rc/dy);
  int Nz_lim = 1+int(Rc/dz);    
    
  // Add up the weights for each point.
  // Get Gauss-Lagrange integration points and weights for the 3D integrations
  gsl_integration_glfixed_table *tr = gsl_integration_glfixed_table_alloc(Ngauss_);
  double  gauss_p[Ngauss_];
  double  gauss_w[Ngauss_];
  for(int i=0;i<Ngauss_;i++) gsl_integration_glfixed_point(0, 1, i, gauss_p+i, gauss_w+i, tr);

  cout << "Ngauss_ = " << Ngauss_ << endl;
    
  double global_factor = dx*dx*dy*dy*dz*dz/(6*6*6);

  long Nmax = (((Nx_lim+1)*(Nx_lim+1+1)*(Nx_lim+1+2))/6) + (((Nx_lim+1)*(Nx_lim+1))/2) + Nx_lim+1;
  int chunk = 1000; 
  size_t steps_completed = 0;
  size_t total_steps = Nmax/chunk;

  vector<double> w2(Nmax,0.0);

  long p = 0;
#pragma omp parallel shared( chunk, w2, gauss_p, gauss_w ) private(p)
  {
  size_t local_count = 0;    
#pragma omp for  
  for(p=0;p<Nmax;p++)
    {
      double d = pow(3*p+sqrt(9*p*p-1.0/27.0),1.0/3.0);
      int ix = (p == 0 ? 0 : int(d+(1.0/(3*d))-1));
      int iy = (p == 0 ? 0 : int(0.5*sqrt(8*(p-ix*(ix+1)*(ix+2)/6)+1)-0.5));
      while(iy > ix)	      
	{ix++; iy = (p == 0 ? 0 : int(0.5*sqrt(8*(p-ix*(ix+1)*(ix+2)/6)+1)-0.5));}      
      int iz = (p == 0 ? 0 : p-((ix*(ix+1)*(ix+2))/6)-((iy*(iy+1))/2));
      while(iz > iy)
	{ iy++; iz = (p == 0 ? 0 : p-((ix*(ix+1)*(ix+2))/6)-((iy*(iy+1))/2));}
            
      double sum = 0.0;	      
      for(int Gx=0; Gx < Ngauss_; Gx++)
	{
	  double x  = gauss_p[Gx];
	  double wx = gauss_w[Gx];
	  for(int Gy=0; Gy < Ngauss_; Gy++)
	    {
	      double y  = gauss_p[Gy];
	      double wy = gauss_w[Gy];
	      for(int Gz=0; Gz < Ngauss_; Gz++)
		{
		  double z  = gauss_p[Gz];
		  double wz = gauss_w[Gz];
		  
		  for(int i=-1;i<=1;i++)
		    for(int j=-1;j<=1;j++)
		      for(int k=-1;k<=1;k++)
			for(int a=-1;a<=1;a+=2)
			  for(int b=-1;b<=1;b+=2)
			    for(int c=-1;c<=1;c+=2)
			      {
				double fx = wx*(1-x)*((1-x-x*x)*(i == 0 ? 2 : -1)+2-3*i*a*x);
				double fy = wy*(1-y)*((1-y-y*y)*(j == 0 ? 2 : -1)+2-3*j*b*y);
				double fz = wz*(1-z)*((1-z-z*z)*(k == 0 ? 2 : -1)+2-3*k*c*z);
				
				//double r = sqrt((ix+x*a+i)*(ix+x*a+i)*dx*dx+(iy+y*b+j)*(iy+y*b+j)*dy*dy+(iz+z*c+k)*(iz+z*c+k)*dz*dz);
				//double watt = v.Watt(r);

				double r2 = (ix+x*a+i)*(ix+x*a+i)*dx*dx+(iy+y*b+j)*(iy+y*b+j)*dy*dy+(iz+z*c+k)*(iz+z*c+k)*dz*dz;
				double watt = v.Watt2(r2);
				
				double w = global_factor*fx*fy*fz*watt;    
				
				sum += w;
			      }
		}
	    }
	}
      w2[p] = sum;	    
      
      if(local_count++ % chunk == chunk-1)
	{
#pragma omp atomic	    
	  ++steps_completed;
	  if (steps_completed % 100 == 1)
	    {
#pragma omp critical
	      cout << '\r';
	      std::cout << "\tProgress: " << steps_completed << " of " << total_steps << " (" << std::fixed << std::setprecision(1) << (100.0*steps_completed/total_steps) << "%)";
	      cout.flush();
	    }
	}
    }
  }
  w_att_.Real().zeros();
  a_vdw_ = 0.0;

  for(int ix = -Nx_lim;ix<=Nx_lim; ix++)
    for(int iy = -Ny_lim;iy<=Ny_lim; iy++)
      for(int iz = -Nz_lim;iz<=Nz_lim; iz++)
	{
	  // get the value of the integrated potential at this point
	  int nx = abs(ix);
	  int ny = abs(iy);
	  int nz = abs(iz);

	  if(ny > nx) swap(nx,ny);
	  if(nz > nx) swap(nx,nz);
	  if(nz > ny) swap(ny,nz);
	  long p = ((nx*(nx+1)*(nx+2))/6) + ((ny*(ny+1))/2) + nz;
	  double w = w2[p];

	  // and get the point it contributes to
	  int jx = ix;
	  int jy = iy;
	  int jz = iz;
	  s1_.getDensity().putIntoBox(jx,jy,jz);	  
	  long pos = jz+Nz*(jy+Ny*jx);

	  w_att_.Real().IncrementBy(pos,w/(dx*dy*dz));
	  a_vdw_ += w;
	}
  cout << " - Finished!" << endl;
  // In real usage, we calculate sum_R1 sum_R2 rho1 rho2 w(R1-R2). For a uniform system, this becomes
  // rho^2 N sum_R w(R). Divide by the volume and by rho^2 to get the vdw constant gives sum_R w(R)/(dx*dy*dz)
  a_vdw_ /= (dx*dy*dz); 

  cout << std::defaultfloat << std::setprecision(6);
    
  cout << myColor::GREEN;
  cout << "/////  Finished.  " << endl;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << myColor::RESET << endl;
    
  /// Dump the weights
  ofstream of(ss.str().c_str(), ios::binary | ios::app);  
  w_att_.Real().save(of);

}


void Interaction_Linear_Interpolation::generateWeights(Potential1 &v, stringstream &ss, Log& log)
{    
  const Density &density = s1_.getDensity();
      
  // The lattice
  long Nx = density.Nx();
  long Ny = density.Ny();
  long Nz = density.Nz();

  double dx = density.getDX();
  double dy = density.getDY();
  double dz = density.getDZ();
  
  // Set up the mean field potential
  // We need to take into account the whole contribution of the potential out to its cutoff of Rc.
  // This may mean going beyond nearest neighbors in certain conditions.
  // We also compute the vdw parameter at the same time.
      
  double Rc = v.getRcut();

  int Nx_lim = 1+int(Rc/dx);
  int Ny_lim = 1+int(Rc/dy);
  int Nz_lim = 1+int(Rc/dz);    
    
  // Add up the weights for each point.
  double global_factor = dx*dx*dy*dy*dz*dz/(120*120*120);

  long Nmax = (((Nx_lim+1)*(Nx_lim+1+1)*(Nx_lim+1+2))/6) + (((Nx_lim+1)*(Nx_lim+1))/2) + Nx_lim+1;
  int chunk = 1000; 
  size_t steps_completed = 0;
  size_t total_steps = Nmax/chunk;

  vector<double> w2(Nmax,0.0);
  int vv[5] = {1,26,66,26,1};
  
  long p = 0;
#pragma omp parallel shared( chunk, w2, vv) private(p)
  {
    size_t local_count = 0;
    long pold = -1;
    int ix,iy,iz;
#pragma omp for  
    for(p=0;p<Nmax;p++)
      {
	if(pold > 0 && (p-pold) < 100) // normally, the difference is just 1
	  {
	    for(int i=0;i<p-pold;i++)
	      {
		iz++;
		if(iz > iy) {iy++; iz=0;}
		if(iy > ix) {ix++; iy=0;}
	      }
	  } else {
	  double d = pow(3*p+sqrt(9*p*p-1.0/27.0),1.0/3.0);
	  ix = (p == 0 ? 0 : int(d+(1.0/(3*d))-1));
	  iy = (p == 0 ? 0 : int(0.5*sqrt(8*(p-ix*(ix+1)*(ix+2)/6)+1)-0.5));
	  while(iy > ix)	      
	    {ix++; iy = (p == 0 ? 0 : int(0.5*sqrt(8*(p-ix*(ix+1)*(ix+2)/6)+1)-0.5));}      
	  iz = (p == 0 ? 0 : p-((ix*(ix+1)*(ix+2))/6)-((iy*(iy+1))/2));
	  while(iz > iy)
	    { iy++; iz = (p == 0 ? 0 : p-((ix*(ix+1)*(ix+2))/6)-((iy*(iy+1))/2));}
	}
	pold = p;
	// check
	if(iy > ix || iz > iy || ix*(ix+1)*(ix+2)+3*iy*(iy+1)+6*iz != 6*p)
	  throw std::runtime_error("Bad indices in generateWeights");

      
	for(int i=0;i<5;i++)
	  for(int j=0;j<5;j++)
	    for(int k=0;k<5;k++)
	      {
		double r2 = (ix+i-2)*(ix+i-2)*dx*dx+(iy+j-2)*(iy+j-2)*dy*dy+(iz+k-2)*(iz+k-2)*dz*dz;
		w2[p] += global_factor*vv[i]*vv[j]*vv[k]*v.Watt2(r2);
	      }
	if(local_count++ % chunk == chunk-1)
	  {
#pragma omp atomic	    
	    ++steps_completed;
	    if (steps_completed % 100 == 1)
	      {
#pragma omp critical
		cout << '\r';
		std::cout << "\tProgress: " << steps_completed << " of " << total_steps << " (" << std::fixed << std::setprecision(1) << (100.0*steps_completed/total_steps) << "%)";
		cout.flush();
	      }
	  }
      }
  }
  w_att_.Real().zeros();
  a_vdw_ = 0.0;
  
  for(int ix = -Nx_lim-1;ix<=Nx_lim+1; ix++)
    {
      cout << '\r';
      std::cout << "\tProgress: " << ix+Nx_lim+1 << " of " << 2*(Nx_lim+1) << " (" << std::fixed << std::setprecision(1) << (100.0*(ix+Nx_lim+1)/(2*(Nx_lim+1))) << "%)";
      cout.flush();
      
      for(int iy = -Ny_lim-1;iy<=Ny_lim+1; iy++)
	for(int iz = -Nz_lim-1;iz<=Nz_lim+1; iz++)
	  {	  
	    // get the value of the integrated potential at this point
	    int nx = abs(ix);
	    int ny = abs(iy);
	    int nz = abs(iz);

	    if(ny > nx) swap(nx,ny);
	    if(nz > nx) swap(nx,nz);
	    if(nz > ny) swap(ny,nz);
	    long p = ((nx*(nx+1)*(nx+2))/6) + ((ny*(ny+1))/2) + nz;
	    double w = w2[p];
	    // and get the point it contributes to
	    long pos = s1_.getDensity().get_PBC_Pos(ix,iy,iz);
	    w_att_.Real().IncrementBy(pos,w/(dx*dy*dz));
	    a_vdw_ += w;
	  }
    }
  cout << " - Finished!" << endl;
  // In real usage, we calculate sum_R1 sum_R2 rho1 rho2 w(R1-R2). For a uniform system, this becomes
  // rho^2 N sum_R w(R). Divide by the volume and by rho^2 to get the vdw constant gives sum_R w(R)/(dx*dy*dz)
  a_vdw_ /= (dx*dy*dz); 

  cout << std::defaultfloat << std::setprecision(6);
    
  cout << myColor::GREEN;
  cout << "/////  Finished.  " << endl;
  cout << "///////////////////////////////////////////////////////////" << endl;
  cout << myColor::RESET << endl;
    
}
