/* This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
 * To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
 *
 * Author: James F. Lutsko
 * www.lutsko.com
 */

#ifndef __LUTSKO__INTERACTION__
#define __LUTSKO__INTERACTION__

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <complex>
#include <complex.h>

#include "Species.h"
#include "Log.h"
#include "myColor.h"

/**
  *  @brief This class encapsulates the interaction between two species (or one species with itself)
  */  


class Interaction
{
 public:
 Interaction(Species &s1, Species &s2, Potential1 &v, double kT, Log &log) : s1_(s1), s2_(s2)
  {
    log << "Calculating mean field potential ... " << endl;
    const Density &density = s1.getDensity();
      
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

    w_att_.initialize(Nx,Ny,Nz);      
    a_vdw_ = 0.0;

    long ntot = long(2*Rc/dx) * long(2*Rc/dy) * long(2*Rc/dz);
    long count = 0;
    
    for(double x = -Rc; x <= Rc; x += dx)
      {
	long nx = x/dx;
	if(nx >= Nx || nx < 0) nx -= Nx*int(nx/Nx);
	while(nx >=  Nx) nx -= Nx;
	while(nx <   0) nx += Nx;

	for(double y = -Rc; y <= Rc; y += dy)
	  {
	    long ny = y/dy; 
	    if(ny >= Ny || ny < 0) ny -= Ny*int(ny/Ny);
	    while(ny >=  Ny) ny -= Ny;
	    while(ny <   0) ny += Ny;
	    
	    for(double z = -Rc; z <= Rc; z += dz)	    
	      {
		long nz = z/dz; 
		if(nz >= Nz || nz < 0) nz -= Nz*int(nz/Nz);

		while(nz >=  Nz) nz -= Nz;
		while(nz <   0) nz += Nz;

		long pos = nz+Nz*(ny+Ny*nx);

		double r2 = x*x+y*y+z*z;
		double w = v.Watt(sqrt(r2))/kT;
		//		double w = -exp(-v.V(sqrt(r2))/kT);
		//		if(r2 > 1.03*1.03) w += 1; 

		w *= dx*dy*dz;
		
		a_vdw_ += w;
		w_att_.Real().IncrementBy(pos,w);

		if(count % 10000 == 0)
		  {if(count > 0) cout << '\r'; cout << "\t" << int(double(count)*100.0/ntot) << "% finished: " << count << " out of " << ntot; cout.flush();}
		count++;
		
	      }	   
	  }
      }
    log << endl;

    // Now save the FFT of the field  
    w_att_.do_real_2_fourier();     
  }

 Interaction(Species &s1, Species &s2, Potential1 &v, double kT, Log &log, string &pointsFile) : s1_(s1), s2_(s2)
  {
    const Density &density = s1.getDensity();
    long Nx = density.Nx();
    long Ny = density.Ny();
    long Nz = density.Nz();

    w_att_.initialize(Nx,Ny,Nz);      
    a_vdw_ = 0.0;
  
    // First, try to read the weights from a file
    
    stringstream ss1;
    ss1 << "weights_" << s1_.getSequenceNumber() << "_" << s2_.getSequenceNumber() << ".dat";
  
    bool readWeights = true;
    
    ifstream in(ss1.str().c_str(), ios::binary);
    if(!in.good())
      {readWeights = false;     cout << myColor::GREEN  << "\n" <<  "Could not open file with potential kernal: it will be generated" << myColor::RESET << endl;}
    else {
      string buf;
      getline(in,buf);

      stringstream ss2(buf);
      int nx, ny, nz;
      double Rc;
      ss2 >> nx >> ny >> nz >> Rc;

      if(nx != Nx)
	{readWeights = false; cout << "\n" <<  "Mismatch in Nx: expected " << Nx << " but read " << nx <<  endl;}
      if(ny != Ny)
	{readWeights = false; cout << "\n" <<  "Mismatch in Ny: expected " << Ny << " but read " << ny <<  endl;}
      if(nz != Nz)
	{readWeights = false; cout << "\n" <<  "Mismatch in Nz: expected " << Nz << " but read " << nz <<  endl;}
      if(fabs(v.getRcut()-Rc) > 1e-8*(v.getRcut()+Rc))
	{readWeights = false; cout << "\n" <<  "Mismatch in cutoff: generating weights: expected " << v.getRcut() << " but read " << Rc << endl;}

      getline(in,buf);
      stringstream ss3(buf);
      if(ss3.str().compare(pointsFile) != 0)
	{readWeights = false; cout << "\n" <<  "Mismatch in points file: expected " << pointsFile << " but read " << ss3.str() <<  endl;}
    }

    if(readWeights)
      {
	w_att_.Real().load(in);
	a_vdw_ = w_att_.Real().accu();	
      } else {
      generateWeights(pointsFile,v, ss1, log, kT);
    }
    // Now generate the FFT of the field  
    w_att_.do_real_2_fourier();     
  }    



  void generateWeights(string &pointsFile, Potential1 &v, stringstream &ss, Log& log, double kT)
  {    
    log << "Calculating mean field potential ... " << endl;
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

    cout << myColor::RESET << endl;


    // Add up the weights for each point.
    // We throw in all permutations (iperm loop) of the integration points and all reflections (is loop)
    // to establish as much symmetry as possible. Probably unnecessary.    
    cout << endl;
    cout << myColor::GREEN;
    cout << "///////////////////////////////////////////////////////////" << endl;
    cout << "/////  Generating integration points for sphere volume " << endl;
    cout << myColor::RESET << endl;

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
		  
		    double watt = v.Watt(r)/kT; 
		  
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
    ofstream of(ss.str().c_str(), ios::binary);

    of.flags (std::ios::scientific);
    of.precision (std::numeric_limits<double>::digits10 + 1);

    of << Nx << " " << Ny << " " << Nz << " " << Rc << endl;
    of << pointsFile << endl;

    cout << "Writing: " << Ny << " " << Nz << " " << Rc << endl;
    cout << "Writing: " << pointsFile << endl;

    w_att_.Real().save(of);
  }
  
  
  // Note that the matrix w already contains a factor of dV
  double getInteractionEnergyAndForces()
  {
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

  double checkCalc(int jx, int jy, int jz)
  {
    const Density &density1 = s1_.getDensity();
    
    long Ntot = density1.Ntot();
    double dV = density1.dV();
    
    int Nx = density1.Nx();
    int Ny = density1.Ny();
    int Nz = density1.Nz();


    
    //    int jx = 0;
    //    int jy = 66;
    //    int jz = 14;

    /*
      double EE = 0.0;
    
      for(int jx=0; jx<Nx; jx++)
      for(int jy=0; jy<Ny; jy++)
      for(int jz=0; jz<Nz; jz++)
      {
    */
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
    //	    EE += dd*density1.getDensity(sj);
    //	      }
    //	    cout << "Direct calc of dF[" << jx << "," << jy << "," << jz << "] = " << dd*dV*dV  << endl;
    //    EE *= 0.5*dV*dV;
    //    cout << setprecision(20) << "E        = " << E << endl;
    //    cout << setprecision(20) << "E-direct = " << EE << endl;
    //    cout << setprecision(20) << "diff     = " << E - EE << endl;
    
    return dd*dV*dV;
  }

  double getW(long s)  { return w_att_.Real().get(s);}

  double Mu(const vector<double> &x, int species) const
  {
    double mu = 0.0;

    if(s1_.getSequenceNumber() == species)
      mu += a_vdw_*x[s2_.getSequenceNumber()];

    if(s2_.getSequenceNumber() == species)
      mu += a_vdw_*x[s1_.getSequenceNumber()];    

    return mu;
  }

  double Fhelmholtz(const vector<double> &x) const
  {
    return a_vdw_*x[s1_.getSequenceNumber()]*x[s2_.getSequenceNumber()];
  }  

  double getVDWParameter() const { return a_vdw_;}
  
 private:
  Species &s1_;
  Species &s2_;

  double a_vdw_;
  DFT_FFT w_att_;  

};

#endif // __LUTSKO__INTERACTION__
