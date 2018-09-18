
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>
#include <unistd.h>

#include <gsl/gsl_integration.h>

#include <unordered_map>

using namespace std;

extern char   __BUILD_DATE;
extern char   __BUILD_NUMBER;

#ifdef USE_OMP
#include <omp.h>
#endif

#include <mgl2/mgl.h>

#include "Integrator.h"
#include "Grace.h"
#include "options.h"
#include "TimeStamp.h"


#include "DFT.h"
#include "Droplet.h"

#include "Minimizer.h"

#include <gsl/gsl_sf_legendre.h>


double dx;

void translate(Droplet& theDensity)
{
  // do translation
  int Nx = theDensity.Nx(); 
  int Ny = theDensity.Ny(); 
  int Nz = theDensity.Nz();

  unsigned count = 0;

  double q6av = 0.0;
  
  vector< vector<double> > pos;

  for(int i = 0; i < Nx; ++i)
    for(int j = 0; j < Ny; ++j)
      for(int k = 0; k < Nz; ++k)
	{
	  int ip = (i+1 == Nx ? 0 : i+1);
	  int im = (i-1 == -1 ? Nx-1 : i-1);

	  int jp = (j+1 == Ny ? 0 : j+1);
	  int jm = (j-1 == -1 ? Ny-1 : j-1);

	  int kp = (k+1 == Nz ? 0 : k+1);
	  int km = (k-1 == -1 ? Nz-1 : k-1);


	  double d = theDensity.getDensity(i,j,k);

	  if(d < 3.5) continue;

	  bool bstop = false;
	  
	  for(int a = -1; a <= 1 && !bstop; a++)
	    {
	      double ia = i+a;
	      if(ia < 0)    ia = Nx-1;
	      if(ia > Nx-1) ia = 0;

	      for(int b = -1; b <= 1 && !bstop; b++)
		{
		  double jb = j+b;
		  if(jb < 0)    jb = Ny-1;
		  if(jb > Ny-1) jb = 0;

		  for(int c = -1; c <= 1 && !bstop; c++)
		    {
		      double kc = k+c;
		      if(kc < 0)    kc = Nz-1;
		      if(kc > Nz-1) kc = 0;		  

		      if(abs(a)+abs(b)+abs(c) < 1) continue;
		      if(abs(a)+abs(b)+abs(c) > 1) continue;
		      
		      double dd = theDensity.getDensity(ia,jb,kc);

		      if(dd > d) bstop = true;
		    }
		}
	    }
	  if(bstop) continue;

	  count++;
	  //	  cout << "Particle " << count << " found: " << endl;
	  vector<double> p;
	  p.push_back(i*dx);
	  p.push_back(j*dx);
	  p.push_back(k*dx);

	  pos.push_back(p);	  
	}

  double rcut = 2;
  
  for(long i = 0;i<pos.size();i++)
    {
      int nbonds = 0;
      double Q6[6]; for(int s=0;s<6;s++) Q6[s] = 0.0;
      
      for(long j = i+1;j<pos.size();j++)
	{
	  double dx[3];
	  dx[0] = pos[i][0] - pos[j][0];
	  dx[1] = pos[i][1] - pos[j][1];
	  dx[2] = pos[i][2] - pos[j][2];

	  double r = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);

	  if(r > rcut) continue;
	  
	  double x = dx[2]/r;  // this is cos(theta)
		  
	  for(int m=0;m<6;m++)
	    Q6[m] += (m == 0 ? 1 : 2)*gsl_sf_legendre_sphPlm(6,m,x);
		  
	  nbonds++;
		  
	}

      if(nbonds > 0)
	for(int m=0;m<6;m++) Q6[m] /= nbonds;
      
      double q6 = 0;
      
      for(int m=0;m<6;m++)
	q6 += Q6[m]*Q6[m];
      q6 = sqrt(q6*4*M_PI/(2*6+1));

      q6av += q6;
      
      //      cout << q6 << endl;
    } 

  cout << "num = " << pos.size() << " q6 av = " << (pos.size() == 0 ? 0 : q6av/pos.size()) << endl;

	  
	  


}

int main(int argc, char** argv)
{
  cout << "Build date  : " <<  (unsigned long) &__BUILD_DATE << endl;
  cout << "Build number: " << (unsigned long) & __BUILD_NUMBER << endl;


  double L[3] = {10,10,10};
  int PointsPerHardSphere = 5;
  int nCores = 1;

  double R = -1;
  double zPos = 0;

  double density = 0.8;

  double alpha = 0.01;
  int MaxIts = 100;

  double epsWall = 1;
  double sigWall = 1;

  double kT = 1;

  double eps = 1;
  double sigma = 1;
  double rcut = 1;

  double sigma_conj_grad = 1;

  string pointsFile("..//SS31-Mar-2016//ss109.05998");
  string outfile("dump.dat");
  string infile;
  double Natoms = -1;

  double Mixing = 0;

  Options options;

  options.addOption("nCores", &nCores);
  options.addOption("BulkDensity", &density);
  options.addOption("PointsPerHardSphere", &PointsPerHardSphere);

  options.addOption("kT", &kT);

  options.addOption("eps", &eps);
  options.addOption("sigma", &sigma);
  options.addOption("rcut", &rcut);

  options.addOption("Lx", L);
  options.addOption("Ly", L+1);
  options.addOption("Lz", L+2);

  options.addOption("R", &R);
  options.addOption("zPosition", &zPos);

  options.addOption("alpha", &alpha);
  options.addOption("MaxIts", &MaxIts);

  options.addOption("sigWall", &sigWall);
  options.addOption("epsWall", &epsWall);

  options.addOption("sigma_conj_grad", &sigma_conj_grad);

  options.addOption("OutputFile", &outfile);
  options.addOption("IntegrationPointsFile", &pointsFile);
  options.addOption("InputFile", &infile);
  options.addOption("Natoms", &Natoms);
  options.addOption("MixingParameter", &Mixing);
  options.read(argc, argv);

  dx = 1.0/PointsPerHardSphere;

  Droplet theDensity(dx, L, PointsPerHardSphere, infile, 0, 0);

  if(! infile.empty())
    theDensity.readDensity(infile.c_str());
  else throw std::runtime_error("Could not read input density");

  translate(theDensity);


  return 1;
}
