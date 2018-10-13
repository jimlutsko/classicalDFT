
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>
#include <unistd.h>
#include <algorithm>

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



void translate(Droplet& theDensity, double dx, double rneighbor, double minPeak, double &nparticles, double &q6, double &q6median)
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

	  if(d < minPeak) continue;

	  bool bstop = false;
	  
	  for(int a = -1; a <= 1.5 && !bstop; a++)
	    {
	      double ia = i+a;
	      if(ia < 0)    ia = Nx-1;
	      if(ia > Nx-1) ia = 0;

	      for(int b = -1; b <= 1.5 && !bstop; b++)
		{
		  double jb = j+b;
		  if(jb < 0)    jb = Ny-1;
		  if(jb > Ny-1) jb = 0;

		  for(int c = -1; c <= 1.5 && !bstop; c++)
		    {
		      double kc = k+c;
		      if(kc < 0)    kc = Nz-1;
		      if(kc > Nz-1) kc = 0;		  

		      if(abs(a)+abs(b)+abs(c) < 1) continue;
		      //		      if(abs(a)+abs(b)+abs(c) > 1) continue;
		      
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


  vector<double> AllQ6;
  
  for(long i = 0;i<pos.size();i++)
    {
      int nbonds = 0;
      double Q6[6]; for(int s=0;s<6;s++) Q6[s] = 0.0;

      //      cout << i << " " << rneighbor << endl;
      
      for(long j = i+1;j<pos.size();j++)
	{
	  double dx[3];
	  dx[0] = pos[i][0] - pos[j][0];
	  dx[1] = pos[i][1] - pos[j][1];
	  dx[2] = pos[i][2] - pos[j][2];

	  double r = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);

	  if(r > rneighbor) continue;
	  //	  cout << "\t " << r << endl;	  
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
      AllQ6.push_back(q6);
    } 

  nparticles = pos.size();

  q6 = (pos.size() == 0 ? 0 : q6av/pos.size());

  if(AllQ6.size() > 0)
    {
      sort(AllQ6.begin(), AllQ6.end());  
      q6median = AllQ6[AllQ6.size()/2];
    } else q6median = 0;

  //  if(AllQ6.size() > 0)
  //    cout << "\t" << pos.size() << " " << AllQ6.size() << " " << AllQ6[0] << " " << AllQ6.back() << " " << AllQ6[AllQ6.size()/2] << endl;
  //  else cout << "0 0 0" << endl;
}

int main(int argc, char** argv)
{
  cout << "Build date  : " <<  (unsigned long) &__BUILD_DATE << endl;
  cout << "Build number: " << (unsigned long) & __BUILD_NUMBER << endl;


  double L[3] = {10,10,10};
  int PointsPerHardSphere = 5;

  double rneighbor = 2;
  double threshold = 5;
  
  string infile;
  double Natoms = -1;

  Options options;

  options.addOption("PointsPerHardSphere", &PointsPerHardSphere);

  options.addOption("Lx", L);
  options.addOption("Ly", L+1);
  options.addOption("Lz", L+2);

  options.addOption("InputFile", &infile);

  options.addOption("NeighborThreshold", &rneighbor);
  options.addOption("MinPeak", &threshold);
  
  options.read(argc, argv);

  ofstream outfile("nout.dat");
  options.write(outfile);
  outfile << "#i\tnpos\tq6\tq6median" << endl;

  double dx = 1.0/PointsPerHardSphere;


  for(int i=0;i<1000;i++)
    {
      stringstream infile_name;
      infile_name << "archive_" << std::internal << std::setfill('0') << std::setw(4) << i << ".dat";

      string infile = infile_name.str();
      
      try{
	Droplet theDensity(dx, L, PointsPerHardSphere, infile, 0.0, 0.0);
	theDensity.readDensity(infile.c_str());
	double npos = 0;
	double q6   = 0;
	double q6median = 0;
	
	translate(theDensity, dx, rneighbor, threshold, npos, q6, q6median);       
	cout << i << " " << npos << " " << q6 << " " << q6median << endl;
	outfile << i << "\t" << npos << "\t" << q6 << "\t" << q6median << endl;
      } catch (...) {
	cout << "Coult not open file ... " << infile << endl;
	break;
      }
    }

  return 1;
}
