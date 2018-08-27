
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
#include "Droplet_Vis.h"
//#include "LJ.h"

#include "Minimizer.h"
//#include "Picard.h"

#include "visit_writer.h"

static   mglData a;
static mglGraph gr_;


void translate(Droplet_VIS& theDensity)
{
  // do translation
  int NX = theDensity.Nx(); //*0.75; //2;
  int NY = theDensity.Ny(); //*0.75; ///2;
  int NZ = theDensity.Nz();

  // I don't understand why it has to be plus 1 ...
  int dims[] = {NX+1, NY+1, NZ+1};

  int nvars = 1;
  int vardims[] = {1};
  int centering[] = {0};
  const char *varnames[] = {"density"};

  unsigned Nmax = NX;
  Nmax *= NY;
  Nmax *= NZ;
  
  float *density = new float[Nmax];
  float *vars[] = {(float *)density};

  /* Create zonal variable */
  unsigned pos = 0;
  for(int k = 0; k < NZ; ++k)
    for(int j = 0; j < NY; ++j)
      for(int i = 0; i < NX; ++i)
	{
	  density[pos] = (theDensity.getDensity(i,j,k));
	  pos++;

	  //	  density[pos] = theDensity.getDensity(i+NX/6,j+NY/6,k);	  
	  //	  density[k][j][i] = theDensity.getDensity(i+NX/6,j+NY/6,k);
	  //  	density[k][j][i] = theDensity.getDensity(i+NX/2,j+NY/2,k);
	}
  cout << "pos 0 = " << density[0] << endl;
  /* Use visit_writer to write a regular mesh with data. */
  write_regular_mesh("density.vtk", 0, dims, nvars, vardims,
		     centering, varnames, vars);
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

  double dx = 1.0/PointsPerHardSphere;

  //  L[2] += int(rcut) + 1; // expand volume to make sure there are no interactions;

  Droplet_VIS theDensity(dx, L, PointsPerHardSphere, R, zPos, sigWall, epsWall/kT); 

  // Determine coexistence to give us a baseline for the densities
  double xliq_coex = 1.0;
  double xgas_coex = 1.0;

  theDensity.initialize(xliq_coex,xgas_coex);

  if(! infile.empty())
    theDensity.readDensity(infile.c_str());
  else throw std::runtime_error("Could not read input density");

  translate(theDensity);


  /*

  double cbmin = 0.0;
  double cbmax = 1.0;
  cout << "Enter cb min: "; cin >> cbmin;
  cout << "Enter cb max: "; cin >> cbmax;

  gr_.RunThr();

  theDensity.initialize_2D_data(a);

  for(double fac = 0.4; fac < 1.0; fac += 0.001)
    {
    theDensity.setCut(fac);

    theDensity.fill_2D_data(a);


    //clear the window
    gr_.Clf();	


    // basic window formating
    gr_.Box();  
    gr_.Alpha(false);
    gr_.SetRange('c', cbmin, cbmax);

    // set the color scheme
    //	  gr_.Dens(a,"kRryw");
    //	  gr_.Dens(a,"UbcyqR");
    gr_.Dens(a,"kBbcw");

  // Write a title
    char str[48];	
    snprintf(str,48,"fac  = %lf",fac);
    gr_.Puts(mglPoint(0,1.1),str);

    // update the window
    gr_.Update();	

    //    cout << "Enter cb min: "; cin >> cbmin;
    //    cout << "Enter cb max: "; cin >> cbmax;

    //    usleep(1000000);

    } ;

  */

  return 1;
}
