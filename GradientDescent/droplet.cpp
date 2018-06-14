
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

using namespace std;

extern char   __BUILD_DATE;
extern char   __BUILD_NUMBER;

#ifdef USE_OMP
#include <omp.h>
#endif


#include "Grace.h"
#include "options.h"
#include "TimeStamp.h"


#include "Droplet.h"

#include "StringMethod.h"

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

  double kT = 1;

  double eps = 1;
  double sigma = 1;
  double rcut = 1;

  double sigma_conj_grad = 1;

  string pointsFile("..//SS31-Mar-2016//ss109.05998");
  string outfile("dump.dat");
  string infile;
  double Natoms = -1;

  bool showGraphics = true;
  double forceLimit = -1;

  int Nimages = 10;
  bool bRestart = false;
  int Jclimber = -1;
  bool freeEnd = false;
  
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

  options.addOption("sigma_conj_grad", &sigma_conj_grad);
  options.addOption("ForceTerminationCriterion",&forceLimit);

  options.addOption("OutputFile", &outfile);
  options.addOption("IntegrationPointsFile", &pointsFile);
  options.addOption("InputFile", &infile);
  options.addOption("Natoms", &Natoms);
  options.addOption("ShowGraphics", &showGraphics);

  options.addOption("Nimages",&Nimages);
  options.addOption("Restart",&bRestart);
  options.addOption("Jclimber",&Jclimber);
    options.addOption("FreeEnd",&freeEnd);
  
  options.read(argc, argv);

  ofstream log("log.dat");
  TimeStamp ts;
  log << "# " << ts << endl;
  log << "#=================================" << endl << "#" << endl;
  log << "#Input parameters:" << endl <<  "#" << endl;
  options.write(log);
  log << "#=================================" << endl;
  log.close();

  double dx = 1.0/PointsPerHardSphere;

#ifdef USE_OMP    
  omp_set_dynamic(0);
  omp_set_num_threads(nCores);

  int fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif

  Droplet finalDensity(dx, L, PointsPerHardSphere, R, zPos); 
  Potential1 potential(sigma, eps, rcut);
  DFT_VDW<RSLT> dft(finalDensity,potential,pointsFile,kT);

  // Determine coexistence to give us a baseline for the densities
  double xliq_coex = 0.001;
  double xgas_coex = 0.6;
  try{
    dft.coexistence(xliq_coex, xgas_coex);
  } catch(...) { cout << "No coexistence found ... " << endl;}
  
  double xs1,xs2;
  dft.spinodal(xs1,xs2);
  cout << "Spinodal: " << xs1 << " " << xs2 <<  " so Nspinodal = " << xs1*finalDensity.getVolume() << endl;
  cout << "Coexistence: " << xgas_coex << " " << xliq_coex << endl;

  // set the thermodynamic state
  double omega_coex = dft.Omega(xliq_coex);
  double mu         = dft.Mu(xliq_coex);
  finalDensity.initialize(xliq_coex,xgas_coex);

  //  double NN = 297.0; //finalDensity.getNtotal();
  double NN = finalDensity.getNumberAtoms();
  cout << "NN = " << NN << endl;

  if(! infile.empty())
    finalDensity.readDensity(infile.c_str());

  if(Natoms > 0) NN = Natoms;
  else NN = finalDensity.getNumberAtoms();

  cout << "Final NN = " << finalDensity.getNumberAtoms() << endl;

  cout << "Hard sphere diameter  = " << dft.HSD() << endl;
  cout << "Coexisting densities  = " << xliq_coex << " " << xgas_coex << endl;  
  cout << "Chemical potential/kT = " << mu << endl;
  cout << "beta * Grand free energy per unit volume = " << omega_coex << endl;
  string s("log.dat");


  double avDensity = NN/finalDensity.getVolume();

  cout << "Average density = " << avDensity << endl;


  ofstream log1("log.dat");
  
  dft.setEtaMax(1.0-1e-8);


  Grace grace(800,600,1);
  bool bFixedBoundaries = true;


  int ret = system("rm string_graph.agr");


  // lower all densities on the boundaries a small amount
  int Nx = finalDensity.Nx();
  int Ny = finalDensity.Ny();
  int Nz = finalDensity.Nz();

  for(int ix=0;ix<Nx;ix++)
    for(int iy=0;iy<Ny;iy++)
      {
	finalDensity.set_Density_Elem(ix,iy,0,finalDensity.getDensity(ix,iy,0)*0.90);
	finalDensity.set_Density_Elem(ix,iy,Nz-1,finalDensity.getDensity(ix,iy,Nz-1)*0.90);
      }
  for(int ix=0;ix<Nx;ix++)
    for(int iz=1;iz<Nz-1;iz++)
      {
	finalDensity.set_Density_Elem(ix,0,iz,finalDensity.getDensity(ix,0,iz)*0.90);
	finalDensity.set_Density_Elem(ix,Ny-1,iz,finalDensity.getDensity(ix,Ny-1,iz)*0.90);
      }

  for(int iy=1;iy<Ny-1;iy++)
    for(int iz=1;iz<Nz-1;iz++)
      {
	finalDensity.set_Density_Elem(0,iy,iz,finalDensity.getDensity(0,iy,iz)*0.90);
	finalDensity.set_Density_Elem(Nx-1,iy,iz,finalDensity.getDensity(Nx-1,iy,iz)*0.90);
      }

  

  
  DDFT ddft(dft,finalDensity,bFixedBoundaries,&grace,showGraphics);
  ddft.initialize();

  string slog("log.dat");

  //  ddft.test_solv_tridiag();
  
  ddft.run(slog);

  return 1;
}
