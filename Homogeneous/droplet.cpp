
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


#include "DFT.h"
#include "Droplet.h"
#include "Minimizer.h"

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
  bool showGraphics = true;
  double forceLimit = -1;
  double scale = 1;

  double dt = 1e-3;
  double dtMax = dt*10;
  double alpha_start = 0.01;
  
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

    options.addOption("Scale", &scale);

  options.addOption("alpha", &alpha);
  options.addOption("MaxIts", &MaxIts);

  options.addOption("sigWall", &sigWall);
  options.addOption("epsWall", &epsWall);

  options.addOption("sigma_conj_grad", &sigma_conj_grad);
  options.addOption("ForceTerminationCriterion",&forceLimit);


  options.addOption("OutputFile", &outfile);
  options.addOption("IntegrationPointsFile", &pointsFile);
  options.addOption("InputFile", &infile);
  options.addOption("Natoms", &Natoms);
  options.addOption("MixingParameter", &Mixing);
  options.addOption("ShowGraphics", &showGraphics);

  options.addOption("TimeStep", &dt);
  options.addOption("TimeStepMax", &dtMax);

  options.addOption("AlphaStart", &alpha_start);

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

  ////// Create potential

  //  L[2] += int(rcut) + 1; // expand volume to make sure there are no interactions;

  Droplet theDensity(dx, L, PointsPerHardSphere, R, zPos, sigWall, epsWall/kT); 
  Potential1 potential(sigma, eps, rcut);
  DFT_VDW<RSLT> dft(theDensity,potential,pointsFile,kT);

  // Determine coexistence to give us a baseline for the densities
  double xliq_coex = 0.8;
  double xgas_coex = 1e-4;
  dft.coexistence(xgas_coex, xliq_coex);

  // set the thermodynamic state
  double omega_coex = dft.Omega(xliq_coex);
  double mu         = dft.Mu(xliq_coex);
  theDensity.initialize(xliq_coex,xgas_coex*scale);

  //  double NN = 297.0; //theDensity.getNtotal();
  double NN = theDensity.getNumberAtoms();
  cout << "NN = " << NN << endl;

  if(! infile.empty())
    theDensity.readDensity(infile.c_str());

  if(Natoms > 0) NN = Natoms;
  else NN = theDensity.getNumberAtoms();

  double xs1,xs2;
  dft.spinodal(xs1,xs2);

  cout << "Hard sphere diameter  = " << dft.HSD() << endl;
  cout << "Coexisting densities  = " << xliq_coex << " " << xgas_coex << endl;  
  cout << "Chemical potential/kT = " << mu << endl;
  cout << "beta * Grand free energy per unit volume = " << omega_coex << endl;
  cout << "Spinodals: " << xs1 << " " << xs2 << endl;
  cout << "Av density = " << NN/(L[0]*L[1]*L[2]) << endl;
  cout << "Homogeneous beta*F = " << L[0]*L[1]*L[2]*dft.Fhelmholtz(NN/(L[0]*L[1]*L[2])) << endl;
  cout << "#Nx Ny Nz = " << theDensity.Nx() << " " << theDensity.Ny() << " " << theDensity.Nz() << endl;

  ofstream log1("log.dat",ios::app);
  log1 << "#=================================" << endl << "#" << endl;
  log1 << "#NN = " << NN << endl;
  log1 << "#Hard sphere diameter  = " << dft.HSD() << endl;
  log1 << "#Coexisting densities  = " << xliq_coex << " " << xgas_coex << endl;  
  log1 << "#Chemical potential/kT = " << mu << endl;
  log1 << "#beta * Grand free energy per unit volume = " << omega_coex << endl;
  //  log1 << "#Spinodals: " << xs1 << " " << xs2 << endl;
  log1 << "#Av density = " << NN/(L[0]*L[1]*L[2]) << endl;
  log1 << "#Homogeneous beta*F = " << L[0]*L[1]*L[2]*dft.Fhelmholtz(NN/(L[0]*L[1]*L[2])) << endl;
  log1 << "#Nx Ny Nz = " << theDensity.Nx() << " " << theDensity.Ny() << " " << theDensity.Nz() << endl;
  log1 << "#=================================" << endl;
  log1.close();


  //  dft.setEtaMax(1.0-1e-8);
  
  string s("log.dat");

  fireMinimizer minimizer(dft, theDensity, NN);
  minimizer.setForceTerminationCriterion(forceLimit);
  minimizer.setTimeStep(dt);
  minimizer.setTimeStepMax(dtMax);
  minimizer.setAlphaStart(alpha_start);
  minimizer.run(s);

  /*
  
  ConjugateGradients2_Fixed cj(dft, theDensity, NN,showGraphics);
  cj.set_sigma_conj_grad_secant(sigma_conj_grad);
  cj.set_mixing_parameter(Mixing);
  cj.setForceTerminationCriterion(forceLimit);
  
  cj.run(s);
    
  */
  
  string dump(outfile.c_str());
  theDensity.writeDensity(dump);

  return 1;
}
