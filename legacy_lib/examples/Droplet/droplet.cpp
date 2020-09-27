#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

using namespace std;

#ifdef USE_OMP
#include <omp.h>
#endif

#include <gsl/gsl_sf_lambert.h>

#include "options.h"
#include "TimeStamp.h"

#include "DFT.h"
#include "Droplet.h"
#include "Minimizer.h"
#include "Log.h"
#include "myColor.h"


int main(int argc, char** argv)
{
  double L[3] = {10,10,10};
  int PointsPerHardSphere = 5;
  int nCores = 6;

  double R = -1;
  double zPos = 0;

  double density_inside_1 = 0.01;
  double density_outside_1 = 0.8;
  double eta2 = 0.2;

  double hsd1 = -1;
  double hsd2 = -1;
  
  double epsWall = 1;
  double sigWall = 1;

  double kT = 1;

  double eps1   = 1;
  double sigma1 = 1;
  double rcut1  = 3;

  double eps2   = 1;
  double sigma2 = 1;
  double rcut2  = 3;

  double eps_surf   = 1;
  double sigma_surf = 1;
  double rcut_surf  = 3;

  string pointsFile("..//SS31-Mar-2016//ss109.05998");
  string outfile("dump.dat");
  string infile;
  
  double forceLimit = 1e-4;
  double dt = 1e-3;
  double dtMax = 1;
  double alpha_start = 0.01;

  bool showGraphics = true;
  double Nfixed = -1;
  
  double Asurf = 0;
  double rho_surf = -1;
  
  Options options;

  options.addOption("nCores", &nCores);
  options.addOption("Interior_Density", &density_inside_1);
  options.addOption("Exterior_Density", &density_outside_1);

  options.addOption("Nfixed", &Nfixed);  

  options.addOption("PointsPerHardSphere", &PointsPerHardSphere);

  options.addOption("kT", &kT);

  options.addOption("eps1",   &eps1);
  options.addOption("sigma1", &sigma1);
  options.addOption("rcut1",  &rcut1);

  options.addOption("Lx", L);
  options.addOption("Ly", L+1);
  options.addOption("Lz", L+2);

  options.addOption("OutputFile", &outfile);
  options.addOption("IntegrationPointsFile", &pointsFile);

  options.addOption("ForceTerminationCriterion",&forceLimit);
  options.addOption("TimeStep", &dt);
  options.addOption("TimeStepMax", &dtMax);
  options.addOption("AlphaStart", &alpha_start);

  options.addOption("InputFile", &infile);
  
  options.addOption("ShowGraphics", &showGraphics);
  
  options.read(argc, argv);

  Log log("log.dat");
  log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;

  log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;

  options.write(log);
  log <<  myColor::GREEN << "=================================" <<  myColor::RESET << endl;

  double dx = 1.0/PointsPerHardSphere;

  if(L[0] < 0) L[0] = dx;
  if(L[1] < 0) L[1] = dx;
  if(L[2] < 0) L[2] = dx;  
  
#ifdef USE_OMP    
  omp_set_dynamic(0);
  omp_set_num_threads(nCores);

  int fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());
  log << "omp max threads = " << omp_get_max_threads() << endl;
#endif

  
  //////////////////////////////////////
  ////// Create potential && effective hsd

  LJ potential1(sigma1, eps1, rcut1);

  if(hsd1 < 0) hsd1 = potential1.getHSD(kT);
  
  /////////////////////////////////////
  // Create density objects
  Droplet theDensity1(dx, L, hsd1);

  /////////////////////////////////////
  // Create the species objects
  
  FMT_Species species1(theDensity1,hsd1,pointsFile);

  Interaction i1(species1,species1,potential1,kT,log,pointsFile);

  /////////////////////////////////////
  // Create the hard-sphere object
  RSLT fmt;
  
  /////////////////////////////////////
  // DFT object

  DFT dft(&species1);

  dft.addHardCoreContribution(&fmt);  
  dft.addInteraction(&i1);
  
  /////////////////////////////////////////////////////
  // Report
  log << "Density 1 = " << density_inside_1 << endl;
  log << "HSD 1 = " << hsd1 << endl;
  
  /////////////////////////////////////////////////////
  // The thermodynamics

  vector<double> densities;
  densities.push_back(density_inside_1);

  double mu1 = dft.Mu(densities,0);
  
  log << "Omega/(V kT) = " << dft.Omega(densities) << endl;  
  log << "mu1 = " << mu1 << endl;

  theDensity1.initialize(density_inside_1, density_outside_1, 5);
  if(Nfixed > 0) theDensity1.scaleTo(Nfixed);

  ///////////////////////////////////////////////
  // Fix the mass of the surfactant species.

  double N = theDensity1.getNumberAtoms();
  species1.setFixedMass(N);
  log << "N fixed at " << N << endl;

  species1.setChemPotential(0.0); //mu1);
  
  if(! infile.empty())
    theDensity1.readDensity(infile.c_str());

  log <<  myColor::GREEN << "=================================" <<  myColor::RESET << endl;
  
  //fireMinimizer_Mu minimizer(dft,log);
  fireMinimizer2 minimizer(dft,log);
  minimizer.setForceTerminationCriterion(forceLimit);
  minimizer.setTimeStep(dt);
  minimizer.setTimeStepMax(dtMax);
  minimizer.setAlphaStart(alpha_start);
  minimizer.setAlphaFac(1.0);
  minimizer.run();

  double Natoms = theDensity1.getNumberAtoms();
  double Omega = minimizer.getF();
  double dOmega = Omega-L[0]*L[1]*L[2]*dft.Omega(densities);
  double SurfaceTension = dOmega/(L[0]*L[1]);

  log << "dft.Omega = " << dft.Omega(densities) << endl;  
  log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
  log << "Final Omega: " << Omega << endl;
  log << "Excess Omega = " << dOmega << endl;
  log << "Surface Tension = " << SurfaceTension << endl;
    
  return 1;
}
