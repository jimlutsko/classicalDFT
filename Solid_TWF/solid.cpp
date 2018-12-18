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

#include "options.h"
#include "TimeStamp.h"
#include "Integrator.h"

#include "DFT.h"
#include "Solid.h"
#include "Minimizer.h"

int main(int argc, char** argv)
{
  cout << "Build date  : " <<  (unsigned long) &__BUILD_DATE << endl;
  cout << "Build number: " << (unsigned long) & __BUILD_NUMBER << endl;

  int PointsPerSigma = 5;
  int nCores = 6;

  double density = 0.8;
  double kT = 1;

  double eps   = 1;
  double sigma = 1;
  double rcut  = 3;

  string pointsFile("..//SS31-Mar-2016//ss109.05998");
  string infile;

  double forceLimit = 1e-4;
  double dt = 1e-3;
  double dtMax = 1;
  double alpha_start = 0.01;

  int ncopy = 3;
  double Nvac = 1e-5;
  int Npoints = 100;

  double InitializationScaleFac = 0.7;
  double InitializationInverseWidth = 50;
  
  Options options;

  options.addOption("nCores", &nCores);
  options.addOption("PointsPerSigma", &PointsPerSigma);

  options.addOption("kT", &kT);

  options.addOption("Npoints", &Npoints);
  options.addOption("Ncopy", &ncopy);
  options.addOption("Nvac", &Nvac);

  options.addOption("InitializationScaleFac", &InitializationScaleFac);
  options.addOption("InitializationInverseWidth", &InitializationInverseWidth);
  
  options.addOption("eps", &eps);
  options.addOption("sigma", &sigma);
  options.addOption("rcut", &rcut);

  options.addOption("IntegrationPointsFile", &pointsFile);

  options.addOption("ForceTerminationCriterion",&forceLimit);
  options.addOption("TimeStep", &dt);
  options.addOption("TimeStepMax", &dtMax);
  options.addOption("AlphaStart", &alpha_start);

  options.addOption("InputFile", &infile);
  
  options.read(argc, argv);

  ofstream log1("log.dat");
  TimeStamp ts;
  log1 << "# " << ts << endl;
  log1 << "#=================================" << endl << "#" << endl;
  log1 << "#Input parameters:" << endl <<  "#" << endl;
  options.write(log1);
  log1 << "#=================================" << endl;
  log1.close();

#ifdef USE_OMP    
  omp_set_dynamic(0);
  omp_set_num_threads(nCores);

  int fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif

  //////////////////////////////////////
  ////// Set up the solid

  double dx            = 1.0/PointsPerSigma;
  int    NumberOfCells = ncopy*ncopy*ncopy;
  double Natoms        = 4*(1-Nvac)*NumberOfCells;
  double L[3]          = {Npoints*dx, Npoints*dx, Npoints*dx};  
  double alatt         = L[0]/ncopy;
  double Density       = 4*(1-Nvac)/(alatt*alatt*alatt);

  cout << "dx            = " << dx << endl;
  cout << "NumberOfCells = " << NumberOfCells << endl;
  cout << "Natoms        = " << Natoms << endl;
  cout << "L[0]          = " << L[0] << endl;
  cout << "alatt         = " << alatt << endl;
  cout << "Density       = " << Density << endl;

  cout << "Check: Natoms/a^3 = " << Natoms/(NumberOfCells*alatt*alatt*alatt) << endl;


  ofstream log("log.dat");
  log << "#Density = " << Density << endl;
  log << "#alatt   = " << alatt << endl;
  log << "#L[0]    = " << L[0] << endl;
  log << "#dx      = " << dx << endl;

  if(L[0] < 2*rcut)
    {
      cout << "L[0] = " << L[0] << " 2*rcut = " << 2*rcut << endl;
      throw std::runtime_error("Simulation cell is too small");
    }
  
  //////////////////////////////////////
  ////// Create potential && effective hsd

  tWF potential(sigma, eps, rcut);
  double hsd = potential.getHSD(kT);

  /////////////////////////////////////
  // Create density object
  
  Solid theDensity(dx, L, alatt, ncopy);
  theDensity.initialize(InitializationInverseWidth, InitializationScaleFac);

  string sd("dump.dat");
  theDensity.writeDensity(sd);
  
  /////////////////////////////////////
  // Create species object

  VDW_Species species(theDensity,hsd,pointsFile, potential, kT);

  
  /////////////////////////////////////
  // DFT object

  DFT_VDW<RSLT> dft(theDensity.Nx(), theDensity.Ny(), theDensity.Nz());
  dft.addSpecies(&species,kT);
  
  if(! infile.empty())
    theDensity.readDensity(infile.c_str());

  string title("Image");
  string file("image.png");
  theDensity.doDisplay(title,file);
  
  cout << "Hard sphere diameter  = " << species.getHSD() << endl;
  //  cout << "Coexisting densities  = " << xliq_coex << " " << xgas_coex << endl;  
  //  cout << "Chemical potential(xliq)/kT = " << mu << endl;
  //  cout << "Chemical potential(xgas)/kT = " << dft.Mu(xgas_coex) << endl;
  //  cout << "beta * Grand free energy per unit volume = " << omega_coex << endl;
  
  string s("log.dat");

  species.setChemPotential(0.0);

  
  fireMinimizer minimizer(dft, Natoms);
  minimizer.setForceTerminationCriterion(forceLimit);
  minimizer.setTimeStep(dt);
  minimizer.setTimeStepMax(dtMax);
  minimizer.setAlphaStart(alpha_start);
  minimizer.setAlphaFac(1.0);
  minimizer.run(s);

   return 1;
}
