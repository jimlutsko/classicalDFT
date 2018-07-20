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
#include "Periodic.h"
#include "Minimizer.h"
/*
Grace *g = NULL;

void display(const Density &density, int calls, double F)
{
  if(g == NULL) return;

  g->deleteDataSet(0);

  for(int i= 0; i < density.Nz(); i++)
    {
      double x = density.getZ(i);
      g->addPoint(x,density.getDensity(density.Nx()/2,density.Ny()/2,i),0);
    }

  stringstream title;
  title << "calls = " << calls <<  " F = " << F;
	
  g->setTitle(title.str().c_str());
  g->redraw();
}

void display1(ConjugateGradients2 &cg)
{
  if(g == NULL) return;

  g->deleteDataSet(0);

  const Density &density = cg.getDensity();

  for(int i= 0; i < density.Nz(); i++)
    {
      double x = density.getZ(i);
      g->addPoint(x,density.getDensity(density.Nx()/2,density.Ny()/2,i),0);
    }

  stringstream title;
  title << "calls = " << cg.getCalls() <<  " F = " << cg.getF();
	
  g->setTitle(title.str().c_str());
  g->redraw();
}
*/
int main(int argc, char** argv)
{
  cout << "Build date  : " <<  (unsigned long) &__BUILD_DATE << endl;
  cout << "Build number: " << (unsigned long) & __BUILD_NUMBER << endl;

  double L[3] = {10,10,10};
  int PointsPerHardSphere = 5;
  int nCores = 6;

  double R = -1;
  double zPos = 0;

  double density = 0.8;

  double alpha = 0.01;
  int MaxIts = 100;

  double epsWall = 1;
  double sigWall = 1;

  double kT = 1;

  double eps   = 1;
  double sigma = 1;
  double rcut  = 3;

  string pointsFile("..//SS31-Mar-2016//ss109.05998");
  string outfile("dump.dat");
  string infile;

  double forceLimit = 1e-4;
  double dt = 1e-3;
  double dtMax = 1;
  double alpha_start = 0.01;

  bool showGraphics = true;
  
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

  options.addOption("OutputFile", &outfile);
  options.addOption("IntegrationPointsFile", &pointsFile);

  options.addOption("ForceTerminationCriterion",&forceLimit);
  options.addOption("TimeStep", &dt);
  options.addOption("TimeStepMax", &dtMax);
  options.addOption("AlphaStart", &alpha_start);

  options.addOption("ShorGraphics", &showGraphics);
  
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


  //////////////////////////////////////
  ////// Create potential

  Potential1 potential(sigma, eps, rcut);

  /////////////////////////////////////
  // Create density object

  Grace *g = NULL;

  if(showGraphics) g = new Grace();
  
  Periodic theDensity(dx, L, g);


  /////////////////////////////////////
  // DFT object

  DFT_VDW<RSLT> dft(theDensity,potential,pointsFile,kT);

  
  double xliq_coex;
  double xgas_coex;
  dft.coexistence(xliq_coex, xgas_coex);
  if(xgas_coex > xliq_coex) swap(xgas_coex, xliq_coex);
  
  double xliq_spin;
  double xgas_spin;
  dft.spinodal(xgas_spin, xliq_spin) ;
  if(xgas_spin > xliq_spin) swap(xgas_spin, xliq_spin);
  
  double omega_coex = dft.Omega(xliq_coex);
  double mu         = dft.Mu(xliq_coex);
  double mu_coex    = mu;
  
  theDensity.initialize(xliq_coex,xgas_coex);

  cout << "Hard sphere diameter  = " << dft.HSD() << endl;
  cout << "Coexisting densities  = " << xliq_coex << " " << xgas_coex << endl;  
  cout << "Chemical potential(xliq)/kT = " << mu << endl;
  cout << "Chemical potential(xgas)/kT = " << dft.Mu(xgas_coex) << endl;
  cout << "beta * Grand free energy per unit volume = " << omega_coex << endl;


  string s("log.dat");

  fireMinimizer_Mu minimizer(dft, theDensity, mu);
  minimizer.setForceTerminationCriterion(forceLimit);
  minimizer.setTimeStep(dt);
  minimizer.setTimeStepMax(dtMax);
  minimizer.setAlphaStart(alpha_start);
  minimizer.run(s);

  double Natoms = theDensity.getNumberAtoms();
  double Omega = minimizer.getF();
  double dOmega = Omega - omega_coex*theDensity.getVolume();
  double SurfaceTension = dOmega/(2*theDensity.Lx()*theDensity.Ly());
  //  SurfaceTension = 1;

  cout << "Final Omega: " << Omega << endl;
  cout << "Excess Omega = " << dOmega << endl;
  cout << "Surface Tension = " << SurfaceTension << endl;

  ofstream log1(s.c_str(),ios::app);
  log1 << "#=================================" << endl << "#" << endl;
  log1 << "#Final Omega: " << Omega << endl;
  log1 << "#Excess Omega = " << dOmega << endl;
  log1 << "#Surface Tension = " << SurfaceTension << endl;

  cout << endl << "Critical Radius" << endl << "Natoms\txgas\txliq\tR" << endl;
  log1 << endl << "#Critical Radius" << endl << "#Natoms\txgas\txliq\tR" << endl;
  for(int i=0;i<100;i++)
    {
      double xgas = xgas_coex+(xgas_spin-xgas_coex)*i/100.0;
      double m = dft.Mu(xgas);
      double xliq = dft.findLiquidFromMu(mu,mu_coex, xliq_coex);
      double V = theDensity.getVolume();

      double domega = dft.Fhelmholtz(xliq)-mu*xliq-dft.Fhelmholtz(xgas)+mu*xgas;
      double R = 2*SurfaceTension/fabs(domega);
      cout << xgas*V << "\t" << xgas << "\t" << xliq << "\t" << R << endl;
      log1 <<"#" <<   xgas*V << "\t" << xgas << "\t" << xliq << "\t" << R << endl;
    }
  log1.close();  
  g->pause();
  g->close();
  return 1;
}
