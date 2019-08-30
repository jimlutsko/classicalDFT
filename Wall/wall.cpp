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

#include <gsl/gsl_sf_lambert.h>

#include "Grace.h"
#include "options.h"
#include "TimeStamp.h"

#include "DFT.h"
#include "Wall.h"
#include "Minimizer.h"


int main(int argc, char** argv)
{
  cout << "Build date  : " <<  (unsigned long) &__BUILD_DATE << endl;
  cout << "Build number: " << (unsigned long) & __BUILD_NUMBER << endl;

  double L[3] = {10,10,10};
  int PointsPerHardSphere = 5;
  int nCores = 6;

  double R = -1;
  double zPos = 0;

  double eta1 = 0.2;
  double eta2 = 0.2;

  double hsd1 = 1;
  double hsd2 = 1;
  
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

  double xmin = 0;
  double xmax = 4;

  double ymin = 0;
  double ymax = 1;  

  
  double forceLimit = 1e-4;
  double dt = 1e-3;
  double dtMax = 1;
  double alpha_start = 0.01;

  bool showGraphics = true;

  double Asurf = 0;
  double rho_surf = -1;
  
  Options options;

  options.addOption("nCores", &nCores);
  options.addOption("BulkEta1", &eta1);
  options.addOption("BulkEta2", &eta2);

  options.addOption("HSD1", &hsd1);
  options.addOption("HSD2", &hsd2);
  
  options.addOption("PointsPerHardSphere", &PointsPerHardSphere);

  options.addOption("kT", &kT);

  options.addOption("ASurfactant", &Asurf);
  options.addOption("RhoSurfactant", &rho_surf);

  options.addOption("eps1",   &eps1);
  options.addOption("sigma1", &sigma1);
  options.addOption("rcut1",  &rcut1);

  options.addOption("eps2",   &eps2);
  options.addOption("sigma2", &sigma2);
  options.addOption("rcut2",  &rcut2);

  options.addOption("eps_surf",   &eps_surf);
  options.addOption("sigma_surf", &sigma_surf);
  options.addOption("rcut_surf",  &rcut_surf);

  options.addOption("Lx", L);
  options.addOption("Ly", L+1);
  options.addOption("Lz", L+2);

  options.addOption("R", &R);
  options.addOption("zPosition", &zPos);

  options.addOption("sigWall", &sigWall);
  options.addOption("epsWall", &epsWall);

  options.addOption("OutputFile", &outfile);
  options.addOption("IntegrationPointsFile", &pointsFile);

  options.addOption("ForceTerminationCriterion",&forceLimit);
  options.addOption("TimeStep", &dt);
  options.addOption("TimeStepMax", &dtMax);
  options.addOption("AlphaStart", &alpha_start);

  options.addOption("InputFile", &infile);
  
  options.addOption("ShorGraphics", &showGraphics);

  options.addOption("Xmin", &xmin);
  options.addOption("Xmax", &xmax);
  options.addOption("Ymin", &ymin);
  options.addOption("Ymax", &ymax);

  
  options.read(argc, argv);

  ofstream log1("log.dat");
  TimeStamp ts;
  log1 << "# " << ts << endl;
  log1 << "#=================================" << endl << "#" << endl;
  log1 << "#Input parameters:" << endl <<  "#" << endl;
  options.write(log1);
  log1 << "#=================================" << endl;


  double dx = 1.0/PointsPerHardSphere;

  if(L[0] < 0) L[0] = dx;
  if(L[1] < 0) L[1] = dx;
  if(L[2] < 0) L[2] = dx;  

  
#ifdef USE_OMP    
  omp_set_dynamic(0);
  omp_set_num_threads(nCores);

  int fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif

  Grace *g = (showGraphics ? new Grace() : NULL);

  g->setLimits(xmin,xmax,ymin,ymax);
  

  //////////////////////////////////////
  ////// Create potential && effective hsd

  LJ potential1(sigma1, eps1, rcut1);
  LJ potential2(sigma2, eps2, rcut2);

  if(hsd1 < 0) hsd1 = potential1.getHSD(kT);
  if(hsd2 < 0) hsd2 = potential2.getHSD(kT);

  double density1 = 6*eta1/(M_PI*hsd1*hsd1*hsd1);
  double density2 = 6*eta2/(M_PI*hsd2*hsd2*hsd2);

  cout << "Density 1 = " << density1 << endl;
  cout << "Density 2 = " << density2 << endl;

  log1 << "#Density 1 = " << density1 << endl;
  log1 << "#Density 2 = " << density2 << endl;

  cout << "HSD 1 = " << hsd1 << endl;
  cout << "HSD 2 = " << hsd2 << endl;

  log1 << "#HSD 1 = " << hsd1 << endl;
  log1 << "#HSD 2 = " << hsd2 << endl;


  
  /////////////////////////////////////
  // Create density objects
  Wall theDensity1(dx, L, g, 0.001, 1, hsd1);
  Wall theDensity2(dx, L, g, 0.001, 1, hsd2);

  /////////////////////////////////////
  // Create species objects

  FMT_Species species1(theDensity1,hsd1,pointsFile);
  FMT_Species species2(theDensity2,hsd2,pointsFile);

  /////////////////////////////////////
  // Create FMT object
  
  WhiteBearI fmt;
  
  /////////////////////////////////////
  // DFT object

  DFT dft(&species1);
  dft.addSpecies(&species2);
  dft.addHardCoreContribution(&fmt);

  
  /////////////////////////////////////////////////////
  // The thermodynamics

  vector<double> densities;
  densities.push_back(density1);
  densities.push_back(density2);

  double mu1 = dft.Mu(densities,0);
  double mu2 = dft.Mu(densities,1);

  cout << "Omega/(V kT) = " << dft.Omega(densities) << endl;
  cout << "mu1 = " << mu1 << endl;
  cout << "mu2 = " << mu2 << endl;
  log1 << "#Omega/(V kT) = " << dft.Omega(densities) << endl;  
  log1 << "#mu1 = " << mu1 << endl;
  log1 << "#mu2 = " << mu2 << endl;

  theDensity1.initialize(density1, density1);
  theDensity2.initialize(density2, density2);

  double N = theDensity2.getNumberAtoms();

  ///////////////////////////////////////////////
  // Fix the mass of the surfactant species.
  //species2.setFixedMass(N);

  
  
  /////////////////Check the thermodynamics
  double eps = 1e-4;

  densities[0] = density1*(1+eps);
  double fp = dft.Fhelmholtz(densities);

  densities[0] = density1*(1-eps);
  double fm = dft.Fhelmholtz(densities);

  densities[0] = density1;

  cout << " mu1 = " << mu1 << " numeric: " << (fp-fm)/(2*eps*density1) << endl;

  densities[1] = density2*(1+eps);
  fp = dft.Fhelmholtz(densities);

  densities[1] = density2*(1-eps);
  fm = dft.Fhelmholtz(densities);

  densities[1] = density2;
  
  cout << " mu2 = " << mu2 << " numeric: " << (fp-fm)/(2*eps*density2) << endl;
  

  
  species1.setChemPotential(mu1);
  species2.setChemPotential(mu2);
  
  if(! infile.empty())
    theDensity1.readDensity(infile.c_str());

  log1.close();
    
  string s("log.dat");

  fireMinimizer_Mu minimizer(dft);
  minimizer.setForceTerminationCriterion(forceLimit);
  minimizer.setTimeStep(dt);
  minimizer.setTimeStepMax(dtMax);
  minimizer.setAlphaStart(alpha_start);
  minimizer.setAlphaFac(1.0);
  minimizer.run(s);

  double Natoms = theDensity1.getNumberAtoms();
  double Omega = minimizer.getF();
  double dOmega = Omega-L[0]*L[1]*(L[2]-2)*dft.Omega(densities);
  double SurfaceTension = dOmega/(L[0]*L[1]);

  cout << "dft.Omega = " << dft.Omega(densities) << endl;
  cout << "Final Omega: " << Omega << endl;
  cout << "Excess Omega = " << dOmega << endl;
  cout << "Surface Tension = " << SurfaceTension << endl;

  ofstream log2(s.c_str(),ios::app);
  log2 << "#dft.Omega = " << dft.Omega(densities) << endl;  
  log2 << "#=================================" << endl << "#" << endl;
  log2 << "#Final Omega: " << Omega << endl;
  log2 << "#Excess Omega = " << dOmega << endl;
  log2 << "#Surface Tension = " << SurfaceTension << endl;


  cout << "dF1.inf_norm() = " <<  species1.getDF().inf_norm() << " dF1.inf_norm()/density_.dV() = " << species1.getDF().inf_norm()/theDensity1.dV() << endl;
  cout << "dF2.inf_norm() = " <<  species2.getDF().inf_norm() << " dF2.inf_norm()/density_.dV() = " << species2.getDF().inf_norm()/theDensity1.dV() << endl;

  cout << "species1.get_convergence_monitor() = " << species1.get_convergence_monitor() << endl;
  cout << "species2.get_convergence_monitor() = " << species2.get_convergence_monitor() << endl;  


  
  char cc;
  cin >> cc;

  ///////////////// CHECK the derivatives


  g->deleteDataSet(0);
  g->deleteDataSet(1);
  g->deleteDataSet(2);
  g->deleteDataSet(3);

  // We need this call because currently, the forces are dF/dx and not dF/d_rho
  double F = dft.calculateFreeEnergyAndDerivatives(false);
  double epsx = 1e-4;

  DFT_Vec dF1 = species1.getDF();
  DFT_Vec dF2 = species2.getDF();

  int ix = theDensity1.Nx()/2;
  int iy = theDensity1.Ny()/2;

  for(int iz=0;iz<theDensity1.Nz();iz++)
    {
      long p = theDensity1.pos(ix,iy,iz);
      
      double d = theDensity1.getDensity(p);

      if(d < 1e-6) continue;
      
      theDensity1.set_Density_Elem(p,d*(1+epsx));
      double FP = dft.calculateFreeEnergyAndDerivatives(false);

      theDensity1.set_Density_Elem(p,d*(1-epsx));      
      double FM = dft.calculateFreeEnergyAndDerivatives(false);

      theDensity1.set_Density_Elem(p,d);      

      g->addPoint(iz, dF1.get(p), 0);
      g->addPoint(iz, (FP-FM)/(2*epsx*d) , 1);

      p = theDensity2.pos(ix,iy,iz);
      
      d = theDensity2.getDensity(p);

      if(d < 1e-6) continue; // eliminates regions where V->inf and rho->zero
      
      theDensity2.set_Density_Elem(p,d*(1+epsx));
      FP = dft.calculateFreeEnergyAndDerivatives(false);

      theDensity2.set_Density_Elem(p,d*(1-epsx));      
      FM = dft.calculateFreeEnergyAndDerivatives(false);

      theDensity2.set_Density_Elem(p,d);      

      g->addPoint(iz, dF2.get(p), 2);
      g->addPoint(iz, (FP-FM)/(2*epsx*d) , 3);

      g->redraw(1);      

    }

  g->pause();


  g->close();
  return 1;
}
