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
#include "Integrator.h"

#include "DFT.h"
#include "Periodic.h"
#include "Minimizer.h"


void solve(DFT_VDW<RSLT> & dft_, double xliq, double xvap, double rho_surf_, double A, Grace *g);

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

  double Asurf = 0;
  double rho_surf = -1;
  
  Options options;

  options.addOption("nCores", &nCores);
  options.addOption("BulkDensity", &density);
  options.addOption("PointsPerHardSphere", &PointsPerHardSphere);

  options.addOption("kT", &kT);

  options.addOption("ASurfactant", &Asurf);
  options.addOption("RhoSurfactant", &rho_surf);

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

  options.addOption("InputFile", &infile);
  
  options.addOption("ShorGraphics", &showGraphics);
  
  options.read(argc, argv);

  ofstream log1("log.dat");
  TimeStamp ts;
  log1 << "# " << ts << endl;
  log1 << "#=================================" << endl << "#" << endl;
  log1 << "#Input parameters:" << endl <<  "#" << endl;
  options.write(log1);
  log1 << "#=================================" << endl;
  log1.close();

  double dx = 1.0/PointsPerHardSphere;

#ifdef USE_OMP    
  omp_set_dynamic(0);
  omp_set_num_threads(nCores);

  int fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif


  //////////////////////////////////////
  ////// Create potential

  LJ potential(sigma, eps, rcut);

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

  theDensity.setDFT(&dft);
  
  theDensity.initialize(xliq_coex,xgas_coex);

  if(! infile.empty())
    theDensity.readDensity(infile.c_str());


  
  cout << "Hard sphere diameter  = " << dft.HSD() << endl;
  cout << "Coexisting densities  = " << xliq_coex << " " << xgas_coex << endl;  
  cout << "Chemical potential(xliq)/kT = " << mu << endl;
  cout << "Chemical potential(xgas)/kT = " << dft.Mu(xgas_coex) << endl;
  cout << "beta * Grand free energy per unit volume = " << omega_coex << endl;


  //  for(Asurf = 0.01; Asurf < 8; Asurf += 0.01)
    {
  
  dft.setSurfactant(rho_surf,Asurf);

  string s("log.dat");

  //  solve(dft, xliq_coex, xgas_coex, rho_surf, Asurf, g);
  //  exit(0);
  
  fireMinimizer_Mu minimizer(dft, theDensity, mu);
  minimizer.setForceTerminationCriterion(forceLimit);
  minimizer.setTimeStep(dt);
  minimizer.setTimeStepMax(dtMax);
  minimizer.setAlphaStart(alpha_start);
  minimizer.setAlphaFac(1.0);
  minimizer.run(s);

  double Natoms = theDensity.getNumberAtoms();
  double Omega = minimizer.getF();
  double dOmega = Omega - omega_coex*theDensity.getVolume();
  double SurfaceTension = dOmega/(2*theDensity.Lx()*theDensity.Ly());

  cout << "Final Omega: " << Omega << endl;
  cout << "Excess Omega = " << dOmega << endl;
  cout << "Surface Tension = " << SurfaceTension << endl;

  ofstream log2(s.c_str(),ios::app);
  log2 << "#=================================" << endl << "#" << endl;
  log2 << "#Final Omega: " << Omega << endl;
  log2 << "#Excess Omega = " << dOmega << endl;
  log2 << "#Surface Tension = " << SurfaceTension << endl;
    
  /*
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
  */
  /*  
  DFT_Vec dF; dF.zeros(theDensity.Ntot());
  DFT_Vec dF0; dF0.zeros(theDensity.Ntot());
  
  double F = dft.calculateFreeEnergyAndDerivatives(theDensity, mu, dF);
  cout << "Recalculated F = " << F << endl;

  
  eps = 0.0001;
  
  
  for(int L=0;L<theDensity.Nz(); L++)
    {
      double x = theDensity.getDensity(0,0,L);

      theDensity.set_Density_Elem(0,0,L,x*(1+eps));
      double Fp0 = dft.calculateFreeEnergyAndDerivatives(theDensity, mu*0, dF0);
      double Fp = dft.calculateFreeEnergyDerivativesSurf(theDensity, Asurf, rho_surf, dF);

      theDensity.set_Density_Elem(0,0,L,x*(1-eps));
      double Fm0 = dft.calculateFreeEnergyAndDerivatives(theDensity, mu*0, dF0);
      double Fm = dft.calculateFreeEnergyDerivativesSurf(theDensity, Asurf, rho_surf, dF);

      theDensity.set_Density_Elem(0,0,L,x);
      double F0 = dft.calculateFreeEnergyAndDerivatives(theDensity, mu*0, dF0);
      F = dft.calculateFreeEnergyDerivativesSurf(theDensity, Asurf, rho_surf, dF);


      cout << "Numeric = " << (Fp0-Fm0)/(2*eps*x) << " analy = " << dF0.get(L) << endl;
      cout << "Numeric Surf= " << (Fp-Fm)/(2*eps*x) << " analy = " << dF.get(L) << endl;
  
  
      cout << "Surfactant F = " << F << endl << endl;
     
    }
  */
  log2.close();
    }
  g->pause();
  g->close();
  return 1;
}

static double omega0;
static double rho_surf;
static double B;
static DFT_VDW<RSLT> *dft;
static double mu;

double S1(double x, void*)
{
  double arg = dft->Fhelmholtz(x)-mu*x-omega0-rho_surf;
  double f   = arg+(1/B)*gsl_sf_lambert_W0(B*rho_surf*exp(-B*arg));
  return 1.0/sqrt(f);
}

void solve(DFT_VDW<RSLT> & dft_, double xliq, double xvap, double rho_surf_, double A, Grace *g)
{
  dft = &dft_;
  rho_surf = rho_surf_;
  omega0 = dft_.Omega(xliq);
  mu = dft_.Mu(xliq);

  B = (M_PI*M_PI/18)*A*pow(dft_.HSD(),6);

  cout << "alpha = " << B << endl;
  
  Integrator1 I1(S1);

  double x0 = (xliq+xvap)/2;
  
  double z = 0;
  double dx = (xliq-1e-4-x0)/100;

  vector<double> zs;
  vector<double> xs;


  for(double x = x0; x-dx > xvap; x -= dx)
    {
      z -= I1.integrateFinite(x-dx,x);
      zs.push_back(z);
      xs.push_back(x-dx);      
    }
  std::reverse(std::begin(zs), std::end(zs));
  std::reverse(std::begin(xs), std::end(xs));

  z = 0;
  zs.push_back(z);
  xs.push_back(x0);
  
  for(double x = x0; x+dx < xliq; x += dx)
    {
      z += I1.integrateFinite(x,x+dx);
      zs.push_back(z);
      xs.push_back(x+dx);      
    }

  for(int i=0;i<xs.size();i++)
    g->addPoint(zs[i], xs[i]);

  g->redraw();

  return;  
}
