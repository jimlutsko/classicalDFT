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
#include "Log.h"

#include <gsl/gsl_poly.h>

double a = 0;
double aa = 0;
double d = 0;
double x2 = 0;

double P(double x)
{
  double e = M_PI*x*d*d*d/6;
  double P = x*(1+e+e*e-e*e*e)*pow(1.0-e,-3)+ a*x*x;
  return P;
}

double f(double x)
{
  double e = M_PI*x*d*d*d/6;
  double f = log(x)-1 + e*(4.0-3.0*e)*pow(1-e,-2) + a*x;
  return f;
}

double f1(double x)
{
  // first, get the pressure multiplied by the hs factor
  double e = M_PI*x*d*d*d/6;
  double aa = a*6/(M_PI*d*d*d);
  
  double P1 = (M_PI*d*d*d/6)*P(x);
  
  // second, the coefficients of the polynomial equation we need to solve
  double c[6] = {-P1,1+3*P1,1-3*P1+aa,1+P1-3*aa,-1+3*aa,-aa};

  //divide out the known root
    for(int i=4;i>=0;i--)
    c[i] += c[i+1]*e;

  // now, find roots and print
  double z[10];  
  gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc(5);
  gsl_poly_complex_solve (c+1, 5, w, z);
  gsl_poly_complex_workspace_free(w);

  double emin = 0;
  double immin = 1e30;
  for(int i = 0; i < 4; i++)
    {
      double e = z[2*i];
      double ei = z[2*i+1];
      double dPde = (1+4*e+4*e*e-4*e*e*e+e*e*e*e)*pow(1-e,-4)+2*aa*e;
      
      if(dPde > 0)
	if(fabs(ei) < immin) { immin = fabs(ei); emin = e;}
    }

  e = emin;
  x2 = 6*emin/(M_PI*d*d*d);
  
  double mu1 = f(x)+P(x)/x;
  double mu2 = f(x2)+P(x2)/x;
  
  cout << "emin = " << emin << " P = " << P(x2) << " dmu = " << mu1-mu2 << endl;
  
  return 0;
}


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
  
  options.addOption("ShowGraphics", &showGraphics);
  
  options.read(argc, argv);

  Log log("log.dat");
  TimeStamp ts;
  log << ts << endl;
  log << "=================================" << endl << " " << endl;
  log << " Input parameters:" << endl <<  " " << endl;
  options.write(log);
  log << " =================================" << endl;


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

  //////////////////////////////////////
  ////// Create potential && effective hsd

  LJ potential1(sigma1, eps1, rcut1);
  LJ potential2(sigma2, eps2, rcut2);

  if(hsd1 < 0) hsd1 = potential1.getHSD(kT);
  if(hsd2 < 0) hsd2 = potential2.getHSD(kT);

  
  /////////////////////////////////////
  // Create density objects

  Grace *g = (showGraphics ? new Grace() : NULL);
  
  Wall theDensity1(dx, L, g, 0.001, 1, hsd1);
  Wall theDensity2(dx, L, g, 0.001, 1, hsd2);

  /////////////////////////////////////
  // Create the species objects
  try {
    Species species1(theDensity1);
    Species species2(theDensity2);

    //    FMT_Species species1(theDensity1,hsd1,pointsFile);
    //    FMT_Species species2(theDensity2,hsd2,pointsFile);    

    Interaction i1(species1,species1,potential1,kT);
    Interaction i2(species2,species2,potential2,kT);

    g->close();
    delete g;


    log << "Hsd = " << hsd1 << endl;
    log << "a vdw = " << i1.getVDWParameter() << endl;
    
    a = i1.getVDWParameter();
    d = hsd1;
    while(1)
      {
	double x;
	cin >> x;
	double r = f1(x);
      }
    exit(0);
  

    /////////////////////////////////////
    // Create the hard-sphere object
    RSLT fmt;
  
  
    ///////////////////////////////////////////////////
    // The surfactant potential
    //  LJ surfactant_potential(sigma_surf, eps_surf, rcut_surf);

    /////////////////////////////////////
    // DFT object
    //  DFT_VDW<RSLT> dft(&species1);

    //  DFT_FMT<RSLT> dft(&species1);

    DFT dft(&species1);
    dft.addSpecies(&species2);
    dft.addHardCoreContribution(&fmt);  
    dft.addInteraction(&i1);
    dft.addInteraction(&i2);

    /////////////////////////////////////////////////////
    // Report
    double density1 = 6*eta1/(M_PI*hsd1*hsd1*hsd1);
    double density2 = 6*eta2/(M_PI*hsd2*hsd2*hsd2);

    log << "Density 1 = " << density1 << endl;
    log << "Density 2 = " << density2 << endl;

    log << "HSD 1 = " << hsd1 << endl;
    log << "HSD 2 = " << hsd2 << endl;

  
    /////////////////////////////////////////////////////
    // The thermodynamics

    vector<double> densities;
    densities.push_back(density1);
    densities.push_back(density2);

    double mu1 = dft.Mu(densities,0);
    double mu2 = dft.Mu(densities,1);

    log << "Omega/(V kT) = " << dft.Omega(densities) << endl;  
    log << "mu1 = " << mu1 << endl;
    log << "mu2 = " << mu2 << endl;

    theDensity1.initialize(density1, density1);
    theDensity2.initialize(density2, density2);

    double N = theDensity2.getNumberAtoms();

    ///////////////////////////////////////////////
    // Fix the mass of the surfactant species.
    //species2.setFixedMass(N);

  
    /*
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
    */

    /*
///////////////// CHECK the derivatives

double F = dft.calculateFreeEnergyAndDerivatives(true);
DFT_Vec dF1 = species1.getDF();
DFT_Vec dF2 = species2.getDF();
double epsx = 1e-4;

int ix = theDensity1.Nx()/2;
int iy = theDensity1.Ny()/2;
double dV = theDensity1.dV();
  
for(int iz=0;iz<theDensity1.Nz();iz++)
{
long p = theDensity2.pos(ix,iy,iz);
      
double d = theDensity2.getDensity(p);

if(d < 1e-6) continue;
      
theDensity2.set_Density_Elem(p,d*(1+epsx));
double FP = dft.calculateFreeEnergyAndDerivatives(true);

theDensity2.set_Density_Elem(p,d*(1-epsx));      
double FM = dft.calculateFreeEnergyAndDerivatives(true);

theDensity2.set_Density_Elem(p,d);      

cout << "iz = " << iz << " dF = " << dF2.get(p) << " " << (FP-FM)/(2*epsx*d) << "  difference: " << fabs(dF2.get(p)-(FP-FM)/(2*epsx*d)) << endl;
}

g->close();
exit(0);
    */
    
  species1.setChemPotential(mu1);
  species2.setChemPotential(mu2);
  
  if(! infile.empty())
    theDensity1.readDensity(infile.c_str());

  fireMinimizer_Mu minimizer(dft,log);
  minimizer.setForceTerminationCriterion(forceLimit);
  minimizer.setTimeStep(dt);
  minimizer.setTimeStepMax(dtMax);
  minimizer.setAlphaStart(alpha_start);
  minimizer.setAlphaFac(1.0);
  minimizer.run();

  double Natoms = theDensity1.getNumberAtoms();
  double Omega = minimizer.getF();
  double dOmega = Omega-L[0]*L[1]*(L[2]-2)*dft.Omega(densities);
  double SurfaceTension = dOmega/(L[0]*L[1]);

  log << "dft.Omega = " << dft.Omega(densities) << endl;  
  log << "=================================" << endl << " " << endl;
  log << "Final Omega: " << Omega << endl;
  log << "Excess Omega = " << dOmega << endl;
  log << "Surface Tension = " << SurfaceTension << endl;
    

  g->redraw();
  g->pause();
  } catch(const std::exception &e) {
    log << " " << endl;
    log << e.what() << endl;
    log << " " << endl;
    log << "ABORTING ..." << endl;
    log << " " << endl;
  }

  g->close();
  delete g;

  return 1;
}
