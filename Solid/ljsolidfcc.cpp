#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
//#include <time.h>

using namespace std;

#ifdef USE_OMP
#include <omp.h>
#endif

#include "options.h"
#include "TimeStamp.h"

#include "DFT.h"
#include "SolidFCC.h"
#include "Minimizer.h"
#include "Log.h"
#include "myColor.h"

void check(Density &density, DFT &dft, Interaction &i);

double L[3] = {10,10,10};
int nCores = 6;

double hsd1 = -1;
  
double kT = 1;

string pointsFile("..//SS31-Mar-2016//ss109.05998");
string outfile("dump.dat");
string infile;
  
double forceLimit = 1e-4;
double dt = 1e-3;
double dtMax = 1;
double alpha_start = 0.01;

bool showGraphics = true;

int maxSteps = -1;

double prefac = 1;
double Nvac = 0;
double dx = 0.1;
double alpha = 50;
double alphaFac = 1;

bool doCalc(double Mu, int Npoints, Potential1 &potential, Log &log, double &F, double &Cvac);
bool doUniform(double Mu, int Npoints, Potential1 &potential1, Log& log, double &Fret, double &Dret);

int main(int argc, char** argv)
{
  double eps1   = 1;
  double sigma1 = 1;
  double rcut1  = 3;

  int Npoints = 1;  
  double Mu = -1.0;
  Options options;

  options.addOption("prefac", &prefac);
  options.addOption("Nval", &Nvac);
  options.addOption("Npoints", &Npoints);
  options.addOption("alpha", &alpha);
  options.addOption("Dx", &dx);
  
  options.addOption("nCores", &nCores);

  options.addOption("kT", &kT);
  options.addOption("Mu",   &Mu);

  
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
  options.addOption("AlphaFac", &alphaFac);
  options.addOption("MaxSteps", &maxSteps);

  options.addOption("InputFile", &infile);
  
  options.addOption("ShowGraphics", &showGraphics);
  
  options.read(argc, argv);

  Log log("log.dat");
  log << myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;

  log << myColor::RED << myColor::BOLD << "Input parameters:" << myColor::RESET << endl <<  "#" << endl;

  options.write(log);
  log <<  myColor::GREEN << "=================================" <<  myColor::RESET << endl;

  
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
  
  ofstream of("scan.dat");
  of.close();
    
  for(Mu = -0.25; Mu > -3; Mu -= 0.25)
    {
      double F;
      double Cvac;
      double D = 1.0;
      bool ret = doUniform(Mu, Npoints, potential1, log, F,D);

      ofstream of("scan.dat", ios::app);      
      of << "#Mu = " << Mu << " Funiform/(kTV) = " << F << " Duniform = " << D << endl;

      D = 1e-6;
      ret = doUniform(Mu, Npoints, potential1, log, F,D);

      of << "#Mu = " << Mu << " Funiform/(kTV) = " << F << " Duniform = " << D << endl;
      of << "#"
	 << setw(15) <<"Npoints"
	 << setw(15) <<"LattPar."
	 << setw(15) <<"FF/(kT V)"
	 << setw(15) <<"Cvacancy"
	 << endl;  
      of.close();      

      bool bstarted = false;

      for(Npoints = 61; Npoints < 140; Npoints++)
	{
	  bool ret = doCalc(Mu, Npoints, potential1, log, F, Cvac);
	  if(ret)
	    {
	      bstarted = true; // found beginning of good range

	      ofstream of("scan.dat", ios::app);      
	      of << setw(15) << Npoints
		 << setw(15) << Npoints*dx
		 << setw(15) << F
		 << setw(15) << Cvac
		 << endl;
	      of.close();
	      
	    }
	  else if(bstarted) break; // exhausted the good range
	}
      ofstream of1("scan.dat", ios::app);
      of1 << "#" << endl;
      of1.close();      
    }

  log << "Finished" << endl;
  
  return 1;

}

bool doUniform(double Mu, int Npoints, Potential1 &potential1, Log& log, double &Fret, double &Dret)
{
  int ncopy = 1;

  /////////////////////////////////////
  // Create objects

  double L[] = {ncopy*Npoints*dx, ncopy*Npoints*dx, ncopy*Npoints*dx};
  SolidFCC theDensity1(dx, L, hsd1);  
  //  theDensity1.initialize2(/*alpha*/ 100, alatt, ncopy, prefac);
  FMT_Species species1(theDensity1,hsd1,pointsFile,Mu,Npoints);
  Interaction i1(species1,species1,potential1,kT,log, pointsFile);
  RSLT fmt;

  DFT dft(&species1);

  dft.addHardCoreContribution(&fmt);  
  dft.addInteraction(&i1);
  theDensity1.initializeUniform(Dret);
  
  log <<  myColor::GREEN << "=================================" <<  myColor::RESET << endl;
  
  fireMinimizer2 minimizer(dft,log);
  minimizer.setForceTerminationCriterion(forceLimit);
  minimizer.setTimeStep(dt);
  minimizer.setTimeStepMax(dtMax);
  minimizer.setAlphaStart(alpha_start);
  minimizer.setAlphaFac(alphaFac);
  minimizer.run(maxSteps); 

  double Natoms = theDensity1.getNumberAtoms();
  double Omega = minimizer.getF();
  //  double dOmega = Omega-L[0]*L[1]*L[2]*dft.Omega(densities);

  log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
  log << "Final Uniform Natoms               : " << Natoms << endl;
  log << "Final Omega                : " << Omega << endl;

  Dret = Natoms/theDensity1.getVolume();
  Fret = Omega/theDensity1.getVolume();
  return true;
}

bool doCalc(double Mu, int Npoints, Potential1 &potential1, Log& log, double &Fret, double &Cvac)
{
  int ncopy = 1;

  /////////////////////////////////////
  // Create objects

  double L[] = {ncopy*Npoints*dx, ncopy*Npoints*dx, ncopy*Npoints*dx};
  SolidFCC theDensity1(dx, L, hsd1);  
  FMT_Species species1(theDensity1,hsd1,pointsFile,Mu,Npoints);
  Interaction i1(species1,species1,potential1,kT,log, pointsFile);
  RSLT fmt;

  DFT dft(&species1);

  dft.addHardCoreContribution(&fmt);  
  dft.addInteraction(&i1);
  
  /////////////////////////////////////////////////////
  // Report
  log << "Lx  = " << L[0] << endl;
  log << "HSD = " << hsd1 << endl;

  //  if(! infile.empty())
  //    theDensity1.readDensity(infile.c_str());

  ///////////////////////////////////////////////
  // Fix the chemical potential

  //  species1.setChemPotential(Mu);
    
  /////////////////////////////////////////////////////
  // The thermodynamics
  /*
  double N = theDensity1.getNumberAtoms();
  
  vector<double> densities;
  densities.push_back(N/(L[0]*L[1]*L[2]));

  double mu1 = dft.Mu(densities,0);
  
  log << "Omega/(V kT) = " << dft.Omega(densities) << endl;  
  log << "mu liq = " << mu1 << endl;
  log << "F liq = " << dft.Fhelmholtz(densities)*L[0]*L[1]*L[2] << endl;
  log << "Omega Liq (at Mu = " << Mu << ") = " << (dft.Fhelmholtz(densities)-Mu*densities[0])*L[0]*L[1]*L[2] << endl;
  log << "a_vdw = " << i1.getVDWParameter() << endl;
  */

  log << " " << endl;
  log << "Determining gaussian approximation" << endl;
  log << " " << endl;
  
  double f = 0;
  double f_old = 0;
  double alf_old = -1;
  double prefac_old = 0;
  bool isDecending = false;
  bool firstIteration = true;
  for(double alf = 30; alf < 1000; alf += 10)
    {
      bool bSuccess = false;
      while(!bSuccess)
	{
	  try {
	    theDensity1.initialize(alf, ncopy, prefac);
	    f = dft.calculateFreeEnergyAndDerivatives(false);
	    bSuccess = true;
	  } catch (...) {
	    prefac *= 0.9999;
	    log << "\t lowering prefac to " << prefac << endl;
	  }
	}
      log << "alf = " << alf << " " << f << endl;

      if(!firstIteration)
	if(f < f_old) isDecending = true;
      
      if(isDecending)
	if(f > f_old) break;

      f_old = f;	        
      alf_old = alf;
      prefac_old = prefac;

      firstIteration = false;
    }
  if(isDecending)
    log << "Found: alf = " << alf_old << " prefac = " << prefac_old << " f = " << f_old << endl;
  else {
    log << "No minimum found with Gaussians" << endl;
    return false;
  }

  theDensity1.initialize(alf_old, ncopy, prefac_old);
  
  //check(theDensity1, dft,i1);

  log <<  myColor::GREEN << "=================================" <<  myColor::RESET << endl;
  
  fireMinimizer2 minimizer(dft,log);
  minimizer.setForceTerminationCriterion(forceLimit);
  minimizer.setTimeStep(dt);
  minimizer.setTimeStepMax(dtMax);
  minimizer.setAlphaStart(alpha_start);
  minimizer.setAlphaFac(alphaFac);
  minimizer.run(maxSteps); 


  double Natoms = theDensity1.getNumberAtoms();
  double Omega = minimizer.getF();
  //  double dOmega = Omega-L[0]*L[1]*L[2]*dft.Omega(densities);

  log <<  myColor::GREEN << "=================================" << myColor::RESET << endl << "#" << endl;
  log << "Final Natoms               : " << Natoms << endl;
  log << "Final Vacancy Concentration: " << (4.0-Natoms)/theDensity1.getVolume() << endl;  
  log << "Final Omega                : " << Omega << endl;

  Fret = Omega/theDensity1.getVolume();
  Cvac = (4.0-Natoms)/theDensity1.getVolume();
  return true;
}


void check(Density &density, DFT &dft, Interaction &i1)
{

  int Nx = density.Nx();
  int Ny = density.Ny();
  int Nz = density.Nz();
  
  // Find largest density
  int jx = 0;
  int jy = 0;
  int jz = 0;
  double dmax = 0;

  for(int ix = 0; ix < Nx; ix++)
    for(int iy = 0; iy < Ny; iy++)
      for(int iz = 0; iz < Nz; iz++)
	if(density.getDensity(ix,iy,iz) > dmax)
	  {
	    dmax = density.getDensity(ix,iy,iz);
	    jx = ix;
	    jy = iy;
	    jz = iz;
	  }


  double  F = dft.calculateFreeEnergyAndDerivatives(false);   
  DFT_Vec dF;
  dF.zeros(dft.getDF(0).size());
  dF.set(dft.getDF(0));
    
  for(int iz = jz; iz < Nz; iz++)
    {
      double eps = 1e-6;

      double d = density.getDensity(jx,jy,iz);

      double test = dF.get(density.pos(jx,jy,iz));


      double test1 = F+eps*d*test+0.5*eps*eps*d*d*i1.getW(0);
      double test2 = F-eps*d*test+0.5*eps*eps*d*d*i1.getW(0);
     
      

      double d1 = d+eps; //*d;
      density.set_Density_Elem(jx,jy,iz,d1);
      double fp = dft.calculateFreeEnergyAndDerivatives(false);
      cout << "Test1 = " << test1 << endl;

      double d2 = d-eps; //*d;
      density.set_Density_Elem(jx,jy,iz,d2);
      double fm = dft.calculateFreeEnergyAndDerivatives(false);
      cout << "Test2 = " << test2 << endl;
      
      density.set_Density_Elem(jx,jy,iz,d);

      double direct = i1.checkCalc(jx,jy,iz);


      
      cout << jx << " " << jy << " " << iz << " " << d << " " << ((fp-fm)/(d1-d2)) << " " << dF.get(density.pos(jx,jy,iz)) << " direct = " << direct << endl;
      
    }
  cout << "Finished " << Nz << endl;
  exit(0);  

}
