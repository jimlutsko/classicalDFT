#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>


using namespace std;

#ifdef USE_OMP
#include <omp.h>
#endif

#include "options.h"
#include "TimeStamp.h"

#include "DFT.h"
#include "VDW1.h"
#include "SolidFCC.h"
#include "Minimizer.h"
#include "Log.h"
#include "myColor.h"

void check(Density &density, DFT &dft, Interaction &i);
double gaussianEval(double alf, double bmu, SolidFCC& theDensity, DFT& dft, double &prefac, Log &log);
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
double alphaFac = 1;

double cutoff = -1;


bool showGraphics = true;

int maxSteps = -1;
int Ngauss = 5;

//double prefac = 1;
double prefac = 0.995;
double dx = 0.1;

bool findGaussian(SolidFCC& theDensity1, double bmu, DFT& dft, Log& log, double &Fgau, double &Agau);
bool findMinimum(SolidFCC& theDensity1, DFT& dft, Log& log, double &Fret, double &Cvac);
bool doUniform(SolidFCC& theDensity1, DFT& dft, Log& log, double &Fret, double &Dret);


//double bmu = 1e-20;

class myMin : public fireMinimizer2
{
public:
  myMin(DFT &dft, ostream &log, double t = 1) :  fireMinimizer2(dft, log) { threshold_ = t;}
  virtual void draw_after()
  {
    cout << "dt_ = " << dt_ << " dt_max_ = " << dt_max_ << " rms_force = " << rms_force_ << " rms_velocity = " << vnorm_ << " fmax = " << fmax_ << endl;    
  }
};


void do_N_first(int Npoints_min, int Npoints_max, double Mu_min, double Mu_max, double Mu_step, bool GaussianOnly, Potential1& potential1, Log &log);
void do_Mu_first(int Npoints_min, int Npoints_max, double Mu_min, double Mu_max, double Mu_step, bool GaussianOnly, Potential1& potential1, Log &log);

int main(int argc, char** argv)
{
  double eps1   = 1;
  double sigma1 = 1;
  double rcut1  = 3;

  int Npoints_min = 120;
  int Npoints_max = 120;    
  double Mu_min  =  0.25;
  double Mu_max  = -0.25;
  double Mu_step = 0.25;

  bool GaussianOnly = false;
  
  Options options;

  //  options.addOption("prefac", &prefac);
  options.addOption("Npoints_min", &Npoints_min);
  options.addOption("Npoints_max", &Npoints_max);
  options.addOption("Dx", &dx);
  options.addOption("Ngauss", &Ngauss);
  
  options.addOption("nCores", &nCores);

  options.addOption("kT", &kT);
  options.addOption("Mu_min",   &Mu_min);
  options.addOption("Mu_max",   &Mu_max);
  options.addOption("Mu_step",  &Mu_step);

  options.addOption("AttractiveCutoff", &cutoff);
  
  options.addOption("eps1",   &eps1);
  options.addOption("sigma1", &sigma1);
  options.addOption("rcut1",  &rcut1);

  options.addOption("GaussianOnly", &GaussianOnly);

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
  log << myColor::GREEN << "#=================================" << myColor::RESET << endl << "#" << endl;

  log << myColor::RED << myColor::BOLD << "#Input parameters:" << myColor::RESET << endl <<  "#" << endl;

  options.write(log);
  log <<  myColor::GREEN << "#=================================" <<  myColor::RESET << endl;

  
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
  potential1.setBH();
  //potential1.set_WCA_limit(0.79);

  if(hsd1 < 0) hsd1 = potential1.getHSD(kT);

  if(cutoff < -1) cutoff = hsd1;
  
  //  do_N_first(Npoints_min, Npoints_max, Mu_min, Mu_max, Mu_step, GaussianOnly, potential1, log);
  do_Mu_first(Npoints_min, Npoints_max, Mu_min, Mu_max, Mu_step, GaussianOnly, potential1, log);

  log << "Finished" << endl;
  
  return 1;

}

void do_Mu_first(int Npoints_min, int Npoints_max, double Mu_min, double Mu_max, double Mu_step, bool GaussianOnly, Potential1& potential1, Log &log)
{
  for(int Npoints = Npoints_min; Npoints < Npoints_max; Npoints++)
    {
      double L[] = {Npoints*dx, Npoints*dx, Npoints*dx};
      SolidFCC theDensity1(dx, L);

      FMT_Species_Analytic species1(theDensity1,hsd1,1);
      Interaction_Linear_Interpolation i1(species1,species1,potential1,kT,log);

      //FMT_Species_Numeric species1(theDensity1,hsd1,pointsFile, 1);
      //Interaction i1(species1,species1,potential1,kT,log,pointsFile);
      //Interaction_Full i1(species1,species1,potential1,kT,log,Ngauss); 
      
      RSLT fmt;

      i1.initialize();
      log << "VDW parameter (potential)   = " << potential1.getVDW_Parameter(kT) << endl;      
      log << "VDW parameter (interaction) = " << 0.5*i1.getVDWParameter() << endl;
      log << "Diff in vdw parameter       = " << potential1.getVDW_Parameter(kT) - 0.5*i1.getVDWParameter() << endl << endl;

      theDensity1.setSpecies(&species1);

      DFT dft(&species1);
      dft.addHardCoreContribution(&fmt);  
      dft.addInteraction(&i1);

      log << "Npoints = " << Npoints << endl;

      double Fgauss = 0;
      double Agauss = -1;
      double Dgauss = 0;
      
      for(double Mu = Mu_min; Mu < Mu_max; Mu += Mu_step)	
	{
	  bool bstarted = false;
	  double f_prev = 0;
	  log << "Mu = " << Mu << " kT = " << kT << endl;
      
	  bool neg = (Mu < 0.0);	  
	  stringstream of_name;
	  of_name << "scan_mu_" << (neg ? '-' : '+') << std::fixed << std::setprecision(2)  << setfill('0') << setw(5) << fabs(Mu) << ".dat";

	  bool file_exists = false;
	  {
	    ifstream f(of_name.str().c_str());
	    file_exists = f.good();
	    f.close();
	  }
	  
	  
	  if(!file_exists) // initialize thermo and files
	    {
	      double xvap = -1;
	      try{
		xvap = dft.XVap_From_Mu(Mu,1.0);
		cout << "DFT find gas: " << xvap << endl;
	      } catch(...) { cout << "DFT find gas: none found" << endl;}
	      double xliq = -1;
	      try{
		xliq = dft.XLiq_From_Mu(Mu,1.0);
		cout << "DFT find liq: " << xliq << endl;
	      } catch(...) { cout << "DFT find liq: none found" << endl;}

	      vector<double> xv(1);
	      
	      ofstream of(of_name.str().c_str());
	      if(xliq > 0)
		{ xv[0] = xliq; of << "#Mu(input) = " << Mu << " Fliq/(kTV) = " << dft.Omega(xv) << " Dliq = " << xliq << " Mu(output) = " << dft.Mu(xv,0) << endl;}
	      if(xvap > 0)
		{ xv[0] = xvap; of << "#Mu(input) = " << Mu << " Fvap/(kTV) = " << dft.Omega(xv) << " Dvap = " << xvap << " Mu(output) = " << dft.Mu(xv,0) << endl;}

	      of << "#"
		 << setw(15) <<"Npoints"
		 << setw(15) <<"Fgau/(kTV)"
		 << setw(15) <<"Agau/(kTV)"		
		 << setw(15) <<"Density"
		 << setw(15) <<"F/(kT V)"
		 << setw(15) <<"Cvacancy"
		 << setw(15) <<"Err"
		 << endl;
	      of.close();
	    } else log << "Scan file \"" << of_name.str() << "\" exists: skipping thermodynamics" << endl;
	  
	  if(Agauss < 1)
	    {
	      // Do Gaussian: note that the minimum is independent of the chemical potential
	      // so we temporarily set it to zero to get the helmholtz part of the free energy
	      species1.setChemPotential(Mu); //1.0);
	      findGaussian(theDensity1, -20, dft, log, Fgauss, Agauss);
	      Dgauss = theDensity1.getNumberAtoms()/theDensity1.getVolume();
	    }
	  // set correct chemical potential
	  species1.setChemPotential(Mu);

	  double F = 0;
	  double Cvac = 0;
	  
	  if(findMinimum(theDensity1, dft, log, F, Cvac))
	    {
	      ofstream of(of_name.str().c_str(), ios::app);      		  
	      of << setw(15) << Npoints
		 << setw(15) << Fgauss - Mu*Dgauss
		 << setw(15) << Agauss	
		 << setw(15) << theDensity1.getNumberAtoms()/pow(Npoints*dx,3.0)
		 << setw(15) << F
		 << setw(15) << Cvac
		 << setw(15) << dft.get_convergence_monitor()
		 << endl;
	      of.close();
	      f_prev = F;
	    }
	}
    }
}


void do_N_first(int Npoints_min, int Npoints_max, double Mu_min, double Mu_max, double Mu_step, bool GaussianOnly, Potential1& potential1, Log &log)
{
  ofstream of("scan.dat");
  of.close();

  for(double Mu = Mu_max; Mu > Mu_min; Mu -= Mu_step)
    {
      bool bstarted = false;
      double f_prev = 0;
      log << "Mu = " << Mu << " kT = " << kT << endl;
      
      double bmu = Mu/kT;
      
      for(int Npoints = Npoints_min; Npoints < Npoints_max; Npoints++)
	{
	  double L[] = {Npoints*dx, Npoints*dx, Npoints*dx};
	  SolidFCC theDensity1(dx, L);

	  double F;
	  double Cvac;
	  
	  FMT_Species_Analytic species1(theDensity1,hsd1,Mu);
	  //FMT_Species_Numeric species1(theDensity1,hsd1,pointsFile, Mu);

	  //Interaction i1(species1,species1,potential1,kT,log, pointsFile);
	  //Interaction_Full i1(species1,species1,potential1,kT,log,2);
	  Interaction_Linear_Interpolation i1(species1,species1,potential1,kT,log);

	  RSLT fmt;

	  log << "VDW(potential)   = " << potential1.getVDW_Parameter(kT) << endl;
	  i1.initialize();
	  log << "VDW(interaction) = " << 0.5*i1.getVDWParameter() << endl;

	  theDensity1.setSpecies(&species1);

	  DFT dft(&species1);
	  dft.addHardCoreContribution(&fmt);  
	  dft.addInteraction(&i1);

	  log << "Npoints = " << Npoints << endl;
	  
	  if(Npoints == Npoints_min)
	    {
	      VDW1 vdw(hsd1,0.5*i1.getVDWParameter()); //potential1.getVDW_Parameter(kT));
	      
	      double xliq = vdw.findLiquidFromMu(Mu, 1.0);
	      double xvap = vdw.findVaporFromMu(Mu, 1.0);

	      ofstream of("scan.dat", ios::app);      		  
	      of << "#Mu = " << Mu << " Fliq/(kTV) = " << -vdw.pressure(xliq) << " Dliq = " << xliq << endl;
	      if(xvap > 0)
		of << "#Mu = " << Mu << " Fvap/(kTV) = " << -vdw.pressure(xvap) << " Dvap = " << xvap << endl;		
	      of << "#"
		 << setw(15) <<"Npoints"
		 << setw(15) <<"Fgau/(kTV)"
		 << setw(15) <<"Agau/(kTV)"		
		 << setw(15) <<"Density"
		 << setw(15) <<"F/(kT V)"
		 << setw(15) <<"Cvacancy"
		 << setw(15) <<"Err"
		 << endl;  
	      of.close();
	    }
	  double Fgauss, Agauss;
	  if(findGaussian(theDensity1, bmu, dft, log, Fgauss, Agauss))	  
	    if(GaussianOnly || findMinimum(theDensity1, dft, log, F, Cvac))
	      {
		bstarted = true; // found beginning of good range

		if(GaussianOnly) {F = Fgauss; Cvac = 0;}
				
		ofstream of("scan.dat", ios::app);      
		of << setw(15) << Npoints
		   << setw(15) << Fgauss		    
		   << setw(15) << Agauss	
		   << setw(15) << theDensity1.getNumberAtoms()/pow(Npoints*dx,3.0)
		   << setw(15) << F
		   << setw(15) << Cvac
		   << setw(15) << dft.get_convergence_monitor()
		   << endl;
		of.close();

		if(Npoints > Npoints_min &&  F > f_prev) break;
		f_prev = F;
	      }
	    else if(bstarted) break; // exhausted the good range
	}
      ofstream of1("scan.dat", ios::app);
      of1 << "#" << endl;
      of1.close();      
    }
}



bool doUniform(SolidFCC& theDensity1, DFT& dft, Log& log, double &Fret, double &Dret) 
{

  /////////////////////////////////////
  // Create objects

  theDensity1.initializeUniform(Dret);
  
  log <<  myColor::GREEN << "#=================================" <<  myColor::RESET << endl;
  
  myMin minimizer(dft,log);
  minimizer.setForceTerminationCriterion(forceLimit);
  minimizer.setTimeStep(dt);
  minimizer.setTimeStepMax(dtMax);
  minimizer.setAlphaStart(alpha_start);
  minimizer.setAlphaFac(alphaFac);
  minimizer.run(maxSteps); 

  double Natoms = theDensity1.getNumberAtoms();
  double Omega = minimizer.getF();
  //  double dOmega = Omega-L[0]*L[1]*L[2]*dft.Omega(densities);

  log <<  myColor::GREEN << "#=================================" << myColor::RESET << endl << "#" << endl;
  log << "Final Uniform Natoms               : " << Natoms << endl;
  log << "Final Omega                : " << Omega << endl;

  Dret = Natoms/theDensity1.getVolume();
  Fret = Omega/theDensity1.getVolume();
  return true;
}

double gaussianEval(double alf, double bmu, SolidFCC& theDensity, DFT& dft, double &prefac, Log&log)
{
  double f; 
  bool bSuccess = false;
  while(!bSuccess)
    {
      try {
	theDensity.initialize(alf, 1, prefac, exp(bmu));
	f = dft.calculateFreeEnergyAndDerivatives(false);
	bSuccess = true;
      } catch (...) {
	prefac *= 0.999999;
	log << "\t lowering prefac to " << prefac << endl;
      }
    }
  return f;
}

bool findGaussian(SolidFCC& theDensity1, double bmu, DFT& dft, Log& log, double &Fgau, double &Agau)
{
  /////////////////////////////////////////////////////
  // Report
  log << "Lx  = " << L[0] << endl;
  log << "HSD = " << hsd1 << endl;

  log << "# " << endl;
  log << "Determining gaussian approximation" << endl;
  log << "# " << endl;

  double f = 0;
  double f_old = 0;
  double f_old_2 = 0;  
  double alf_old = -1;
  double alf_old_2 = -1;
  double prefac_old = 0;
  bool isDecending = false;
  bool firstIteration = true;
  for(double alf = 30; alf < 10000; alf += 10)
    {
      f = gaussianEval(alf, bmu, theDensity1, dft, prefac,log);
      log << "alf = " << alf << " " << f << endl;
      stringstream title;
      title << "Gaussian: alf = " << alf << " f = " << f;
      string t = title.str();
      string file("dummy.dat");
      theDensity1.doDisplay(t, file, 0);

      
      if(!firstIteration)
	if(f < f_old) isDecending = true;
      
      //if(isDecending)
      //      	if(f > f_old) break;
      if(alf > 231) exit(0);

      f_old_2 = f_old;
      f_old = f;
      
      alf_old_2 = alf_old;
      alf_old = alf;

      prefac_old = prefac;

      firstIteration = false;
    }
  if(isDecending & alf_old_2 > 0)  // Min is between alf_old_2 and alf
    {
      double R = 0.61803399;
      double C = 1.0-R;

      double alf = alf_old + 10;
      
      double x0 = alf_old_2;
      double x3 = alf;
      double x1 = alf_old;
      double x2 = alf_old+C*(alf-alf_old);

      double f1 = gaussianEval(x1, bmu, theDensity1, dft, prefac,log);
      double f2 = gaussianEval(x2, bmu, theDensity1, dft, prefac,log);

      while(fabs(x3-x0) > 1e-3)
	{
	  if(f2 < f1) { x0 = x1; x1 = x2; x2 = R*x1+C*x3;
	    f1 = f2; f2 = gaussianEval(x2, bmu, theDensity1, dft, prefac,log);}
	  else { x3 = x2; x2 = x1; x1 = R*x2+C*x0;
	    f2 = f1; f1 = gaussianEval(x1, bmu, theDensity1, dft, prefac,log);}

	  log << "alf = " << x2 << " " << f2 << endl;
	  stringstream title;
	  title << "Gaussian: alf = " << x2 << " f = " << f2;
	  string t = title.str();
	  string file("dummy.dat");
	  theDensity1.doDisplay(t, file, 0);	  
	}
      if(f1 < f2) {Agau = x1; Fgau = f1;}
      else { Agau = x2; Fgau = f2;}

      theDensity1.initialize(Agau, 1, prefac); //_old);      
      log << "Found: alf = " << Agau << " prefac = " << prefac_old << " f = " << Fgau << endl;

      Fgau /= theDensity1.getVolume();
    } else {
    log << "No minimum found with Gaussians" << endl;
    return false;
  }

  return true;
}

bool findMinimum(SolidFCC& theDensity1, DFT& dft, Log& log, double &Fret, double &Cvac)
{

  log <<  myColor::GREEN << "#=================================" <<  myColor::RESET << endl;
  
  myMin minimizer(dft,log,1);
  minimizer.setForceTerminationCriterion(forceLimit);
  minimizer.setTimeStep(dt);
  minimizer.setTimeStepMax(dtMax);
  minimizer.setAlphaStart(alpha_start);
  minimizer.setAlphaFac(alphaFac);
  minimizer.run(maxSteps);

  double Natoms = theDensity1.getNumberAtoms();
  double Omega = minimizer.getF();
  //  double dOmega = Omega-L[0]*L[1]*L[2]*dft.Omega(densities);

  log <<  myColor::GREEN << "#=================================" << myColor::RESET << endl << "#" << endl;
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

/*
bool doCalc(SolidFCC& theDensity1, DFT& dft, Log& log, double &Fret, double &Cvac)
{
  /////////////////////////////////////////////////////
  // Report
  log << "Lx  = " << L[0] << endl;
  log << "HSD = " << hsd1 << endl;

  log << " " << endl;
  log << "Determining gaussian approximation" << endl;
  log << " " << endl;

  double f = 0;
  double f_old = 0;
  double f_old_2 = 0;  
  double alf_old = -1;
  double alf_old_2 = -1;
  double prefac_old = 0;
  bool isDecending = false;
  bool firstIteration = true;
  for(double alf = 30; alf < 10000; alf += 10)
    {
      f = gaussianEval(alf, bmu, theDensity1, dft, prefac,log);
      log << "alf = " << alf << " " << f << endl;

      if(!firstIteration)
	if(f < f_old) isDecending = true;
      
      if(isDecending)
	if(f > f_old) break;

      f_old_2 = f_old;
      f_old = f;
      
      alf_old_2 = alf_old;
      alf_old = alf;

      prefac_old = prefac;

      firstIteration = false;
    }
  if(isDecending & alf_old_2 > 0)  // Min is between alf_old_2 and alf
    {
      double R = 0.61803399;
      double C = 1.0-R;

      double alf = alf_old + 10;
      
      double x0 = alf_old_2;
      double x3 = alf;
      double x1 = alf_old;
      double x2 = alf_old+C*(alf-alf_old);

      double f1 = gaussianEval(x1, bmu, theDensity1, dft, prefac,log);
      double f2 = gaussianEval(x2, bmu, theDensity1, dft, prefac,log);

      while(fabs(x3-x0) > 1e-3)
	{
	  if(f2 < f1) { x0 = x1; x1 = x2; x2 = R*x1+C*x3;
	    f1 = f2; f2 = gaussianEval(x2, bmu, theDensity1, dft, prefac,log);}
	  else { x3 = x2; x2 = x1; x1 = R*x2+C*x0;
	    f2 = f1; f1 = gaussianEval(x1, bmu, theDensity1, dft, prefac,log);}
	}
      if(f1 < f2) {alf_old = x1; f_old = f1;}
      else { alf_old = x2; f_old = f2;}

      log << "Found: alf = " << alf_old << " prefac = " << prefac_old << " f = " << f_old << endl;
    } else {
    log << "No minimum found with Gaussians" << endl;
    return false;
  }

  theDensity1.initialize(alf_old, 1, prefac_old, exp(bmu));
  
  //check(theDensity1, dft,i1);

  log <<  myColor::GREEN << "#=================================" <<  myColor::RESET << endl;
  
  myMin minimizer(dft,log);
  minimizer.setForceTerminationCriterion(forceLimit);
  minimizer.setTimeStep(dt);
  minimizer.setTimeStepMax(dtMax);
  minimizer.setAlphaStart(alpha_start);
  minimizer.setAlphaFac(alphaFac);
  minimizer.run(maxSteps); 

  double Natoms = theDensity1.getNumberAtoms();
  double Omega = minimizer.getF();
  //  double dOmega = Omega-L[0]*L[1]*L[2]*dft.Omega(densities);

  log <<  myColor::GREEN << "#=================================" << myColor::RESET << endl << "#" << endl;
  log << "Final Natoms               : " << Natoms << endl;
  log << "Final Vacancy Concentration: " << (4.0-Natoms)/theDensity1.getVolume() << endl;  
  log << "Final Omega                : " << Omega << endl;

  Fret = Omega/theDensity1.getVolume();
  Cvac = (4.0-Natoms)/theDensity1.getVolume();
  return true;
}
*/
