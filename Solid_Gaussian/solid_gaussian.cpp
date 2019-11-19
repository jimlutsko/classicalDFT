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

#include "DFT.h"
#include "VDW1.h"
#include "Solid.h"
#include "Log.h"

bool findGaussian(SolidFCC& theDensity1, DFT& dft, Log& log, double Amax, double &Fgau, double &Agau, double &prefac);

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

  int nCores = 6;

  double kT = 1;
  double dx = 0.1;
  
  string pointsFile("..//SS31-Mar-2016//ss109.05998");
  
  Options options;

  options.addOption("Npoints_min", &Npoints_min);
  options.addOption("Npoints_max", &Npoints_max);
  options.addOption("Dx", &dx);
  
  options.addOption("nCores", &nCores);

  options.addOption("kT", &kT);
  options.addOption("Mu_min",   &Mu_min);
  options.addOption("Mu_max",   &Mu_max);
  options.addOption("Mu_step",  &Mu_step);
  
  options.addOption("eps1",   &eps1);
  options.addOption("sigma1", &sigma1);
  options.addOption("rcut1",  &rcut1);

  options.addOption("IntegrationPointsFile", &pointsFile);

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
  ////// Create potential
  LJ potential1(sigma1, eps1, rcut1);

  
  for(kT = 0.9; kT < 2.0; kT += 0.1)
    {      
      double hsd1 = potential1.getHSD(kT);

      
      // Part I: Starting with the smallest lattice parameter (number of cells),
      //  find F as the lattice parameter increases and stop when a minimum has been detected
      //  Store the values so that a minimum can be estimated below.

      vector<int> Np;
      vector<double> Natoms;
      vector<double> F;
      vector<double> alf;

      double Amax = 1000;
      double prefac = 1;
      
      for(int Npoints = Npoints_min; Npoints < Npoints_max; Npoints++)
	{
	  double L[] = {Npoints*dx, Npoints*dx, Npoints*dx};
	  SolidFCC theDensity1(dx, L);  

	  double Mu = 0.0;  	  
	  FMT_Species species1(theDensity1,hsd1,pointsFile,Mu,Npoints);
	  Interaction i1(species1,species1,potential1,kT,log, pointsFile);
	  RSLT fmt;

	  DFT dft(&species1);
	  dft.addHardCoreContribution(&fmt);  
	  dft.addInteraction(&i1);

	  log << "Npoints = " << Npoints << endl;

	  double Fgauss;
	  double Agauss;
	  if(findGaussian(theDensity1, dft, log, Amax, Fgauss, Agauss,prefac))	  
	    {
	      Np.push_back(Npoints);
	      F.push_back(Fgauss);
	      alf.push_back(Agauss);
	      Natoms.push_back(theDensity1.getNumberAtoms());

	      Amax = Agauss + 10;
	    } else break;
	}

      // Part II: For a range of chemical potentials, determine the
      // minimum for the solid (omega = f - mu n and both f and n depend on npoints).
      // Compare solid free energy to fluid and find the value of mu where they are equal (linear interpolation).
      
      VDW1 vdw(hsd1,potential1.getVDW_Parameter(kT));

      // Loop over phases
      for(int J=0;J<2;J++)
	{
	  // For linear interpolation of dOmega
	  double mua, mub, d_omega_a, d_omega_b, x_sol_a, x_sol_b;
	  d_omega_a = d_omega_b = 0;
  
	  for(double Mu = Mu_max; Mu > Mu_min; Mu -= Mu_step)
	    {
	      double xliq = (J == 0 ? vdw.findLiquidFromMu(Mu, 1.0) : vdw.findVaporFromMu(Mu, 1.0));
	      if(xliq < 0) continue;
	  
	      double omega_liquid = -vdw.pressure(xliq);

	      // Find solid
	      double oa, ob, oc;
	      double na, nb, nc;
	      double da, db, dc;
	      double aa, ab, ac;
	      da = db = dc = -1;

	      for(int i=0;i<Np.size();i++)
		{
		  double np = Np[i];
		  double V = (np*dx)*(np*dx)*(np*dx);

		  na = nb; nb = nc; nc = np;
		  oa = ob; ob = oc; oc = F[i] - Mu*Natoms[i]/V;
		  da = db; db = dc; dc = Natoms[i]/V;
		  aa = ab; ab = ac; ac = alf[i];
		  if(na > 0 && oc > ob) break;	  
		}
	      if(na > 0 && ob < oa && ob < oc)
		{
		  // Estimate minimum
		  double nmin = nb - 0.5*(nb-na)*(oa-oc)/(2*ob-oa-oc);
		  double omin = ((nmin-nb)*(nmin-nc)/((na-nb)*(na-nc)))*oa+((nmin-na)*(nmin-nc)/((nb-na)*(nb-nc)))*ob+((nmin-na)*(nmin-nb)/((nc-na)*(nc-nb)))*oc;
		  double dmin = ((nmin-nb)*(nmin-nc)/((na-nb)*(na-nc)))*da+((nmin-na)*(nmin-nc)/((nb-na)*(nb-nc)))*db+((nmin-na)*(nmin-nb)/((nc-na)*(nc-nb)))*dc;
		  double amin = ((nmin-nb)*(nmin-nc)/((na-nb)*(na-nc)))*aa+((nmin-na)*(nmin-nc)/((nb-na)*(nb-nc)))*ab+((nmin-na)*(nmin-nb)/((nc-na)*(nc-nb)))*ac;

		  d_omega_a    = d_omega_b;   d_omega_b    = omin-omega_liquid;
		  mua    = mub;   mub    = Mu;
		  x_sol_a  = x_sol_b; x_sol_b  = dmin;

		  cout << "Mu = " << Mu << " omega_sol = " << omin << " omega_liquid = " << omega_liquid << " df = " << omin-omega_liquid << endl;
		} else continue; // no solid identified ... 
	      if(d_omega_a*d_omega_b < 0)
		{
		  double Mu = mua - (mub-mua)*(d_omega_a/(d_omega_b-d_omega_a));      
		  double xsol = x_sol_a + (x_sol_b-x_sol_a)*((Mu-mua)/(mub-mua));      
		  double xliq = (J == 0 ? vdw.findLiquidFromMu(Mu, 1.0) : vdw.findVaporFromMu(Mu, 1.0));
		  cout << "Mucoex = " << Mu << " xliq = " << xliq << " xsol = " << xsol << endl;

		  ofstream of("coex.dat",ios::app);
		  of << kT << " ";	      
		  of << xliq << " " << xsol;
		  of << endl;
		  of.close();	      
		}
      	    
	    }
	}
    }
  return 1;
  
}

double gaussianEval(double alf, SolidFCC& theDensity, DFT& dft, double &prefac, Log&log)
{
  double f; 
  bool bSuccess = false;
  while(!bSuccess)
    {
      try {
	theDensity.initialize(alf, 1, prefac);
	f = dft.calculateFreeEnergyAndDerivatives(false);
	bSuccess = true;
      } catch (...) {
	prefac *= 0.99999;
	log << "\t lowering prefac to " << prefac << endl;
      }
    }
  return f;
}

bool findGaussian(SolidFCC& theDensity1, DFT& dft, Log& log, double Amax, double &Fgau, double &Agau, double &prefac)
{
  /////////////////////////////////////////////////////
  // Report
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
  for(double alf = 30; alf < Amax; alf += 10)
    {
      f = gaussianEval(alf, theDensity1, dft, prefac,log);
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

      double f1 = gaussianEval(x1, theDensity1, dft, prefac,log);
      double f2 = gaussianEval(x2, theDensity1, dft, prefac,log);

      while(fabs(x3-x0) > 1e-3)
	{
	  if(f2 < f1) { x0 = x1; x1 = x2; x2 = R*x1+C*x3;
	    f1 = f2; f2 = gaussianEval(x2, theDensity1, dft, prefac,log);}
	  else { x3 = x2; x2 = x1; x1 = R*x2+C*x0;
	    f2 = f1; f1 = gaussianEval(x1, theDensity1, dft, prefac,log);}
	}
      if(f1 < f2) {Agau = x1; Fgau = f1;}
      else { Agau = x2; Fgau = f2;}

      theDensity1.initialize(Agau, 1, prefac_old);      
      log << "Found: alf = " << Agau << " prefac = " << prefac_old << " f = " << Fgau << " Natoms = " << theDensity1.getNumberAtoms() << endl;

      Fgau /= theDensity1.getVolume();
    } else {
    log << "No minimum found with Gaussians" << endl;
    return false;
  }

  return true;
}

