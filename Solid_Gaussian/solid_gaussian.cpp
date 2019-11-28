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

#include <math.h>
#include <nlopt.h>

#include "options.h"

#include "DFT.h"
#include "VDW1.h"
#include "Solid.h"
#include "Log.h"

bool findGaussian(SolidFCC& theDensity1, DFT& dft, Log& log, double Amax, double &Fgau, double &Agau, double &prefac);
bool Picard(SolidFCC& theDensity1, DFT& dft, Log& log, double mu);
bool nlopt(SolidFCC& theDensity1, DFT& dft, Log& log1, double &F, double &alf, double &prefac, double &low);
  
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

  double kTmin = 1;
  double kTmax = 1;
  double dkT   = 1;
    
  double dx = 0.1;
  
  string pointsFile("..//SS31-Mar-2016//ss109.05998");
  
  Options options;

  options.addOption("Npoints_min", &Npoints_min);
  options.addOption("Npoints_max", &Npoints_max);
  options.addOption("Dx", &dx);
  
  options.addOption("nCores", &nCores);

  options.addOption("kTmin", &kTmin);
  options.addOption("kTmax", &kTmax);
  options.addOption("dkT", &dkT);
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
  //  LJ potential1(sigma1, eps1, rcut1);
  tWF potential1(sigma1, eps1, rcut1);

  
  for(double kT = kTmin; kT < kTmax+dkT/2; kT += dkT)
    {
      log << "kT = " << kT << endl;
      double hsd1 = potential1.getHSD(kT);
      
      // Part I: Starting with the smallest lattice parameter (number of cells),
      //  find F as the lattice parameter increases and stop when a minimum has been detected
      //  Store the values so that a minimum can be estimated below.

      vector<int> Np;
      vector<double> Natoms;
      vector<double> F;
      vector<double> alf;

      double Amax = 10000;

double Fgauss;
	  double Agauss = 700;
double Lgauss = 1e-4;
double prefac = 0.995;


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

	  log << " " << endl;
	  log << "Npoints = " << Npoints << endl;


//nlopt(theDensity1, dft, log);

//if(findGaussian(theDensity1, dft, log, Amax, Fgauss, Agauss,prefac))
if(nlopt(theDensity1, dft, log, Fgauss, Agauss,prefac, Lgauss))

	    {
	      Np.push_back(Npoints);
	      F.push_back(Fgauss);
	      alf.push_back(Agauss);
	      Natoms.push_back(theDensity1.getNumberAtoms());

//	      Picard(theDensity1, dft, log, -5);	      
	      Amax = Agauss + 10;
	    } else break;
	}

      // Part II: For a range of chemical potentials, determine the
      // minimum for the solid (omega = f - mu n and both f and n depend on the lattice parameter).
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
	      int na, nb, nc;
	      double da, db, dc;
	      double aa, ab, ac;
	      da = db = dc = -1;
	      na = nb = nc = -1;
	      
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

		  log << "Mu = " << Mu << " omega_sol = " << omin << " omega_liquid = " << omega_liquid << " df = " << omin-omega_liquid << " alf = " << amin << " density = " << dmin << " npts = " << nmin << endl;
		} else {log << "Mu = " << Mu << " no solid found: na = " << na << endl; continue; } // no solid identified ... 
	      if(d_omega_a*d_omega_b < 0)
		{
		  double Mu = mua - (mub-mua)*(d_omega_a/(d_omega_b-d_omega_a));      
		  double xsol = x_sol_a + (x_sol_b-x_sol_a)*((Mu-mua)/(mub-mua));      
		  double xliq = (J == 0 ? vdw.findLiquidFromMu(Mu, 1.0) : vdw.findVaporFromMu(Mu, 1.0));
		  log << "Mucoex = " << Mu << " xliq = " << xliq << " xsol = " << xsol << endl;

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

double x_alf;
SolidFCC *x_Density;
DFT *x_DFT;
Log *x_log;

void deriv(const double alpha, const double prefac, const double lower, const DFT_Vec &df, double *grad)
{
  double a_latt = x_Density->Lx();
    
  double atoms[4][3];
  
  atoms[0][0] = 0.0;
  atoms[0][1] = 0.0;
  atoms[0][2] = 0.0;
  
  atoms[1][0] = a_latt/2;
  atoms[1][1] = a_latt/2;
  atoms[1][2] = 0.0;
  
  atoms[2][0] = a_latt/2;
  atoms[2][1] = 0.0;
  atoms[2][2] = a_latt/2;
  
  atoms[3][0] = 0.0;
  atoms[3][1] = a_latt/2;
  atoms[3][2] = a_latt/2;
  
  double dV = x_Density->dV();

  int Nx = x_Density->Nx();
  int Ny = x_Density->Ny();
  int Nz = x_Density->Nz();

  grad[0] = grad[1] = grad[2] = 0.0;
  
  for(int i=0;i<Nx;i++)
    for(int j=0;j<Ny;j++)
      for(int k=0;k<Nz; k++)
	{
	  double x = x_Density->getX(i);
	  double y = x_Density->getY(j);
	  double z = x_Density->getZ(k);
	  
	  double dsum = 0;
	  double d0 = 0;
	  double d1 = 0;

	  long s = x_Density->pos(i,j,k);
	  double g = df.get(s);
	  
	  for(int l=0;l<4;l++)	     
	    {
	      double dx = fabs(x-atoms[l][0]); if(dx > a_latt/2) dx -= a_latt;
	      double dy = fabs(y-atoms[l][1]); if(dy > a_latt/2) dy -= a_latt;
	      double dz = fabs(z-atoms[l][2]); if(dz > a_latt/2) dz -= a_latt;
	      
	      double r2 = dx*dx+dy*dy+dz*dz;
	      double v = pow(alpha/M_PI,1.5)*exp(-alpha*r2);
	      dsum += v*prefac;
	      d0 += g*((1.5/alpha)-r2)*v*prefac;
	      d1 += g*v;
	    }
	  if(dsum < lower)
	    grad[2] += g;
	  else {
	    grad[0] += d0;
	    grad[1] += d1;
	  }
	}
}

double myfunc(unsigned n, const double *x, double *grad, void *data)
{
  double f = 0.0;

  try {

    x_Density->initialize(10000*x[0], 1, exp(x[1]), exp(-1000*x[2]));
    f = x_DFT->calculateFreeEnergyAndDerivatives(false);

    if(grad) {deriv(10000*x[0],exp(x[1]), exp(-1000*x[2]), x_DFT->getDF(0),grad);

    grad[0] *= 10000;
    grad[1] *=  exp(x[1]);
    grad[2] *= -1000*exp(-1000*x[2]);
    }
    
  } catch(...) { f = HUGE_VAL;}

  cout << "alf = " << 10000*x[0] << " pref = " << exp(x[1]) << " lower = " << exp(-1000*x[2]) << " f = " << f;
  if(grad)cout << " grad: " << grad[0] << " " << grad[1] << " " << grad[2] << endl;
  cout << endl;
  return f;
}




bool nlopt(SolidFCC& theDensity1, DFT& dft, Log& log1, double &F, double &alf, double &prefac, double &low)
{
  x_Density = &theDensity1;
  x_DFT = &dft;
  x_log = &log1;

  double lb[3] = {0.003,-HUGE_VAL, 0};
  double ub[3] = {HUGE_VAL,0.0, HUGE_VAL};
  
  //nlopt_opt opt = nlopt_create(NLOPT_LN_COBYLA,3);
  //       nlopt_opt opt = nlopt_create(NLOPT_LN_BOBYQA,3);
  //    nlopt_opt opt = nlopt_create(NLOPT_LN_PRAXIS,2);
  nlopt_opt opt = nlopt_create(NLOPT_LN_SBPLX,3);
  //    nlopt_opt opt = nlopt_create(NLOPT_LD_MMA,3);
    // nlopt_opt opt = nlopt_create(NLOPT_LD_LBFGS,3);
  //nlopt_opt opt = nlopt_create(NLOPT_LD_TNEWTON_PRECOND_RESTART,3);  

     
   nlopt_set_lower_bounds(opt, lb);
   nlopt_set_upper_bounds(opt, ub);   
   nlopt_set_min_objective(opt, myfunc, NULL);

   nlopt_set_xtol_rel(opt, 1e-4);
   nlopt_set_ftol_abs(opt,1e-4);
     
   double x[3] = { 0.08, 0.995, 1 };  /* `*`some` `initial` `guess`*` */

   x[0] = alf/10000;   
   x[1] = log(prefac); //log(0.995);
   x[2] = -log(low)/1000; //-log(1e-4)/1000;
   
   double minf; /* `*`the` `minimum` `objective` `value,` `upon` `return`*` */
   if (nlopt_optimize(opt, x, &minf) < 0)
     throw std::runtime_error("nlopt failed!");

   F = minf/theDensity1.getVolume();



   nlopt_destroy(opt);

   alf    = 10000*x[0];   
   prefac = exp(x[1]);
   low    = exp(-1000*x[2]);

   log1 << "Alf = " << alf << " prefac = " << prefac << " low = " << low << " f = " << F << endl;

   
   return true;
}
  


bool findGaussian(SolidFCC& theDensity1, DFT& dft, Log& log, double Amax, double &Fgau, double &Agau, double &prefac)
{
  /////////////////////////////////////////////////////
  // Report
  double f = 0;
  double f_old = 0;
  double f_old_2 = 0;  
  double alf_old = -1;
  double alf_old_2 = -1;
  double prefac_old = 0;
  bool isDecending = false;
  bool firstIteration = true;
  //  for(double alf = 30; alf < Amax; alf += 10)
  for(double alf = 2000; alf < Amax; alf += 10)
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

bool Picard(SolidFCC& theDensity1, DFT& dft, Log& log, double mu)
{
  double amix = 0.01;
  double dV = theDensity1.dV();
  
  do {
    dft.calculateFreeEnergyAndDerivatives(true);

    cout << "Err = " << dft.get_convergence_monitor() << endl;
    
    auto &f = dft.getDF(0);
    
    for(long i=0;i<theDensity1.Ntot();i++)
      {
	double d = theDensity1.getDensity(i);
	double F = f.get(i)/dV;
	d = (1-amix)*d+amix*exp(mu-F);

	theDensity1.set_Density_Elem(i, d);
      }
  } while(1);


}
