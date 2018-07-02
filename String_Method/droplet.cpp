
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


#include "Droplet.h"

#include "StringMethod.h"



int main(int argc, char** argv)
{
  cout << "Build date  : " <<  (unsigned long) &__BUILD_DATE << endl;
  cout << "Build number: " << (unsigned long) & __BUILD_NUMBER << endl;


  double L[3] = {10,10,10};
  int PointsPerHardSphere = 5;
  int nCores = 1;

  double R = -1;
  double zPos = 0;

  double density = 0.8;

  double alpha = 0.01;
  int MaxIts = 100;

  double kT = 1;

  double eps = 1;
  double sigma = 1;
  double rcut = 1;

  double sigma_conj_grad = 1;

  string pointsFile("..//SS31-Mar-2016//ss109.05998");
  string outfile("dump.dat");
  string infile;
  double Natoms = -1;

  bool showGraphics = true;
  double forceLimit = -1;

  int Nimages = 10;
  bool bRestart = false;
  int Jclimber = -1;
  bool freeEnd = false;

  string ddft_type;
  
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

  options.addOption("sigma_conj_grad", &sigma_conj_grad);
  options.addOption("ForceTerminationCriterion",&forceLimit);

  options.addOption("OutputFile", &outfile);
  options.addOption("IntegrationPointsFile", &pointsFile);
  options.addOption("InputFile", &infile);
  options.addOption("Natoms", &Natoms);
  options.addOption("ShowGraphics", &showGraphics);

  options.addOption("Nimages",&Nimages);
  options.addOption("Restart",&bRestart);
  options.addOption("Jclimber",&Jclimber);
  options.addOption("FreeEnd",&freeEnd);

  options.addOption("DDFT_Type",&ddft_type);
  
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

  Droplet finalDensity(dx, L, PointsPerHardSphere, R, zPos); 
  Potential1 potential(sigma, eps, rcut);
  DFT_VDW<RSLT> dft(finalDensity,potential,pointsFile,kT);

  // Determine coexistence to give us a baseline for the densities
  double xliq_coex = 0.001;
  double xgas_coex = 0.6;
  try{
    dft.coexistence(xliq_coex, xgas_coex);
  } catch(...) { cout << "No coexistence found ... " << endl;}
  
  double xs1,xs2;
  dft.spinodal(xs1,xs2);
  cout << "Spinodal: " << xs1 << " " << xs2 <<  " so Nspinodal = " << xs1*finalDensity.getVolume() << endl;
  cout << "Coexistence: " << xgas_coex << " " << xliq_coex << endl;

  // set the thermodynamic state
  double omega_coex = dft.Omega(xliq_coex);
  double mu         = dft.Mu(xliq_coex);

  // Begin intialization of images

  finalDensity.initialize(xliq_coex,xgas_coex);

  if(! infile.empty())
    finalDensity.readDensity(infile.c_str());

  double NN = finalDensity.getNumberAtoms();
  if(Natoms > 0) NN = Natoms;

  cout << "Final NN = " << finalDensity.getNumberAtoms() << endl;

  cout << "Hard sphere diameter  = " << dft.HSD() << endl;
  cout << "Coexisting densities  = " << xliq_coex << " " << xgas_coex << endl;  
  cout << "Chemical potential/kT = " << mu << endl;
  cout << "beta * Grand free energy per unit volume = " << omega_coex << endl;
  string s("log.dat");


  double avDensity = NN/finalDensity.getVolume();

  cout << "Average density = " << avDensity << endl;

  // For closed system:
  double bav = 0;
  double nav = 0;

  int Nx = finalDensity.Nx();
  int Ny = finalDensity.Ny();
  int Nz = finalDensity.Nz();

  double bmax = 0;
  double bmin = 100;
  
  for(int ix=0;ix<Nx;ix++)
    for(int iy=0;iy<Ny;iy++)
      {
	//finalDensity.set_Density_Elem(ix,iy,0,finalDensity.getDensity(ix,iy,0)*fac);
	bav += finalDensity.getDensity(ix,iy,0);
	bmax = max(bmax,finalDensity.getDensity(ix,iy,0));
	bmin = min(bmin,finalDensity.getDensity(ix,iy,0));
	nav++;
      }
  for(int ix=0;ix<Nx;ix++)
    for(int iz=1;iz<Nz-1;iz++)
      {
	//	finalDensity.set_Density_Elem(ix,0,iz,finalDensity.getDensity(ix,0,iz)*fac);
	//	finalDensity.set_Density_Elem(ix,Ny-1,iz,finalDensity.getDensity(ix,Ny-1,iz)*fac);
	bav += finalDensity.getDensity(ix,0,iz);
	bmax = max(bmax,finalDensity.getDensity(ix,0,iz));
	bmin = min(bmin,finalDensity.getDensity(ix,0,iz));
	nav++;
      }

  for(int iy=1;iy<Ny-1;iy++)
    for(int iz=1;iz<Nz-1;iz++)
      {
	//	finalDensity.set_Density_Elem(0,iy,iz,finalDensity.getDensity(0,iy,iz)*fac);
	//	finalDensity.set_Density_Elem(Nx-1,iy,iz,finalDensity.getDensity(Nx-1,iy,iz)*fac);
	bav += finalDensity.getDensity(0,iy,iz);
	bmax = max(bmax,finalDensity.getDensity(0,iy,iz));
	bmin = min(bmin,finalDensity.getDensity(0,iy,iz));
	nav++;	
      }
  bav /= nav;

  ofstream log1("log.dat");
  
  cout << "bmax = " << bmax << endl;
  cout << "bmin = " << bmin << endl;
  cout << "bav  = " << bav  << endl;
  
  double mu_boundary = dft.Mu(bav);
  cout << "Boundary Mu = " <<   mu_boundary << endl;
  log1 << "#Boundary Mu = " <<   mu_boundary << endl;

  vector<Density*> Images(Nimages); //,finalDensity);
  double Rmax = 3.0;

  for(int J=0;J<Nimages;J++)
    {
      Droplet *d = new Droplet(dx, L, PointsPerHardSphere, R, zPos);

      if(J == 0)
	{
	  // For closed system:
	  d->initialize(bav,bav,&finalDensity);
	  // For open system:
	  // d->initialize(avDensity,avDensity);
	} else d->initializeTo(*((Droplet*) Images.front()),finalDensity,1-J*1.0/(Images.size()-1));

      Images[J] = d;
    }

  for(int i=0;i<Images.size();i++)
    {
      cout << "Image " << i << " has N = " << Images[i]->getNumberAtoms() << " atoms" << endl;
      log1 << "#Image " << i << " has N = " << Images[i]->getNumberAtoms() << " atoms" << endl;
    }

  dft.setEtaMax(1.0-1e-8);


  DDFT *ddft;
  if(ddft_type.compare("CLOSED") == 0) {ddft = new DDFT_IF(dft,finalDensity,NULL,showGraphics);;}
  else if(ddft_type.compare("OPEN_SIMPLE") == 0) {ddft = new DDFT_IF(dft,finalDensity,NULL,showGraphics); ddft->setFixedBoundary();}
  else if(ddft_type.compare("OPEN_SIMPLE_MODIFIED") == 0) {ddft = new DDFT_IF(dft,finalDensity,NULL,showGraphics); ddft->setFixedBoundary(); ddft->setModified();}
  else if(ddft_type.compare("OPEN_SIN") == 0) {ddft = new DDFT_IF_Open(dft,finalDensity,bav,NULL,showGraphics); }
  else {
    cout << "DDFT type must be defined: current value is \"" << ddft_type << "\" which is unknown"  << endl;
    log1<< "DDFT type must be defined: current value is \"" << ddft_type << "\" which is unknown"  << endl;
  }

  log1.close();    
  
  //DDFT_IF_Open ddft(dft,finalDensity,bav,NULL,showGraphics);  
  //  DDFT_IF ddft(dft,finalDensity,NULL,showGraphics);
  // For closed system:
  //  ddft.setFixedBoundary();
  
  ddft->initialize();

  ddft->setTimeStep(1e-2);
  ddft->set_tolerence_fixed_point(1e-4);
  ddft->set_max_time_step(1e-2);

  Grace *grace = NULL;
  if(showGraphics)
    grace = new Grace(800,600,2);

  StringMethod theString(*ddft, Images, grace, freeEnd);

  if(Jclimber > 0) theString.setClimbingImage(Jclimber);

  theString.setMu(mu_boundary);
  
  if(remove("string_graph.agr"))
    cout << "Error removing string_graph.agr" << endl;
  else cout << "Removed string_graph.agr" << endl;

  
  if(bRestart)
    {
      string file("..//archive");
      theString.read(file);
    }
  
  theString.run(s);
  grace->pause();
  grace->close();
  delete grace;
  delete ddft;

  return 1;
}
