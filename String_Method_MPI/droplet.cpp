
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

// required MPI include file  
#include <mpi.h>

#include "Grace.h"
#include "options.h"
#include "TimeStamp.h"


#include "Droplet.h"

#include "StringMethod_MPI.h"

#define MASTER 0        /* task ID of master task */

int main(int argc, char** argv)
{
  
  MPI_Status status;
  int	taskid,	        /* task ID - also used as seed number */
    numtasks,       /* number of tasks */
    rc;   /* return code */

  /* Obtain number of tasks and task ID */
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

  if(taskid == MASTER) cout << "MPI: numtasks = " << numtasks << endl;
  //  if(taskid == MASTER) cout << "MPI: taskid = " << taskid << endl;
  
  if(taskid == MASTER) cout << "Build date  : " <<  (unsigned long) &__BUILD_DATE << endl;
  if(taskid == MASTER) cout << "Build number: " << (unsigned long) & __BUILD_NUMBER << endl;


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

  double TimeStepMax = -1;
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
  options.addOption("TimeStepMax",&TimeStepMax);
  
  options.read(argc, argv, taskid == MASTER);

  if(taskid == MASTER)
    {
      ofstream log("log.dat");
      TimeStamp ts;
      log << "# " << ts << endl;
      log << "#=================================" << endl << "#" << endl;
      log << "#Input parameters:" << endl <<  "#" << endl;
      options.write(log);
      log << "#=================================" << endl;
      log.close();
    }

  if(TimeStepMax < 0) throw std::runtime_error("Must set TimeStepMax in input file ... aborting");
      
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
  } catch(...) { if(taskid == MASTER) cout << "No coexistence found ... " << endl;}
  
  double xs1,xs2;
  dft.spinodal(xs1,xs2);
  if(taskid == MASTER)
    {
      cout << "Spinodal: " << xs1 << " " << xs2 <<  " so Nspinodal = " << xs1*finalDensity.getVolume() << endl;
      cout << "Coexistence: " << xgas_coex << " " << xliq_coex << endl;
    }

  // set the thermodynamic state
  double omega_coex = dft.Omega(xliq_coex);
  double mu         = dft.Mu(xliq_coex);

  // Begin intialization of images

  finalDensity.initialize(xliq_coex,xgas_coex);

  double mu_boundary = 0;
  double bav = 0;
  if(taskid == MASTER)
    {
  
      if(! infile.empty())
	finalDensity.readDensity(infile.c_str());

      double NN = finalDensity.getNumberAtoms();

      cout << "Final NN = " << finalDensity.getNumberAtoms() << endl;
      cout << "Hard sphere diameter  = " << dft.HSD() << endl;
      cout << "Coexisting densities  = " << xliq_coex << " " << xgas_coex << endl;  
      cout << "Chemical potential/kT = " << mu << endl;
      cout << "beta * Grand free energy per unit volume = " << omega_coex << endl;
      
      double avDensity = NN/finalDensity.getVolume();
      
      cout << "Average density = " << avDensity << endl;

      // For closed system:
      double nav = 0;

      int Nx = finalDensity.Nx();
      int Ny = finalDensity.Ny();
      int Nz = finalDensity.Nz();

      for(int ix=0;ix<Nx;ix++)
	for(int iy=0;iy<Ny;iy++)
	  {
	    //finalDensity.set_Density_Elem(ix,iy,0,finalDensity.getDensity(ix,iy,0)*fac);
	    bav += finalDensity.getDensity(ix,iy,0);
	    nav++;
	  }
      for(int ix=0;ix<Nx;ix++)
	for(int iz=1;iz<Nz-1;iz++)
	  {
	    //	finalDensity.set_Density_Elem(ix,0,iz,finalDensity.getDensity(ix,0,iz)*fac);
	    //	finalDensity.set_Density_Elem(ix,Ny-1,iz,finalDensity.getDensity(ix,Ny-1,iz)*fac);
	    bav += finalDensity.getDensity(ix,0,iz);
	    nav++;
	  }

      for(int iy=1;iy<Ny-1;iy++)
	for(int iz=1;iz<Nz-1;iz++)
	  {
	    //	finalDensity.set_Density_Elem(0,iy,iz,finalDensity.getDensity(0,iy,iz)*fac);
	    //	finalDensity.set_Density_Elem(Nx-1,iy,iz,finalDensity.getDensity(Nx-1,iy,iz)*fac);
	    bav += finalDensity.getDensity(0,iy,iz);
	    nav++;	
	  }
      bav /= nav;

  
      mu_boundary = dft.Mu(bav);
      if(taskid == MASTER)
	{
	  ofstream log1("log.dat", ios::app);
	  cout << "Boundary Mu = " <<   mu_boundary << endl;
	  log1 << "#Boundary Mu = " <<   mu_boundary << endl;
	}

      // broadcast mu_boundary
      for(int jtask = 1; jtask < numtasks; jtask++)
	{
	  MPI_Send(&mu_boundary,1,MPI_DOUBLE,jtask,/*tag*/ 0 ,MPI_COMM_WORLD);
	  MPI_Send(&bav,1,MPI_DOUBLE,jtask,/*tag*/ 0 ,MPI_COMM_WORLD);
	}	  
    } else {
    MPI_Status *stat;
    MPI_Recv(&mu_boundary,1,MPI_DOUBLE,MPI_ANY_SOURCE,/*tag*/ MPI_ANY_TAG ,MPI_COMM_WORLD,stat );
    MPI_Recv(&bav,1,MPI_DOUBLE,MPI_ANY_SOURCE,/*tag*/ MPI_ANY_TAG ,MPI_COMM_WORLD,stat ); 
  }

  Grace *grace = NULL;

  if(taskid == MASTER)
    {
      if(showGraphics)
	grace = new Grace(800,600,2);

      if(!remove("string_graph.agr"))
	cout << "Removed string_graph.agr" << endl;
    }


  // There are numtasks-1 slave tasks available so each one will take some number of images
  // Since the first and last images are fixed, there are only N-2 to allocate
  

  dft.setEtaMax(1.0-1e-8);

  DDFT *ddft;
  if(ddft_type.compare("CLOSED") == 0) {ddft = new DDFT_IF(dft,finalDensity,NULL,showGraphics);;}
  else if(ddft_type.compare("OPEN_SIMPLE") == 0) {ddft = new DDFT_IF(dft,finalDensity,NULL,showGraphics); ddft->setFixedBoundary();}
  else if(ddft_type.compare("OPEN_SIMPLE_MODIFIED") == 0) {ddft = new DDFT_IF(dft,finalDensity,NULL,showGraphics); ddft->setFixedBoundary(); ddft->setModified();}
  else if(ddft_type.compare("OPEN_SIN") == 0) {ddft = new DDFT_IF_Open(dft,finalDensity,bav,NULL,showGraphics); }
  else {
    ofstream log1("log.dat", ios::app);
    cout << "DDFT type must be defined: current value is \"" << ddft_type << "\" which is unknown"  << endl;
    log1<< "DDFT type must be defined: current value is \"" << ddft_type << "\" which is unknown"  << endl;
  }

  
  //  DDFT_IF ddft(dft,finalDensity,NULL,showGraphics);
  ddft->initialize();

  ddft->setTimeStep(1e-2);
  ddft->set_tolerence_fixed_point(1e-4);
  ddft->set_max_time_step(1e-2);

  // For closed system:
  //  ddft.setFixedBoundary();

  long Ntot = finalDensity.Ntot();

  StringMethod_MPI *theString = NULL;
  
  if(taskid == MASTER)
    {
      double Finitial = dft.Fhelmholtz(bav)*finalDensity.getVolume();
      double Ni = bav*finalDensity.getVolume();
      
      cout << "Image 0 " << " F = " << Finitial << " mu = " << mu_boundary << " N = " << Ni  << " Omega = " << Finitial-mu_boundary*Ni << endl;

      double Ffinal = ddft->getF();
      double NN = finalDensity.getNumberAtoms();

      cout << "Image " << Nimages-1 << " F = " << Ffinal << " mu = " << mu_boundary << " N = " << NN  << " Omega = " << Ffinal-mu_boundary*NN << endl;
      
      theString = new StringMethod_MPI_Master(Nimages, finalDensity, bav, Ffinal-mu_boundary*NN, Finitial-Ni*mu_boundary, mu_boundary, grace, freeEnd);
      
      int assigned = 0;
      int chunk = (Nimages-2)/(numtasks-1);
      int left_over = (Nimages-2) - chunk*(numtasks-1);

      double *d = new double[Ntot];
      for(long i=0;i<Ntot;i++) d[i] = finalDensity.getDensity(i);
      
      for(int jtask = 1;jtask < numtasks; jtask++)
	{
	  int todo = chunk + (left_over ? 1 : 0);
	  if(left_over) left_over--;
	  
	  int start_index = assigned+1;
	  
	  MPI_Send(d,           Ntot,MPI_DOUBLE,jtask,/*tag*/ 0 ,MPI_COMM_WORLD);
	  MPI_Send(&start_index,   1,   MPI_INT,jtask,/*tag*/ 0 ,MPI_COMM_WORLD);
	  MPI_Send(&todo,          1,   MPI_INT,jtask,/*tag*/ 0 ,MPI_COMM_WORLD);

	  assigned += todo;
	  ((StringMethod_MPI_Master*) theString)->addTask(todo);	 
	}
      delete d;
    } else { 
    double *final = new double[Ntot];
    int todo;
    int start_index;
    double alf = 1.0/Nimages;
    MPI_Status *stat;
    
    MPI_Recv(final,      Ntot,MPI_DOUBLE,MPI_ANY_SOURCE,/*tag*/ MPI_ANY_TAG ,MPI_COMM_WORLD,stat );     
    MPI_Recv(&start_index,   1,MPI_INT   ,MPI_ANY_SOURCE,/*tag*/ MPI_ANY_TAG ,MPI_COMM_WORLD,stat );
    MPI_Recv(&todo,          1,MPI_INT   ,MPI_ANY_SOURCE,/*tag*/ MPI_ANY_TAG ,MPI_COMM_WORLD,stat ); 

    vector<Density*> Images(todo);
    for(int J=0;J<todo;J++)
      {
	Droplet *drop = new Droplet(dx, L, PointsPerHardSphere, R, zPos);
	double alf = (J+start_index)*1.0/(Nimages-1);	      

	for(long k=0;k<Ntot;k++)
	  drop->set_Density_Elem(k,(1-alf)*bav + alf*final[k]);
	Images[J] = drop;
      }
    delete final;

    theString = new StringMethod_MPI_Slave(*ddft, Images,mu_boundary, taskid, start_index, TimeStepMax);
    //    theString->setMu(mu_boundary);  

  }

  string s("log.dat");


  /*
    if(bRestart)
    {
    string file("..//archive");
    theString.read(file);
    }
  */

  //  theString->setMu(mu);
  
  theString->run(s);
  
  if(taskid == MASTER)
    {
      grace->pause();
      grace->close();
      delete grace;
    }

  delete ddft;
  
  MPI_Finalize();
  return 0;
}
