#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

#include <mgl2/mgl.h>

#ifdef USE_OMP
#include <omp.h>
#endif

using namespace std;

// required MPI include file  
#include <mpi.h>

#include "StringMethod_MPI.h"
#include "spliner.h"


void StringMethod_MPI_Master::run(string& logfile)
{

  ofstream log(logfile.c_str(), ios::app);
  log << "#Initial free energies = ";

  //for(double F: oldF)
  ///    log << F << " ";
  log << endl;
  
  log << "#step_counter\tdFmax\tdFav\tdist_max" << endl;
  log.close();

  //Eliminate existing files, if any ...
  int ret = system("rm image_*.png");

  int Nimages = Images_.size();
  
  for(int J=0;J<Nimages;J++)
    Images_[J].resize(Ntot_);

  for(long i=0;i<Ntot_;i++)
    {
      Images_[0][i] = bav_;
      Images_[Nimages-1][i] = finalDensity_.getDensity(i);
    }

  report(logfile);
  
  do {
    // Collect the densities
    int pos = 1;
    cout << "Waiting for densities ..." << endl;
    
    for(int Im=0;Im<taskList.size();Im++)
      for(int J=0;J<taskList[Im];J++)
	MPI_Recv(Images_[pos++].data(),Ntot_,MPI_DOUBLE,Im+1,/*tag*/ MPI_ANY_TAG ,MPI_COMM_WORLD,MPI_STATUS_IGNORE);


    cout << "... begin interpolation" << endl;
    
    // Do interpolation
    double movement = 1;
    int k = 0;
    do {
	movement = interpolate();
	cout << "Interpolation iteration " << ++k << " has movement = " << movement << endl;
    } while(movement > interpolation_tolerence_); //1e-6);

    // archive and draw images
    cout << "post-process images ..." << endl;
    processImages();

    
    // send the densities back
    cout << "Send densities back ..." << endl;    
    pos = 1;
    for(int Im=0;Im<taskList.size();Im++)
      for(int J=0;J<taskList[Im];J++)
	MPI_Send(Images_[pos++].data(),Ntot_,MPI_DOUBLE,Im+1,/*tag*/ 0 ,MPI_COMM_WORLD);

    // Report
    step_counter_++;    
    report(logfile);

    if(grace_) grace_->redraw(1,0);

    cout << "delta_max_ = " << delta_max_ << " termination_criterion_ = " << termination_criterion_ << endl;

  } while(delta_max_ > termination_criterion_);

}

void StringMethod_MPI_Master::processImages()
{
  int N = Images_.size();

  if(grace_)
    {

      int Nx = finalDensity_.Nx();
      int Ny = finalDensity_.Ny();
      int Nz = finalDensity_.Nz();
  
      for(int J=0;J<N;J++)	
	grace_->deleteDataSet(J,1);
      for(int J=0;J<N;J++)
	{
	  for(int i=(Nx-1)/2;i<Nx;i++)
	    grace_->addPoint(finalDensity_.getX(i), Images_[J][finalDensity_.pos(i,(Ny-1)/2,(Nz-1)/2)],J,1);
	}
    }

  // archive
  string arc("archive");
  archive(arc);    

  for(int J=0;J<N;J++)
    Draw(Images_[J], J, dF_[J]);
}


void StringMethod_MPI_Master::archive(string &filename) const
{
  stringstream ss;
  ss << filename << "_summary.dat";
  ofstream of(ss.str().c_str());
  of << finalDensity_.Nx() << " " << finalDensity_.Ny()  << " " << finalDensity_.Nz()  << endl;
  of << finalDensity_.Lx() << " " << finalDensity_.Ly()  << " " << finalDensity_.Lz()  << endl;
  of << Images_.size() << endl;
  of.close();    


  stringstream ss1;
  std::ios  state(NULL);
  state.copyfmt(ss1);
  
  int J = 0;

  ss1 << filename << "_";
  ss1 << setfill ('0') << std::setw(4);
  ss1 << J;
  ss1.copyfmt(state);
  ss1 <<".dat";  

  DFT_Vec dd;
  dd.set(Images_[J].data(), finalDensity_.Ntot());
  finalDensity_.set(dd);
  string s1 = ss1.str();
  finalDensity_.writeDensity(s1);

  if(!freeEnd_)
    {
      stringstream ss2;
      std::ios  state2(NULL);
      state2.copyfmt(ss2);
      J = Images_.size()-1;

      ss2 << filename << "_";
      ss2 << setfill ('0') << std::setw(4);
      ss2 << J;
      ss2.copyfmt(state2);
      ss2 <<".dat";  

      dd.set(Images_[J].data(), finalDensity_.Ntot());
      finalDensity_.set(dd);
      string s2 = ss2.str();
      finalDensity_.writeDensity(s2);
    }
}


void StringMethod_MPI_Master::report(string &logfile)
{
  // get the stats from the slave processes
  // and process the images as the info comes in


  double dFmax = 0;
  double dFav = 0;
  double delta_max = 0;
  int Nimages = 0;

  //  ofstream status("status.dat");

  stringstream status;
  
  status << "0  0.000 " << dF_[0] << " 0 " << N_[0] << endl; 
  
  for(int Im=0;Im<taskList.size();Im++)
    {
      int Nimage = taskList[Im];

      double dFmax1 = 0;
      double dFav1 = 0;
      double delta_max1 = 0;
               
      double *dF1 = new double[Nimage];			      
      double *dist = new double[Nimage];
      double *N = new double[Nimage];
      double *DT = new double[Nimage];

      MPI_Recv(dF1,  Nimage, MPI_DOUBLE, Im+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(dist, Nimage, MPI_DOUBLE, Im+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(N,    Nimage, MPI_DOUBLE, Im+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(DT,   Nimage, MPI_DOUBLE, Im+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
      MPI_Recv(&dFmax1,     1, MPI_DOUBLE, Im+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&dFav1,      1, MPI_DOUBLE, Im+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&delta_max1, 1, MPI_DOUBLE, Im+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      for(int K=0;K<Nimage;K++)
	{
	  status <<  1+Nimages+K << " " << DT[K] << " " << dF1[K] << " " << dist[K] << " " << N[K] << endl; 
	  dF_[1+Nimages+K] = dF1[K];
	  N_[1+Nimages+K] = N[K];
	}
      Nimages  += Nimage;
      dFmax     = max(dFmax, dFmax1);
      delta_max = max(delta_max, delta_max1);
      dFav     += dFav1*Nimage;

      delete dF1;
      delete dist;
      delete N;
      delete DT;
    }
  if(!freeEnd_)
    {
      int nlast = Images_.size()-1;
      status <<  1+Nimages << " " << "0.0" << " " << dF_[nlast] << " " << "0.0" << " " << N_[nlast] << endl;
    }

  ofstream status_file("status.dat");
  status_file << status.str() << endl;
  status_file.close();      

  dFav /= Nimages;

  /////////////////// Report
  if(grace_)
    Display((step_counter_ == 0 ? 0 : 1), dFmax, dFav);

  cout << "Step = " << step_counter_ << " dFmax = " << dFmax << " dFav = " << dFav << endl;
  cout << endl;
  
  ofstream log(logfile.c_str(), ios::app);
  log.precision(12);
  log << step_counter_ << "\t"
      << dFmax << "\t"
      << dFav << "\t"
      << delta_max << "\t";
  log << endl;
  log.close();

  delta_max_ = delta_max;

  
}

void StringMethod_MPI_Master::Display(int dataSet, double dFmax, double dFav)
{

  int N = dF_.size();
  
  grace_->deleteDataSet(dataSet,0);
  for(int J=0;J<N;J++)
    grace_->addPoint(J,dF_[J],dataSet,0);
  grace_->circle(1,0);
  grace_->redraw(1,0);

  stringstream ss;
  ss << "Step = " << step_counter_ << " dFmax = " << dFmax << " dFav = " << dFav;
  cout << "Setting title to: " << ss.str() << endl;
  grace_->setTitle(ss.str().c_str());
  
  grace_->redraw(1,0);
  if(step_counter_ > 0) grace_->redraw(1,1);
  string s("string_graph.agr");
  grace_->store(s);
}



double StringMethod_MPI_Master::interpolate()
{
  // Calculate the "distance" of each image from the origin

  int Nimages = Images_.size();
  
  vector<double> alpha;
  alpha.push_back(0.0);

  for(int J = 1;J<Nimages;J++)
    {      
      double d = 0;
      for(long i=0;i<Ntot_;i++)
	{
	  double d1 = Images_[J-1][i];
	  double d2 = Images_[J][i];

	  d += (d1-d2)*(d1-d2);
	}
      d = sqrt(d) + alpha.back();
      alpha.push_back(d);
    }


  // Interpolation.
  // First, find the interval we will use for each image
  
  double dL = alpha.back()/(Nimages-1);

  vector<int> Intervals;
  Intervals.push_back(0); // first image is not interpolated
  for(int J = 1;J<Nimages-1; J++) 
    {
      double alpha_local = J*dL;
      bool found = false;
      for(int K=0;K<alpha.size() && !found;K++)
	if(alpha_local >= alpha[K] && alpha_local < alpha[K+1])
	  {
	    Intervals.push_back(K);
	    found = true;
	  }
      if(!found)
	throw std::runtime_error("Could not locate interval for linear interpolation");
    }
  
  int chunk = finalDensity_.Ntot()/10;
  long i;

  double movement = 0.0;
  
#pragma omp parallel for			\
  shared(chunk)                                 \
  private(i)					\
  schedule(static,chunk)			  
  for(i=0;i<Ntot_;i++)
    {
      vector<double> y;

      for(int J = 0;J<Nimages;J++)
	y.push_back(Images_[J][i]);

      for(int J = 1;J<Nimages-1; J++) 
	{
	  double alpha_local = J*dL;
	  int K = Intervals[J];
	  double h = alpha[K+1]-alpha[K];
	  double x = (alpha_local-alpha[K])/h;

	  double dx = Images_[J][i];
	  double sx = dx;
	  Images_[J][i] = x*y[K+1]+(1-x)*y[K];
	  dx -= Images_[J][i];
	  sx += Images_[J][i];
	  dx /= (0.5*sx);
	  movement += dx*dx;
	}
    }
  movement /= Ntot_;
  movement /= (Nimages-2);
  movement = sqrt(movement);
  return movement;
}


void StringMethod_MPI_Master::interpolate_cubic()
{
  // Calculate the "distance" of each image from the origin

  int Nimage = Images_.size();
  
  vector<double> alpha;
  alpha.push_back(0.0);

  for(int J = 1;J<Nimage;J++)
    {      
      double d = 0;
      for(long i=0;i<Ntot_;i++)
	{
	  double d1 = Images_[J-1][i];
	  double d2 = Images_[J][i];

	  d += (d1-d2)*(d1-d2);
	}
      d = sqrt(d) + alpha.back();
      alpha.push_back(d);
    }

  
  // Now, interpolate to get the new images

  double dL = alpha.back()/(Nimage-1);

  int chunk = Ntot_/10;
  long i;
  
  //#pragma omp parallel for			\
  //  shared(chunk)				\
  //  private(i)				\
  //  schedule(static,chunk)			  
  for(i=0;i<Ntot_;i++)
    {
      vector<double> y;

      for(int J = 0;J<Nimage;J++)
	y.push_back(sqrt(Images_[J][i]));

      SplinerVec s(alpha,y);

      for(int J = 1;J<Nimage-1; J++)
	{
	  double d = s.f(J*dL);
	  d = d*d;	  
	  Images_[J][i] = d;
	}
    }
}


void StringMethod_MPI_Slave::run(string& logfile)
{
  int Nimages = string_.size();
  
  for(int J=0;J<Nimages;J++)
      string_copy_[J].set(string_[J]->getDensity());

  report();
   
  cout << "Task " << id_ << " beginning to relax" << endl;
  
  do {
    // Copy string for evaluating distance moved at end of step        
    for(int J=0;J<Nimages;J++)
      string_copy_[J].set(string_[J]->getDensity());

    ///////////////// Relaxation step

    for(int J=0;J<Nimages;J++)
      {
	double dt = DT_[J];

	dt = min(2*dt, Time_Step_Max_);
	ddft_.step_string(dt, *(string_[J]),false);
		
	DT_[J] = dt;
      }
    cout << "Task " << id_ << " finished relaxtion" << endl;
    ///////////////// Interpolation step: send to master and wait to get back

    long Ntot = (Nimages > 0 ? string_[0]->Ntot() : 0);
    
    
    for(int J=0;J<Nimages;J++)
      MPI_Send(string_[J]->getData(), Ntot,MPI_DOUBLE,0,/*tag*/ 0 ,MPI_COMM_WORLD);

    for(int J=0;J<Nimages;J++)
      MPI_Recv(string_[J]->getData(),Ntot,MPI_DOUBLE,0,/*tag*/ MPI_ANY_TAG ,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    
    /////////////////// Compute stats and send to master
    report();
    string arc("archive");
    archive(arc);        
  } while(1);
}

void StringMethod_MPI_Slave::report()
{
    int Nimages = string_.size();
    long Ntot = (Nimages > 0 ? string_[0]->Ntot() : 0);
  
    vector<double> newF(Nimages);
    dFmax_     = 0;
    dFav_      = 0;
    delta_max_ = 0;

    for(int J=0;J<Nimages;J++)
      {
	N_[J] = string_[J]->getNumberAtoms();

	double fmax = 0;
	double fj = ddft_.F_string(*(string_[J]), &fmax);
	newF[J] = fj - N_[J]*mu_;
	double ff = fabs(newF[J]-oldF_[J]);
	dFmax_ = max(ff,dFmax_);
	dFav_ += ff;

	distances_[J] = 0.0;
	double mass = 0.0;
	for(long k=0;k<Ntot;k++)
	  {
	    double d1 = string_copy_[J].get(k);
	    double d2 = string_[J]->getDensity(k);
	    
	    //	    double d = 2*(d1-d2)/(d1+d2);
	    double d = (d1-d2);
	    distances_[J]  += d*d;
	    mass += d1;
	    //	    distances_[J]  = max(d*d, distances_[J]);

	  }
	//	distances_[J]  = string_[J]->dV()*sqrt(distances_[J] )/DT_[J];
	//	distances_[J]  = sqrt(distances_[J] )/(DT_[J]*Ntot);
	distances_[J]  = sqrt(distances_[J] )/(DT_[J]*mass);

	delta_max_ = max(delta_max_, distances_[J] );
	cout << "Image " << J+offset_ << " F = " << fj << " mu = " << mu_ << " N = " << N_[J]  << " Omega = " << newF[J] << endl;
      }

    dFav_ /= Nimages;
    oldF_ = newF;

    MPI_Send(newF.data(),       Nimages,MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    MPI_Send(distances_.data(), Nimages,MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    MPI_Send(N_.data(),         Nimages,MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    MPI_Send(DT_.data(),        Nimages,MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    
    MPI_Send(&dFmax_,     1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&dFav_,      1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&delta_max_, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

}

void StringMethod_MPI_Master::Draw(vector<double> &data, int image_number, double F)
{
  if(!gr_) return;

  // A quick fix: this should be changed back
  //  density.fill_2D_data(data_2D_);

  int Nx = finalDensity_.Nx();
  int Ny = finalDensity_.Ny();
  int Nz = finalDensity_.Nz();
  
  for(int i=0;i<Nx;i++)
    for(int j=0;j<Ny;j++)
      data_2D_.a[i+Nx*j] = data[finalDensity_.pos(i,j, int((Nz-1)/2))]; 	
  
  //clear the window
  gr_->Clf();	

  // basic window formating
  gr_->Box();  
  gr_->Alpha(false);
  gr_->SetRange('c', 0.0, 2.0);

  // set the color scheme
  //	  gr_->Dens(a,"kRryw");
  //	  gr_->Dens(a,"UbcyqR");
  gr_->Dens(data_2D_,"kBbcw");

  // Write a title
  char str[48];	
  snprintf(str,48,"Step = %ld Image = %d F = %lf",step_counter_, image_number, F);
  gr_->Puts(mglPoint(0,1.1),str);

  stringstream ss;
  std::ios  state(NULL);
  state.copyfmt(ss);
  
  ss << "image_";
  ss << setfill ('0') << std::setw(8);
  ss << image_number;
  ss.copyfmt(state);
  ss <<".png";
  
  gr_->WriteFrame(ss.str().c_str());

}


void StringMethod_MPI_Slave::read(string &filename)
{
  for(int J=0;J<string_.size();J++)
    {
      stringstream ss;
      std::ios  state(NULL);
      state.copyfmt(ss);
      
      ss << filename << "_";
      ss << setfill ('0') << std::setw(4);
      ss << J+offset_;
      ss.copyfmt(state);
      ss <<".dat";  

      string s = ss.str();
      string_[J]->readDensity(s.c_str());
    }

}

void StringMethod_MPI_Slave::archive(string &filename) const
{
  for(int J=0;J<string_.size();J++)
    {
      stringstream ss;
      std::ios  state(NULL);
      state.copyfmt(ss);
      
      ss << filename << "_";
      ss << setfill ('0') << std::setw(4);
      ss << J+offset_;
      ss.copyfmt(state);
      ss <<".dat";  

      string s = ss.str();
      string_[J]->writeDensity(s);
    }
}
