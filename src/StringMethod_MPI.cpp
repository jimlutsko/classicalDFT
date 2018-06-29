#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

#include <mgl2/mgl.h>
#include <mgl2/fltk.h>

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

  ofstream log(logfile.c_str());
  log << "#Initial free energies = ";

  //for(double F: oldF)
  ///    log << F << " ";
  //  log << endl;
  
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
  
  do {
    // Collect the densities
    MPI_Status *stat;
    int pos = 1;
    cout << "Reading densities ..." << endl;
    
    for(int Im=0;Im<taskList.size();Im++)
      for(int J=0;J<taskList[Im];J++)
	{
	  MPI_Recv(Images_[pos++].data(),Ntot_,MPI_DOUBLE,Im+1,/*tag*/ MPI_ANY_TAG ,MPI_COMM_WORLD,stat);
	  cout << " Read Im = " << Im << " J = " << J << endl;
	}

    
    // Do interpolation
    for(int k=0;k<5;k++)
      interpolate();

    // archive and draw images
    processImages();

    // send the densities back
    pos = 1;
    for(int Im=0;Im<taskList.size();Im++)
      for(int J=0;J<taskList[Im];J++)
	MPI_Send(Images_[pos++].data(),Ntot_,MPI_DOUBLE,Im+1,/*tag*/ 0 ,MPI_COMM_WORLD);

    // Report
    step_counter_++;    
    report(logfile);

    if(grace_) grace_->redraw(1,0);

  } while(1);

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
  
  for(int J=0;J<Images_.size();J++)
    {
      stringstream ss;
      std::ios  state(NULL);
      state.copyfmt(ss);
      
      ss << filename << "_";
      ss << setfill ('0') << std::setw(4);
      ss << J;
      ss.copyfmt(state);
      ss <<".dat";  

      DFT_Vec tmp;
      tmp.set(Images_[J].data(),finalDensity_.Ntot());

      ofstream of1(ss.str().c_str());
      tmp.save(of1);
      of1.close();
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

  ofstream status("status.dat");

  status << "0  0.000 " << dF_[0] << " 0 " << bav_*finalDensity_.getVolume() << endl; 
  
  MPI_Status *stat;

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

      MPI_Recv(dF1,  Nimage, MPI_DOUBLE, Im+1, 0, MPI_COMM_WORLD, stat);
      MPI_Recv(dist, Nimage, MPI_DOUBLE, Im+1, 0, MPI_COMM_WORLD, stat);
      MPI_Recv(N,    Nimage, MPI_DOUBLE, Im+1, 0, MPI_COMM_WORLD, stat);
      MPI_Recv(DT,   Nimage, MPI_DOUBLE, Im+1, 0, MPI_COMM_WORLD, stat);
            
      MPI_Recv(&dFmax1,     1, MPI_DOUBLE, Im+1, 0, MPI_COMM_WORLD, stat);
      MPI_Recv(&dFav1,      1, MPI_DOUBLE, Im+1, 0, MPI_COMM_WORLD, stat);
      MPI_Recv(&delta_max1, 1, MPI_DOUBLE, Im+1, 0, MPI_COMM_WORLD, stat);

      for(int K=0;K<Nimage;K++)
	{
	  status <<  Nimages+K << " " << DT[K] << " " << dF1[K] << " " << dist[K] << " " << N[K] << endl; 
	  dF_[1+Nimages+K] = dF1[K];
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
  status.close();      

  dFav /= Nimages;

  /////////////////// Report
  if(grace_)
    Display(1, dFmax, dFav);

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
  grace_->redraw(1,1);
  string s("string_graph.agr");
  grace_->store(s);
}



void StringMethod_MPI_Master::interpolate()
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
	if(alpha_local > alpha[K] && alpha_local < alpha[K+1])
	  {
	    Intervals.push_back(K);
	    found = true;
	  }
      if(!found)
	throw std::runtime_error("Could not locate interval for linear interpolation");
    }
  
  int chunk = finalDensity_.Ntot()/10;
  long i;

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
	  Images_[J][i] = x*y[K+1]+(1-x)*y[K];
	}
    }
}

void StringMethod_MPI_Slave::run(string& logfile)
{
  double dtmax = -1;
  int Nimages = string_.size();

  cout << "Task " << id_ << " beginning to relax" << endl;
  
  do {
    // Copy string for evaluating distance moved at end of step        
    for(int J=0;J<Nimages;J++)
      string_copy_[J].set(string_[J]->getDensity());

    ///////////////// Relaxation step

    double self_consistency_threshold = -1e-5;
    
    for(int J=0;J<Nimages;J++)
      {
	double dt = DT_[J];

	if(dtmax < 0) dtmax = dt;

	if(dt < dtmax) dt *= 2;

	ddft_.step_string(dt, *(string_[J]),self_consistency_threshold);
	
	cout << "ID " << id_ << " Image " << J << " dt = " << dt << endl;
		
	DT_[J] = dt;
      }
    cout << "Task " << id_ << " finished relaxtion" << endl;
    ///////////////// Interpolation step: send to master and wait to get back

    long Ntot = (Nimages > 0 ? string_[0]->Ntot() : 0);
    
    
    for(int J=0;J<Nimages;J++)
      MPI_Send(string_[J]->getData(), Ntot,MPI_DOUBLE,0,/*tag*/ 0 ,MPI_COMM_WORLD);

    cout << "Task " << id_ << " is waiting for interpolation ... " << endl;

    MPI_Status *stat;
    for(int J=0;J<Nimages;J++)
      MPI_Recv(string_[J]->getData(),Ntot,MPI_DOUBLE,0,/*tag*/ MPI_ANY_TAG ,MPI_COMM_WORLD,stat);
    
    /////////////////// Compute stats

    vector<double> newF(Nimages);
    dFmax_     = 0;
    dFav_      = 0;
    delta_max_ = 0;

    for(int J=0;J<Nimages;J++)
      {
	N_[J] = string_[J]->getNumberAtoms();

	double fmax = 0;
	newF[J] = ddft_.F_string(*(string_[J]), &fmax) - N_[J]*mu_;
	double ff = fabs(newF[J]-oldF_[J]);
	dFmax_ = max(ff,dFmax_);
	dFav_ += ff;

	distances_[J] = 0.0;
	for(long k=0;k<Nimages;k++)
	  {
	    double d = (string_copy_[J].get(k)-string_[J]->getDensity(k));
	    distances_[J]  += d*d;
	  }
	distances_[J]  = string_[J]->dV()*sqrt(distances_[J] )/DT_[J];

	delta_max_ = max(delta_max_, distances_[J] );

	
	cout << "ID " << id_ << " Image " << J << " F = " << newF[J]  << " N = " << N_[J] << " delta dist = " << distances_[J] << endl;
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
    
    //    string arc("archive");
    //    archive(arc);
    
  
  } while(1);


}
/*
void StringMethod_MPI_Slave::rescale_linear()
{
  // Calculate the "distance" of each image from the origin
  
  vector<double> alpha;
  alpha.push_back(0.0);

  for(int J = 1;J<string_.size();J++)
    {      
      double d = 0;
      for(long i=0;i<string_[0]->Ntot();i++)
	{
	  double d1 = string_[J-1]->getDensity(i);
	  double d2 = string_[J]->getDensity(i);

	  d += (d1-d2)*(d1-d2);
	}
      d = sqrt(d) + alpha.back();
      alpha.push_back(d);
    }


  // Interpolation.
  // First, find the interval we will use for each image
  
  double dL = alpha.back()/(string_.size()-1);

  vector<int> Intervals;
  Intervals.push_back(0); // first image is not interpolated
  for(int J = 1;J<string_.size()-1; J++) 
    {
      double alpha_local = J*dL;
      bool found = false;
      for(int K=0;K<alpha.size() && !found;K++)
	if(alpha_local > alpha[K] && alpha_local < alpha[K+1])
	  {
	    Intervals.push_back(K);
	    found = true;
	  }
      if(!found)
	throw std::runtime_error("Could not locate interval for linear interpolation");
    }

  //  cout << endl;
  //  cout << "Intervals for interpolation: " << endl;
  //  for(int J=0;J<string_.size()-1;J++)
  //    cout << J << " " << Intervals[J] << " - " << Intervals[J]+1 << endl;
  
  int chunk = string_[0]->Ntot()/10;
  long i;

#pragma omp parallel for			\
  shared(chunk)                                 \
  private(i)					\
  schedule(static,chunk)			  
  for(i=0;i<string_[0]->Ntot();i++)
    {
      vector<double> y;

      for(int J = 0;J<string_.size();J++)
	y.push_back(string_[J]->getDensity(i));

      for(int J = 1;J<string_.size()-1; J++) 
	{
	  double alpha_local = J*dL;
	  int K = Intervals[J];
	  double h = alpha[K+1]-alpha[K];
	  double x = (alpha_local-alpha[K])/h;
	  string_[J]->set_Density_Elem(i,x*y[K+1]+(1-x)*y[K]);
	  //	  string_[J]->set_Density_Elem(i,((alpha_local-alpha[K])*y[K+1]+(alpha[K+1]-alpha_local)*y[K])/(alpha[K+1]-alpha[K]));
	}
    }
}
*/


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




  //  ofstream dd("image_dump.dat");

  //  for(int i=0;i<density.Nx();i++)
  //    {
  //      for(int j=0;j<density.Ny();j++)
  //	dd << i << " " << j << " " << data_2D_.a[i+density.Nx()*j] << endl;
  //      cout << endl;
  //    }

  
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
      ss << J;
      ss.copyfmt(state);
      ss <<".dat";  

      string s = ss.str();
      string_[J]->readDensity(s.c_str());
    }

}
