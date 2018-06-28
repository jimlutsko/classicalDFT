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

#include "StringMethod_MPI.h"
#include "spliner.h"

void StringMethod_MPI_Master::Display(vector<double> &F, int dataSet, double dFmax, double dFav)
{

  int N = F.size();
  
  grace_->deleteDataSet(dataSet,0);
  for(int J=0;J<N;J++)
    grace_->addPoint(J,F[J],dataSet,0);
  grace_->circle(1,0);
  grace_->redraw(1,0);

  /*
  int Nx = string_[0]->Nx();
  int Ny = string_[0]->Ny();
  int Nz = string_[0]->Nz();
  
  for(int J=0;J<N;J++)	
    grace_->deleteDataSet(J,1);
  for(int J=0;J<N;J++)
    {
      for(int i=(Nx-1)/2;i<Nx;i++)
	grace_->addPoint(string_[J]->getX(i), string_[J]->getDensity(i,(Ny-1)/2,(Nz-1)/2),J,1);
    }
  */
  stringstream ss;
  ss << "Step = " << step_counter_ << " dFmax = " << dFmax << " dFav = " << dFav;
  cout << "Setting title to: " << ss.str() << endl;
  grace_->setTitle(ss.str().c_str());
  
  grace_->redraw(1,0);
  grace_->redraw(1,1);
  string s("string_graph.agr");
  grace_->store(s);
}


void StringMethod_MPI_Master::run(string& logfile)
{
  for(int i=0;i<Images_.size();i++)
    Images_[i].resize(Ntot_);

  do {
    // Collect the densities
    MPI_Status *stat;
    int pos = 1;
    for(int I=0;I<taskList.size();I++)
      {
	for(int J=0;J<taskList[I];J++)
	  {
	    MPI_Recv(d,Ntot,MPI_DOUBLE,I,/*tag*/ MPI_ANY_TAG ,MPI_COMM_WORLD,stat);
	    Images_[pos] = d;
	  }
      }


  } while(1);



  
  if(grace_) grace_->redraw(1,0);

  ofstream log(logfile.c_str());
  log << "#Initial free energies = ";

  //for(double F: oldF)
  ///    log << F << " ";
  //  log << endl;
  
  log << "#step_counter\tdFmax\tdFav\tdist_max" << endl;
  log.close();

  // Draw images
  int ret = system("rm image_*.png");

  //for(int J=0;J<N;J++)
  //    Draw(*(string_[J]),J, oldF[J]);

    step_counter_++;



}

void StringMethod_MPI_Slave::run(string& logfile)
{
  double dtmax = -1;

  do {
    ///////////////// Execute one step
    

    // Copy string
        
    for(int J=1;J<string_.size()-1;J++)
      string_copy_[J].set(string_[J]->getDensity());
    
    cout << "Relaxing now ... " << endl;
    int Limit = (freeEnd_ ? string_.size() : string_.size()-1);

    double self_consistency_threshold = -1e-5;
    
    for(int J=1;J<Limit;J++)
      {
	double dt = DT_[J];

	if(dtmax < 0) dtmax = dt;

	if(dt < dtmax) dt *= 2;

	ddft_.step_string(dt, *(string_[J]),self_consistency_threshold);
	
	cout << "Image " << J << " dt = " << dt << endl;
		
	DT_[J] = dt;
      }

    cout << "Rescaling ... " << endl;
    for(int k=0;k<5;k++)
      rescale_linear();


    /////////////////// Compute stats

    int N = string_.size();
    
    vector<double> newF(N);
    vector<double> distances(N);
    vector<double> Natoms(N);
    double dFmax = 0;
    double dFav  = 0;
    double fmax_climber = 0;
    double delta_max = 0;
    cout << "string_.size()-1 = " << string_.size()-1 << endl;
    for(int J=0;J<N;J++)
      {
	double fmax = 0;
	newF[J] = ddft_.F_string(*(string_[J]), &fmax) - string_[J]->getNumberAtoms()*mu_;
	double ff = fabs(newF[J]-oldF_[J]);
	dFmax = max(ff,dFmax);
	dFav += ff;

	distances[J] = 0.0;
	for(long k=0;k<string_copy_[J].size();k++)
	  {
	    double d = (string_copy_[J].get(k)-string_[J]->getDensity(k));
	    distances[J]  += d*d;
	  }
	distances[J]  = string_[J]->dV()*sqrt(distances[J] )/DT_[J];

	delta_max = max(delta_max, distances[J] );

	Natoms[J] = string_[J]->getNumberAtoms();
	
	cout << "Image " << J << " F = " << newF[J] - mu_*Natoms[J] << " N = " << Natoms[J] << " delta dist = " << distances[J] << endl;
      }

    dFav /= N;

    /////////////////// Report
    //    if(grace_)
    //      Display(newF,1, dFmax, dFav);

    //    cout << "Step = " << step_counter_ << " dFmax = " << dFmax << " dFav = " << dFav << " fmax_climber = " << fmax_climber << endl;
    //    cout << endl;
    /*
    ofstream log(logfile.c_str(), ios::app);
    log.precision(12);
    log << step_counter_ << "\t"
	<< dFmax << "\t"
	<< dFav << "\t"
	<< delta_max << "\t";
    log << endl;
    log.close();
    
    log.close();

    ofstream status("status.dat");
    for(int J=0;J<N;J++)
      status << J << " " << DT_[J] << " " << newF[J] << " " << distances[J] << " " << Natoms[J] << endl; 
    status.close();
    */

    oldF_ = newF;

    
    string arc("archive");
    archive(arc);    
  
  } while(1);


}

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


void StringMethod_MPI_Slave::Draw(Density& density, int image_number, double F)
{
  if(!gr_) return;

  density.fill_2D_data(data_2D_);



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

void StringMethod_MPI_Slave::archive(string &filename) const
{
  stringstream ss;
  ss << filename << "_summary.dat";
  ofstream of(ss.str().c_str());
  of << string_[0]->Nx() << " " << string_[0]->Ny()  << " " << string_[0]->Nz()  << endl;
  of << string_[0]->Lx() << " " << string_[0]->Ly()  << " " << string_[0]->Lz()  << endl;
  of << string_.size() << endl;
  of.close();
  
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
      string_[J]->writeDensity(s);
      of.close();
    }
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
