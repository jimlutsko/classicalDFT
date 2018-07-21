
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

double range0 = 0;
double range1 = 2;

void Draw(Density& density, int image_number, mglGraph& gr, mglData &data_2D)
{
  density.fill_2D_data(data_2D);

  //clear the window
  gr.Clf();	

  // basic window formating
  gr.Box();  
  gr.Alpha(false);
  gr.SetRange('c', 0.0, 2.0);
    gr.SetRange('c', range0, range1);

  // set the color scheme
    gr.Dens(data_2D,"kRryw");
    //gr.Dens(data_2D,"UbcyqR");
    //gr.Dens(data_2D,"kBbcw");

  // Write a title
  char str[48];	
  snprintf(str,48,"Image = %d",image_number);
  gr.Puts(mglPoint(0,1.1),str);

  stringstream ss;
  std::ios  state(NULL);
  state.copyfmt(ss);
  
  ss <<  "copy_image_";
  ss << setfill ('0') << std::setw(8);
  ss << image_number;
  ss.copyfmt(state);
  ss <<".png";
  
  gr.WriteFrame(ss.str().c_str());

}




int main(int argc, char** argv)
{
  cout << "Build date  : " <<  (unsigned long) &__BUILD_DATE << endl;
  cout << "Build number: " << (unsigned long) & __BUILD_NUMBER << endl;

  cout << "OK" << endl;
  
  double L[3] = {10,10,10};
  int PointsPerHardSphere = 5;
  int nCores = 1;
  string infile;

  Options options;

  options.addOption("nCores", &nCores);
  options.addOption("PointsPerHardSphere", &PointsPerHardSphere);

  options.addOption("Lx", L);
  options.addOption("Ly", L+1);
  options.addOption("Lz", L+2);

  options.addOption("InputFile", &infile);
  
  options.read(argc, argv);

  cout << "OK" << endl;
  
  double dx = 1.0/PointsPerHardSphere;

#ifdef USE_OMP    
  //  omp_set_dynamic(0);
  //  omp_set_num_threads(nCores);

  //  int fftw_init_threads();
  //  fftw_plan_with_nthreads(omp_get_max_threads());
#endif

  double R = 1;
  double zPos = 0;
  Droplet finalDensity(dx, L, PointsPerHardSphere, R, zPos); 


  // Begin intialization of images

  //  finalDensity.initialize(xliq_coex,xgas_coex);



  //      string ss;
  //      
  //      cout << "Enter name of file: " << flush;
  //      cin >> ss;
  
  //      if(!ss.empty())
  //	finalDensity.readDensity(ss.c_str());
      //      else break;

  //      double NN = finalDensity.getNumberAtoms();

  //      cout << "NN = " << NN << endl;

      while(1)
	{
	  cout << "Enter range0: " << flush; cin >> range0;
	  cout << "Enter range1: " << flush; cin >> range1;
	  
	  mglGraph gr;
	  mglData data_2D;

	  //	  for(int i=49;i<50;i++)

	  ofstream of_radius("Radius.dat");
	  of_radius << "#Radius of gyration of the images" << endl;
	  
	  for(int j=0;j<50;j++)
	    {
	      Jglobal = j;
	      stringstream ss;
	      std::ios  state(NULL);
	      state.copyfmt(ss);

	      int i = j;
	      
	      ss << "archive_";
	      ss << setfill ('0') << std::setw(4);
	      ss << i;
	      ss.copyfmt(state);
	      ss <<".dat";  

	      string s = ss.str();

	      finalDensity.readDensity(ss.str().c_str());
	      finalDensity.initialize_2D_data(data_2D);	  	  
	      Draw(finalDensity,j,gr, data_2D);

	      double Rx = 0;
	      double Ry = 0;
	      double Rz = 0;
	      double R2 = 0;
	      double M = 0;

	      double m0 = finalDensity.getDensity(0,0,0);

	      for(int ix=0;ix<finalDensity.Nx(); ix++)
		for(int iy=0;iy<finalDensity.Ny(); iy++)
		  for(int iz=0;iz<finalDensity.Nz(); iz++)
		    {
		      double x = finalDensity.getX(ix);
		      double y = finalDensity.getX(iy);
		      double z = finalDensity.getX(iz);
		      double m = finalDensity.getDensity(ix,iy,iz);

		      m -= m0;
		      
		      Rx += m*x;
		      Ry += m*y;
		      Rz += m*z;

		      R2 += m*(x*x+y*y+z*z);		      
		      M += m;
		    }

	      Rx /= M;
	      Ry /= M;
	      Rz /= M;
	      double RR = Rx*Rx+Ry*Ry+Rz*Rz;
	      R2 = (R2/M)-RR;
	      of_radius << i << " " << sqrt(R2) << endl;
	    }
	  of_radius.close();
	}
  return 1;
}
