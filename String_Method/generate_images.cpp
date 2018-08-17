
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
double Omega;
double N_omega;

bool writeTitle = false;
bool addColorBar = false;
bool firstImageOnly = false;

void Draw(Density& density, int image_number, mglGraph& gr, mglData &data_2D)
{
  density.fill_2D_data(data_2D);

  //clear the window
  gr.Clf();	

  // basic window formating
  gr.Box();
  gr.Alpha(false);
  gr.SetRange('c', range0, range1);

  if(bLogarithmic)
    {
      gr.Dens(data_2D,"kRryw");
      gr.SetTicks('c',-5);
      if(addColorBar)
	if(!firstImageOnly || image_number == 0)
	  gr.Colorbar("<IkRryw");
    }  else {
    gr.Dens(data_2D,"kBbcw");
    if(addColorBar)
      if(!firstImageOnly || image_number == 0)
	gr.Colorbar("<IkBbcw");
  }
  // Write a title
  if(writeTitle)
    {
      char str[72];	
      snprintf(str,72,"Image = %d Omega = %.2lf N = %.2lf",image_number,Omega,N_omega);
      gr.Puts(mglPoint(0,1.1),str);
    }
  
  stringstream ss;
  std::ios  state(NULL);
  state.copyfmt(ss);

  if(bLogarithmic)
    ss <<  "log_image_";
  else ss <<  "image_";
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
  int Nimages = -1;
  
  Options options;

  options.addOption("PointsPerHardSphere", &PointsPerHardSphere);
  options.addOption("Nimages", &Nimages);
  options.addOption("Lx", L);
  options.addOption("Ly", L+1);
  options.addOption("Lz", L+2);
  options.addOption("WriteTitle", &writeTitle);
  options.addOption("ShowColorBar", &addColorBar);
  options.addOption("ShowColorBarFirstImage", &firstImageOnly);
  options.read(argc, argv);



  if(argc < 5)
    throw std::runtime_error("Usage: generate_images <inputfile> bLogarithmic(0 or 1) range_min range_max");

  int a = atoi(argv[2]);

  bLogarithmic = (a == 0 ? false : true);
  
  range0 = atof(argv[3]);
  range1 = atof(argv[4]);

  TimeStamp ts;
  
  ofstream ofl("IMAGE_LOG.txt");
  ofl << "# " << ts << endl;
  ofl << "Input file: " << argv[1] << endl;
  ofl << "bLogarithmic = " << bLogarithmic << endl;
  ofl << "range min = " << range0 << endl;
  ofl << "range max = " << range1 << endl;
  ofl.close();
  
  double dx = 1.0/PointsPerHardSphere;

  double R = 1;
  double zPos = 0;
  Droplet finalDensity(dx, L, PointsPerHardSphere, R, zPos); 

  mglGraph gr(0,400,400);
  mglData data_2D;

  ofstream of_radius("Radius.dat");
  of_radius << "#Radius of gyration of the images" << endl;

  ifstream in_state("status.dat");
  
  for(int j=0;j<Nimages;j++)
    {
      string status;
      getline(in_state,status);
      stringstream ss_status(status);
      double a,b,c;
      ss_status >> a >> b >> Omega >> c >> N_omega;
      
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

      if(j == 0) continue;
	      
      double Rx = 0;
      double Ry = 0;
      double Rz = 0;
      double R2 = 0;
      double M = 0;
      double m_max = 0.0;
      double m0 = finalDensity.getDensity(0,0,0);
      double m2 = 0;
	      
      for(int ix=0;ix<finalDensity.Nx(); ix++)
	for(int iy=0;iy<finalDensity.Ny(); iy++)
	  for(int iz=0;iz<finalDensity.Nz(); iz++)
	    {
	      double x = finalDensity.getX(ix);
	      double y = finalDensity.getX(iy);
	      double z = finalDensity.getX(iz);
	      double m = finalDensity.getDensity(ix,iy,iz);

	      m -= m0;

	      m_max = max(m_max,m);
	      m2 += m*m;
		      
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

      R2 *= 5/3.0;

      double d = 0.0;
      long n = 0;
      for(int ix=0;ix<finalDensity.Nx(); ix++)
	for(int iy=0;iy<finalDensity.Ny(); iy++)
	  for(int iz=0;iz<finalDensity.Nz(); iz++)
	    {
	      double x = finalDensity.getX(ix);
	      double y = finalDensity.getX(iy);
	      double z = finalDensity.getX(iz);
	      double m = finalDensity.getDensity(ix,iy,iz);

	      double rr = (x*x+y*y+z*z);		      

	      if(rr <= R2){d += m; n++;}
		      
	    }

      of_radius << i << " " << sqrt(R2) << " " << d/n << " " << M*finalDensity.dV()/(4*M_PI*pow(R2,1.5)/3) << endl;
    }
  of_radius.close();

  return 1;
}
