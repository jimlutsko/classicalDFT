
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

int main(int argc, char** argv)
{
  cout << "Build date  : " <<  (unsigned long) &__BUILD_DATE << endl;
  cout << "Build number: " << (unsigned long) & __BUILD_NUMBER << endl;


  double L[3] = {10,10,10};
  int PointsPerHardSphere = 5;
  string outfile("dump.dat");
  string infile;
  double Natoms = -1;
  double alpha = 50;
  double baseline = 3e-4;
  
  Options options;

  options.addOption("PointsPerHardSphere", &PointsPerHardSphere);
  options.addOption("Lx", L);
  options.addOption("Ly", L+1);
  options.addOption("Lz", L+2);

  options.addOption("OutputFile", &outfile);
  options.addOption("InputFile", &infile);
  options.addOption("Natoms", &Natoms);
  options.addOption("alpha", &alpha);
  options.addOption("Baseline", &baseline);


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


  ////// Create potential

  Droplet theDensity(dx, L, PointsPerHardSphere, infile,alpha,baseline);

  theDensity.initialize();
  
  string out("snapshot.dat");
  theDensity.writeDensity(out);






  

  return 1;
}
