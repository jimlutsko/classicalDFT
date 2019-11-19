#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>

#include <armadillo>

using namespace std;

#include "Potential1.h"
#include "VDW1.h"
#include "myColor.h"


int main()
{
  double sigma = 1.0;
  double eps   = 1.0;
  double rcut  = 3.0;

  LJ potential(sigma, eps, rcut);
      
  double dT = 0.1;
  for(double kT = 0.1; kT < 2 && dT > 1e-15; kT += dT)
    {
      double hsd = potential.getHSD(kT);
      double avdw = potential.getVDW_Parameter(kT);

      VDW1 vdw(hsd,avdw);
      double x1 = 1e-10;
      double x2 = 1;

      //      try {
	vdw.findCoexistence(x1,x2);
	cout << kT << "\t" << x1 << "\t" << x2 << endl;
	//    } catch (...) { kT -= dT; dT /= 10;}
    }
  
  return 1;
}
