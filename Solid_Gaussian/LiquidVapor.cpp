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
#ifdef USE_GRACE
#include "Grace.h"
#include "LJ.h"
#endif


int main()
{
  double sigma = 1.0;
  double eps   = 1.0;
  double rcut  = 2.0;
  rcut = 1.2;
  
  //  LJ potential(sigma, eps, rcut);
  //  WHDF potential(sigma, eps, rcut);

  rcut = 2.5;
  tWF potential(sigma, eps, rcut);

#ifdef USE_GRACE
      Grace g;

      JZG lj(potential.getRcut(),true);
      
#endif
  
  double dT = 0.1;
  for(double kT = 0.2; kT < 2 && dT > 1e-8; kT += dT)
    {
      double hsd = potential.getHSD(kT);
      double avdw = potential.getVDW_Parameter(kT);

      double avdw1 = potential.getModifiedVDW_Parameter(kT,hsd);

      //      cout << "\t integral = " << potential.getModifiedVDW_Parameter(kT,hsd) << endl;
      //      cout << "\tavdw = " << avdw << endl;
      
      VDW1 vdw(hsd,avdw);
      VDW1 vdw1(hsd,avdw1);

#ifdef USE_GRACE

      stringstream ss;
      ss << "kT = " << kT;
      g.setTitle(ss.str().c_str());
      for(double x = 1e-8; x < 1.2; x += 0.01)
	{
	  g.addPoint(x,vdw.helmholtzPerUnitVolume(x),0);
	  g.addPoint(x,vdw1.helmholtzPerUnitVolume(x),1);
	  g.addPoint(x,x*lj.freeEnergy(x,kT),2);
	}      
      g.redraw();
      g.pause();
      g.deleteDataSet(0);
      g.deleteDataSet(1);
      g.deleteDataSet(2);
#endif




      
      
      double x1 = 1e-10;
      double x2 = 1;

      double xs1,xs2;

      try {
	vdw.spinodal(xs1,xs2);      
	vdw.findCoexistence(x1,x2);
	cout << kT << "\t" << xs1 << "\t" << xs2 << "\t" << x1 << "\t" << x2 << endl;	
      } catch(std::runtime_error &e) {
	//		cout << "\tkT = " << kT << " " << e.what() << endl;
	kT -= dT;
	dT/= 2;
      }      
    }
  
  return 1;
}
