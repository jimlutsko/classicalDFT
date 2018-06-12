#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

//#include <mgl2/mgl.h>
//#include <mgl2/fltk.h>

//#include "util/Grace.h"

using namespace std;


#include "Minimizer.h"

/*
void Picard::initialize()
{
  Minimizer::initialize();

  x_ = density_.getDensity().getData();

  F_ = dft_.calculateFreeEnergyAndDerivatives(density_,mu_, dF_);   

  cout << "Initial value of F = " << F_ << endl;

  calls_ = 0;
  step_counter_ = 0;
}


double Picard::step()
{
  static Grace g;

  double mix = Mixing_;

  vec ystep = x_;

  double dV = density_.dV();

  for(long i=0;i<density_.Ntot();i++)
    {
      double dd = x_[i];

      double dOmega_ex = dF_[i] - log(dd)*dV;

      ystep[i] =  exp(-dOmega_ex/dV);
    }


  bool OK = true;

  vec y = x_;
double Fold = F_;

  do {

for(long i=0;i<x_.size();i++)
  x_[i] = max(SMALL_VALUE, (1-mix)*y[i] + mix*ystep[i]);

    density_.set_density(x_);
    OK = true;


    try {
      F_ = dft_.calculateFreeEnergyAndDerivatives(density_,mu_, dF_);   
    } catch( Eta_Too_Large_Exception &e) {
      mix /= 2;
      OK = false;
     }
if(F_ - Fold > 0.01)
  {
cout << "\tmix = " << mix << " F_ = " << F_ << " Fold = " << Fold << endl;
mix /= 2;
 OK = false;
}

   } while(!OK);

  cout << "\tmix = " << mix << endl;

  g.deleteDataSet(0);
  g.deleteDataSet(1);

  int iy = density_.Ny()/2;
  int iz = density_.Nz()/2;

  for(int i=0;i<density_.Nx();i++)
    {
      long p = density_.pos(i,iy,iz);
      double dd = y(p);
      double d1 = x_(p);

      g.addPoint(i*density_.getDX(),dd,0);
      g.addPoint(i*density_.getDX(),d1,1);
    }
  g.redraw();
  //  g.pause();
  //  g.close();

  return F_;
}




*/
