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

#include "Minimizer.h"

void DDFT::initialize()
{
  Minimizer::initialize();
  successes_ = 0;
}

//void DDFT::reverseForce(DFT_Vec *tangent) 
//{
//  double prod = dF_.dotWith(*tangent);
//  dF_.Increment_And_Scale(*tangent,-2*prod);
//}

void DDFT::Display(double F, double dFmin, double dFmax, double N)
{
  /*
  static int cc = 0;
  
  grace_->deleteDataSet(0);
  grace_->deleteDataSet(1);
  grace_->deleteDataSet(2);
  for(int i=0;i<density_.Nx();i++)
    grace_->addPoint(i,density_.getDensity(i,density_.Ny()/2, density_.Nz()/2),0);
  
  for(int i=0;i<density_.Ny();i++)
    grace_->addPoint(i,density_.getDensity(density_.Nx()/2,i, density_.Nz()/2),1);
  
  for(int i=0;i<density_.Nz();i++)
    grace_->addPoint(i,density_.getDensity(density_.Nx()/2,density_.Ny()/2,i),2);
  
  if(cc == 0)
    for(int i=0;i<density_.Nz();i++)
      grace_->addPoint(i,density_.getDensity(density_.Nx()/2,density_.Ny()/2, i),3);
  
  cc = 1;
  
  grace_->redraw();
  
  stringstream ss;
  ss << "Step = " << step_counter_ << " F = " << F << " dFMin = " << dFmin << " dFmax = " << dFmax << " N = " << N;
  cout << "Setting title to: " << ss.str() << endl;
  grace_->setTitle(ss.str().c_str());
  
  grace_->redraw(1,0);
  string s("string_graph.agr");
  grace_->store(s);
  */
}

  /*
double DDFT::F_string(Density &original_density, double *fmax)
{
  double F = 0;

  DFT_Vec dummy;
  double F =  dft_.calculateFreeEnergyAndDerivatives(original_density,0.0, dummy,false);

  if(fmax) *fmax = dummy.inf_norm()/original_density.dV();;
  return F;
}
  */


