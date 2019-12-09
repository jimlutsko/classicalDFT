#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>

#include <armadillo>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>


using namespace std;

#include "EOS.h"
#include "LJ.h"
#include "myColor.h"
#ifdef USE_GRACE
#include "Grace.h"
#endif


int main()
{
  double kT = 1;
  double density = 0.5;

  JZG eos;
  eos.findCriticalPoint(density,kT);
  cout << "density = " << density << " kT = " << kT << endl;

  Grace g;

  double dT = 0.001;

  double x1 = 0.8;
  double x2 = 0.001;

  double xs1 = 0.7;
  double xs2 = 0.01;
  
  for(double kT1 = 0.6; kT1 < kT; kT1 += dT)
    {
      if(eos.findCoexistence(x1,x2,kT1))
	{
	  g.addPoint(max(x1,x2),kT1,0);
	  g.addPoint(min(x1,x2),kT1,1);
	  g.redraw();
	}

      if(eos.findSpinodal(xs1,xs2,kT1))
	{
	  g.addPoint(max(xs1,xs2),kT1,2);
	  g.addPoint(min(xs1,xs2),kT1,3);
	  g.redraw();
	}      
    }
  g.addPoint(density,kT,0);
  g.addPoint(density,kT,1);

  g.addPoint(density,kT,2);
  g.addPoint(density,kT,3);  
    
  g.redraw();
  g.pause();

  g.close();  
  
  return 1;
}

const EOS *eos;

int criticalPoint(const gsl_vector * v, void *params, gsl_vector * f)
{
  const double x = gsl_vector_get (v, 0);
  const double kT = gsl_vector_get (v, 1);

  const double y0 = eos->dP_dx(x,kT);
  const double y1 = eos->dP_dx2(x,kT);

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);

  return GSL_SUCCESS;
}

int print_state (size_t iter, gsl_multiroot_fsolver * s)
{
  printf ("iter = %3u x = % .3f % .3f "
          "f(x) = % .3e % .3e\n",
          int(iter),
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1));
}

// Solve for dP/dn = 0, d2P/dn2 = 0;
// Now, P = kT*(n + kT phi1) so
// dP/dn = kT*(1 + phi2)
// d2P/dn2 = kT*phi3
bool EOS::findCriticalPoint(double &density, double &kT) const
{
  eos = this;
  
  const size_t n = 2;
  gsl_multiroot_function f = {&criticalPoint, n, NULL};

  gsl_vector *x = gsl_vector_alloc (n);
  gsl_vector_set (x, 0, density);
  gsl_vector_set (x, 1, kT);

  const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
  gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc (T, 2);
  gsl_multiroot_fsolver_set (s, &f, x);

  int status;
  size_t iter = 0;  
  //  print_state (iter, s);

  do {
    iter++;
    status = gsl_multiroot_fsolver_iterate (s);
    
    //    print_state (iter, s);
    
    if (status)   /* check if solver is stuck */
      break;
    
    status =
      gsl_multiroot_test_residual (s->f, 1e-7);
  }
  while (status == GSL_CONTINUE && iter < 1000);

  printf ("status = %s %d\n", gsl_strerror(status),status);

  density = gsl_vector_get(s->x, 0);
  kT      = gsl_vector_get(s->x, 1);
  
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);

  return (status == 0);
}

// Look for equal pressures and equal chemical potentials:
// mu = df/dx = ln(x) + d(x phix)/dx
int coexistence(const gsl_vector * v, void *params, gsl_vector * f)
{
  const double x1 = gsl_vector_get (v, 0);
  const double x2 = gsl_vector_get (v, 1);

  double kT = *((double *) params);
  
  const double y0 = eos->Pressure(x1,kT) - eos->Pressure(x2,kT);
  const double y1 = (kT*log(x1) + eos->phix(x1,kT)+x1*eos->dphix_dx(x1,kT) - kT*log(x2) - eos->phix(x2,kT)-x2*eos->dphix_dx(x2,kT));

  gsl_vector_set (f, 0, y0/(x1-x2));
  gsl_vector_set (f, 1, y1/(x1-x2));

  return GSL_SUCCESS;
}

bool EOS::findCoexistence(double &x1, double &x2, double &kT) const
{
  eos = this;
  
  const size_t n = 2;
  gsl_multiroot_function f = {&coexistence, n, &kT};

  gsl_vector *x = gsl_vector_alloc (n);
  gsl_vector_set (x, 0, x1);
  gsl_vector_set (x, 1, x2);

  const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
  gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc (T, 2);
  gsl_multiroot_fsolver_set (s, &f, x);

  int status;
  size_t iter = 0;  
  //  print_state (iter, s);

  do {
    iter++;
    status = gsl_multiroot_fsolver_iterate (s);
    
    //    print_state (iter, s);
    
    if (status)   /* check if solver is stuck */
      break;    
    status = gsl_multiroot_test_residual (s->f, 1e-7);
  }
  while (status == GSL_CONTINUE && iter < 1000);

  //  printf ("status = %s %d\n", gsl_strerror(status),status);

  //  print_state (iter, s);
  
  x1 = gsl_vector_get(s->x, 0);
  x2 = gsl_vector_get(s->x, 1);
  
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);

  return (status == 0);
}


// Look for dP/dx = 0
int spinodal(const gsl_vector * v, void *params, gsl_vector * f)
{
  const double x1 = gsl_vector_get (v, 0);
  const double x2 = gsl_vector_get (v, 1);

  double kT = *((double *) params);
  
  const double y0 = eos->dP_dx(x1,kT);
  const double y1 = eos->dP_dx(x2,kT);

  gsl_vector_set (f, 0, y0/(x1-x2));
  gsl_vector_set (f, 1, y1/(x1-x2));

  return GSL_SUCCESS;
}

bool EOS::findSpinodal(double &x1, double &x2, double &kT) const
{
  eos = this;
  
  const size_t n = 2;
  gsl_multiroot_function f = {&spinodal, n, &kT};

  gsl_vector *x = gsl_vector_alloc (n);
  gsl_vector_set (x, 0, x1);
  gsl_vector_set (x, 1, x2);

  const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
  gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc(T, 2);
  gsl_multiroot_fsolver_set(s, &f, x);

  int status = 1;
  for(size_t iter = 0;  iter < 1000; iter++)
    {
      if((status = gsl_multiroot_fsolver_iterate(s)))   // check if solver is stuck 
	break;    
      if((status = gsl_multiroot_test_residual(s->f, 1e-7)) != GSL_CONTINUE)   // check if solver is finished
	break;          
    }
  x1 = gsl_vector_get(s->x, 0);
  x2 = gsl_vector_get(s->x, 1);
  
  gsl_multiroot_fsolver_free(s);
  gsl_vector_free(x);

  return (status == 0);
}
