


#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <time.h>
#include <random>

#ifdef USE_OMP
#include <omp.h>
#endif

using namespace std;

#include "Dynamical_Matrix.h"



// This is the method of Boutsidis, Drineas, Kambadur, Kontopoulou and Zouzias
// https://www.boutsidis.org/Boutsidis_LAA2017.pdf
// The idea is quite simple: Let the eigenvalues of A be l_i
// write log(detA) as = log(prod_i(l_i) = sum_i(log(1-(1-l_i)))
//                    = -sum_i(sum_k(l_i^k)/k) =-sum_k(1/k)sum_i(l_i^k)
//                    = -sum_k((1/k)Tr(A^k))
// and the trace is evaluated stochastically.

double Dynamical_Matrix::log_det_1(double lam_max, int num_samples, int order)
{
  double alpha = 7*lam_max;

  int  Nx = get_dimension(0);
  int  Ny = get_dimension(1);
  int  Nz = get_dimension(2);
  long Ntot = long(Nx)*long(Ny)*long(Nz);

  DFT_Vec g(Ntot);
  vector<DFT_FFT> v(1); v[0].initialize(Nx,Ny,Nz);
  
  double log_det = Ntot*log(alpha);


  random_device r;
  int seed = r();      
  mt19937 rng(seed);
  uniform_int_distribution<> distrib(0, 1);      
  
  for(int i=1;i<=num_samples;i++)
    {
      // No real difference between these two ...
      //      g.set_random_normal();
      for(long pos = 0; pos < Ntot; pos++) g.set(pos,(distrib(rng) > 0 ? 1 : -1));      
      v[0].Real().set(g); 
      
      vector<DFT_Vec> result(1,Ntot);

      for(int k=1;k<=order;k++)
	{
	  // Compute v-Mv/alpha
	  v[0].do_real_2_fourier();      
	  matrix_dot_v(v,result,NULL);
	  result[0].MultBy(-1.0/alpha);
	  result[0].IncrementBy(v[0].cReal());

	  v[0].Real().set(result[0]); 
	  
	  log_det -= v[0].cReal().dotWith(g)/(k*num_samples);
	}
    }  
  return log_det;   
}


// This is the method of Han, Maliouto, Avron and Shin from https://arxiv.org/pdf/1606.00942.pdf
// It is similar to the above method, but Chebychev interpolation is used since this is supposed to be more accurate.

double Dynamical_Matrix::log_det_2(double lam_max, double lam_min, int num_samples, int order)
{
  double a = lam_min/(lam_min+lam_max);
  double b = lam_max/(lam_min+lam_max);

  int  Nx = get_dimension(0);
  int  Ny = get_dimension(1);
  int  Nz = get_dimension(2);
  long Ntot = long(Nx)*long(Ny)*long(Nz);

  // first, we need the interpolation coefficients for log(x) 
  vector<double> c(order+1);
  for(int k=0;k<=order;k++)
    {
      double xk = cos(M_PI*(k+0.5)/(order+1));	  
      double Tm2 = 0;
      double Tm1 = 0;
      double T   = 1;
      for(int j=0;j<=order;j++)
	{      
	  if(j == 0) T = 1;
	  else if(j == 1) T = xk;
	  else T = 2*xk*Tm1-Tm2;

	  c[j] += log(0.5*(b-a)*xk+0.5*(b+a))*T;

	  Tm2 = Tm1;
	  Tm1 = T;
	}
    }
  for(int j=0;j<=order;j++) c[j] *= (j == 0 ? 1.0 : 2.0)/(order+1);
  /*
  for(double x = a; x <= b; x += (b-a)/10)
    {
      double f = log(x);
      double y = (2*x-b-a)/(b-a);
      double g = 0;
      double Tm2 = 0;
      double Tm1 = 0;
      for(int j=0;j<=order;j++)
	{
	  double T   = 1;      

	  if(j == 1) T = y;
	  if(j > 1)  T = 2*y*Tm1-Tm2;

	  g += c[j]*T;

	  Tm2 = Tm1;
	  Tm1 = T;	  

	}
      cout << x << " " << f << " " << g << endl;
    }
  */
  
  DFT_Vec w0(Ntot);
  DFT_Vec w1(Ntot);
  DFT_Vec w2(Ntot);
  DFT_Vec  v(Ntot);
  DFT_Vec  u(Ntot);
  
  vector<DFT_FFT> work(1); work[0].initialize(Nx,Ny,Nz);
  
  double log_det = Ntot*log(lam_min+lam_max);

  vector<DFT_Vec> result(1,Ntot);

  random_device r;
  int seed = r();      
  mt19937 rng(seed);
  uniform_int_distribution<> distrib(0, 1);      
  
  for(int i=1;i<=num_samples;i++)
    {
      for(long pos = 0; pos < Ntot; pos++) v.set(pos,(distrib(rng) > 0 ? 1 : -1));
      
      work[0].Real().set(v);
      work[0].do_real_2_fourier();
      matrix_dot_v(work,result,NULL);
      result[0].MultBy(1.0/(lam_min+lam_max)); // This is because we use A/( lam_min+ lam_max)	        
      result[0].MultBy(2.0/(b-a));
      result[0].IncrementBy_Scaled_Vector(v, -(b+a)/(b-a));

      w0.set(v);
      w1.set(result[0]); 
      
      u.set(w0);
      u.MultBy(c[0]);
      u.IncrementBy_Scaled_Vector(w1,c[1]);
      
      for(int k=2;k<=order;k++)
	{
	  work[0].Real().set(w1);
	  work[0].do_real_2_fourier();
	  matrix_dot_v(work,result,NULL);
	  result[0].MultBy(1.0/(lam_min+lam_max)); // This is because we use A/( lam_min+ lam_max)	  
	  result[0].MultBy(4.0/(b-a));
	  result[0].IncrementBy_Scaled_Vector(w1, -2*(b+a)/(b-a));
	  result[0].DecrementBy(w0);
	  w2.set(result[0]); 

	  u.IncrementBy_Scaled_Vector(w2,c[k]);

	  w0.set(w1);
	  w1.set(w2);
	}
      log_det += v.dotWith(u)/num_samples;
    }
  return log_det;   
}
