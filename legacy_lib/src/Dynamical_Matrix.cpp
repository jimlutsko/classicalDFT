


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


void Dynamical_Matrix::set_boundary_points_to_zero(DFT_Vec &v) const
{
  long pos = 0;
  
  do{
    v.set(pos,0.0);
  } while(get_next_boundary_point(pos));
  
}


void Dynamical_Matrix::matrix_dot_v1(const DFT_Vec &v, DFT_Vec &result, void *param) const
{
  vector<DFT_FFT> vwork(1);
  vwork[0].initialize(get_dimension(0), get_dimension(1), get_dimension(2));
  vwork[0].Real().set(v);

  if(is_fixed_boundary()) set_boundary_points_to_zero(vwork[0].Real());
  
  vwork[0].do_real_2_fourier();

  vector<DFT_Vec> rwork(1, result.size());

  matrix_dot_v(vwork,rwork,param);
  result.set(rwork[0]);

  if(is_fixed_boundary()) set_boundary_points_to_zero(result);
}

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
  DFT_Vec v(Ntot);
  DFT_Vec result(Ntot);
  
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
      if(is_fixed_boundary()) set_boundary_points_to_zero(g);

      v.set(g);       

      for(int k=1;k<=order;k++)
	{
	  // Compute v-Mv/alpha
	  matrix_dot_v1(v,result,NULL);
	  result.MultBy(-1.0/alpha);
	  result.IncrementBy(v);

	  v.set(result); 
	  
	  log_det -= v.dotWith(g)/(k*num_samples);
	}
    }  
  return log_det;   
}


// This is the method of Han, Maliouto, Avron and Shin from https://arxiv.org/pdf/1606.00942.pdf
// It is similar to the above method, but Chebychev interpolation is used since this is supposed to be more accurate.
// If do_condition = true, there is a zero eigenvalue associated with vector v0 = (1/sqrt(Ntot),...,1/sqrt(Ntot)) and so we use A->A+lam_max v0v0
double Dynamical_Matrix::log_det_2(double lam_max, double lam_min, int num_samples, int order, bool has_zero_eigenvalue)
{
  double scale = 1.0/(lam_min+lam_max);
  double a = lam_min*scale;
  double b = lam_max*scale;

  double lam_mid = (lam_max+lam_min)/2;

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
  
  long Ntot    = get_Ntot();
  long Nactive = Ntot - (is_fixed_boundary() ? get_Nboundary() : 0);

  DFT_Vec w0(Ntot); w0.zeros();
  DFT_Vec w1(Ntot); w1.zeros();
  DFT_Vec w2(Ntot); w2.zeros();
  DFT_Vec  v(Ntot);  v.zeros();
  DFT_Vec  u(Ntot);  u.zeros();
  DFT_Vec result(Ntot); result.zeros();
  
  double log_det = -Nactive*log(scale);

  random_device r;
  int seed = r();      
  mt19937 rng(seed);
  uniform_int_distribution<> distrib(0, 1);      
  
  for(int i=1;i<=num_samples;i++)
    {
      for(long pos = 0; pos < Ntot; pos++) v.set(pos,(distrib(rng) > 0 ? 1 : -1));
      if(is_fixed_boundary()) set_boundary_points_to_zero(v);
      
      matrix_dot_v1(v,result,NULL);      
      if(is_dynamic())        result.MultBy(-1);
      if(has_zero_eigenvalue) result.add(lam_mid*v.accu()/Ntot);
	
      result.MultBy(scale); // This is because we use A/( lam_min+ lam_max)	        
      result.MultBy(2.0/(b-a));
      result.IncrementBy_Scaled_Vector(v, -(b+a)/(b-a));

      w0.set(v);
      w1.set(result); 
      
      u.set(w0);
      u.MultBy(c[0]);
      u.IncrementBy_Scaled_Vector(w1,c[1]);
      
      for(int k=2;k<=order;k++)
	{
	  matrix_dot_v1(w1,result,NULL);
	  if(is_dynamic())        result.MultBy(-1); 
	  if(has_zero_eigenvalue) result.add(lam_mid*w1.accu()/Ntot);
	  
	  result.MultBy(scale); // This is because we use A/( lam_min+ lam_max)	  
	  result.MultBy(4.0/(b-a));
	  result.IncrementBy_Scaled_Vector(w1, -2*(b+a)/(b-a));
	  result.DecrementBy(w0);
	  w2.set(result);

	  u.IncrementBy_Scaled_Vector(w2,c[k]);

	  w0.set(w1);
	  w1.set(w2);
	}
      log_det += v.dotWith(u)/num_samples;
      if(dm_verbose_)
	cout << "\t log_det = " << (has_zero_eigenvalue ? -log(lam_mid) : 0)  + num_samples*log_det/i << endl;
    }

  if(has_zero_eigenvalue) log_det -= log(lam_mid);

  return log_det;   
}

// This is the method of Han, Maliouto, Avron and Shin from https://arxiv.org/pdf/1606.00942.pdf
// It is similar to the above method, but Chebychev interpolation is used since this is supposed to be more accurate.

double Dynamical_Matrix::log_fabs_det_2(double lam_max, double lam_min, int num_samples, int order, bool has_zero_eigenvalue)
{
  double scale = 1.0/(lam_min*lam_min+lam_max*lam_max);
  double a = lam_min*lam_min*scale;
  double b = lam_max*lam_max*scale;
  
  double lam_mid = lam_max*lam_max*scale;  

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

  long Ntot = get_Ntot();
  long Nactive = Ntot - (is_fixed_boundary() ? get_Nboundary() : 0);
  
  DFT_Vec w0(Ntot);
  DFT_Vec w1(Ntot);
  DFT_Vec w2(Ntot);
  DFT_Vec  v(Ntot);
  DFT_Vec  u(Ntot);
  DFT_Vec result(Ntot);

  //  double log_det = -Nactive*log(scale);
  double log_det = 0;
  
  random_device r;
  int seed = r();      
  mt19937 rng(seed);
  uniform_int_distribution<> distrib(0, 1);      
  
  for(int i=1;i<=num_samples;i++)
    {
      for(long pos = 0; pos < Ntot; pos++) v.set(pos,(distrib(rng) > 0 ? 1 : -1));
      if(is_fixed_boundary()) set_boundary_points_to_zero(v);
	
      matrix_dot_v1(v,result,NULL);
      matrix_dot_v1(result,result,NULL);

      if(has_zero_eigenvalue) result.add(lam_mid*v.accu()/Ntot);      

      result.MultBy(scale); // This is because we use A/( lam_min+ lam_max)	        
      result.MultBy(2.0/(b-a));
      result.IncrementBy_Scaled_Vector(v, -(b+a)/(b-a));

      w0.set(v);
      w1.set(result); 
      
      u.set(w0);
      u.MultBy(c[0]);
      u.IncrementBy_Scaled_Vector(w1,c[1]);
      
      for(int k=2;k<=order;k++)
	{
	  matrix_dot_v1(w1,result,NULL);
	  matrix_dot_v1(result,result,NULL);

	  if(has_zero_eigenvalue) result.add(lam_mid*w1.accu()/Ntot);
	  
	  result.MultBy(scale); // This is because we use A/( lam_min+ lam_max)	  
	  result.MultBy(4.0/(b-a));
	  result.IncrementBy_Scaled_Vector(w1, -2*(b+a)/(b-a));
	  result.DecrementBy(w0);
	  w2.set(result);

	  u.IncrementBy_Scaled_Vector(w2,c[k]);

	  w0.set(w1);
	  w1.set(w2);
	}
      log_det += v.dotWith(u)/num_samples;
      if(dm_verbose_)
	cout << "\t log_det = " << -0.5*Nactive*log(scale) + (has_zero_eigenvalue ? -0.5*log(lam_mid) : 0) + 0.5*num_samples*log_det/i << endl;      
    }

  if(has_zero_eigenvalue) log_det -= log(lam_mid);
  log_det -= Nactive*log(scale);
  
  return log_det/2;   
}
