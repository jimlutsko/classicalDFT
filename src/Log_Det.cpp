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

#include <complex>

#include "myColor.h"

#include "Log_Det.h"

Log_Det::Log_Det(const Dynamical_Matrix& matrix, int order, double lam_max, double lam_min, bool verbose)
  : matrix_(matrix), verbose_(verbose), lam_max_(lam_max), lam_min_(lam_min)
{
  //  if(matrix_.is_dynamic()) { lam_min_ *= lam_min_; lam_max_ *= lam_max_;}
  //  lam_min_ *= lam_min_; lam_max_ *= lam_max_;
  
  a_     = lam_min_;
  b_     = lam_max_;

  cout << "a = " << a_ << " b_ = " << b_ << endl;
  
  get_coefficients(order);
}

void Log_Det::get_coefficients(int order)
{
  c_.resize(order+1,0.0);
  
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

	  c_[j] += log(0.5*(b_-a_)*xk+0.5*(b_+a_))*T;

	  Tm2 = Tm1;
	  Tm1 = T;
	}
    }
  for(int j=0;j<=order;j++) c_[j] *= (j == 0 ? 1.0 : 2.0)/(order+1);
  /*
  for(int k=0;k<=order;k++)
    {
      cout << "c_[" << k << "] = " << c_[k] << endl;
    }
  cout << endl;
  */
}

void Log_Det::matrix_dot_v1(const DFT_Vec &v, DFT_Vec &result) const
{
  //  DFT_Vec r1(v);
  //matrix_.matrix_dot_v1(v,r1,NULL);  // r1 = (gH)v
  //matrix_.matrix_dot_v1(r1,result,NULL); // result = (gH)r1 = (gH)(gH)v

  matrix_.matrix_dot_v1(v,result,NULL);  // r1 = (gH)v

  int s = 1;
  if(matrix_.is_dynamic()) s = -1;
  result.MultBy(s*2.0/(b_-a_));
  result.IncrementBy_Scaled_Vector(v, -(b_+a_)/(b_-a_));
}


double Log_Det::eval_poly(double x) const
{
  //  double y = (2*x/(lam_max_-lam_min_))-(lam_max_+lam_min_)/(lam_max_-lam_min_);
  double y = (2*x/(b_ - a_))-(b_+a_)/(b_-a_);
  double z = 0;

  double T0 = 1;
  double T1 = y;
  z = c_[0]*T0+c_[1]*T1;
  for(int k=2;k<c_.size();k++)
    {
      double T = 2*y*T1-T0;
      z += c_[k]*T;
      
      T0 = T1;
      T1 = T;
    }
  return z;
}

double Log_Det::calculate_log_det(long seed, int num_samples, bool has_zero_eigenvalue, double &variance, vector<double> eigenvalues)
{
  double lam_mid = (lam_max_+lam_min_)/2;

  // test
  for(double x=log(lam_min_); x<=log(lam_max_);x+= (log(lam_max_)-log(lam_min_))/10)
    {
      double xx = exp(x);
      double z  = eval_poly(xx);
      cout << "\tcheck: x = "<< xx << " z=sum c_k*T_k(y) = " << z << " ln(x) = " << log(xx) << endl;
    }
  cout << endl;

  cout << "lam_min = "<< lam_min_ << " lam_max = " << lam_max_ << endl;
  cout << "Has zero eigenvalue = " << has_zero_eigenvalue << endl << endl;
  cout << "Fixed boundary      = " << matrix_.is_fixed_boundary() << endl;
  cout << "Is Dynamic          = " << matrix_.is_dynamic() << endl;

  long Ntot    = matrix_.get_Ntot();
  long Nactive = Ntot - (matrix_.is_fixed_boundary() ? matrix_.get_Nboundary() : 0);

  DFT_Vec w0(Ntot); w0.zeros();
  DFT_Vec w1(Ntot); w1.zeros();
  DFT_Vec w2(Ntot); w2.zeros();
  DFT_Vec  v(Ntot);  v.zeros();
  DFT_Vec  u(Ntot);  u.zeros();
  DFT_Vec result(Ntot); result.zeros();
  
  double log_det = 0;
  double var_log_det = 0;
 
  if(seed <= 0) { random_device r;  seed = r();}
  mt19937 rng(seed);
  uniform_int_distribution<> distrib(0, 1);      

  for(int i=1;i<=num_samples;i++)
    {
      for(long pos = 0; pos < Ntot; pos++) v.set(pos,(distrib(rng) > 0 ? 1 : -1));
      if(matrix_.is_fixed_boundary()) matrix_.set_boundary_points_to_zero(v);

      matrix_dot_v1(v,result);      
      //      if(has_zero_eigenvalue) result.add(lam_mid*v.accu()/Ntot);

      w0.set(v);
      w1.set(result); 

      u.set(w0);
      u.MultBy(c_[0]);
      u.IncrementBy_Scaled_Vector(w1,c_[1]);
      
      for(int k=2;k<c_.size();k++)
	{
	  if(i == 1)
	    {
	      cout << myColor::YELLOW;
	      cout << setprecision(6);      
	      cout << '\r'; cout << "\t First sample: evaluating order " << k
		   << " of " << c_.size() <<  "                         "; cout.flush();
	      cout << myColor::RESET;
	    }

	  // T_n = 2*A*T_{n-1}-T_{n-2}
	  matrix_dot_v1(w1,result);
	  //	  if(has_zero_eigenvalue) result.add(lam_mid*w1.accu()/Ntot);
	  result.MultBy(2);	  
	  result.DecrementBy(w0);
	  w2.set(result);

	  u.IncrementBy_Scaled_Vector(w2,c_[k]);

	  w0.set(w1);
	  w1.set(w2);
	}
      
      double val = v.dotWith(u);
      double val1 = val;
      for(double& x: eigenvalues)
	{
	  val -= eval_poly(x);
	  //	  if(x > 1e-3) val += log(fabs(x));
	  val += log(fabs(x));
	}
      
      log_det     += val;
      var_log_det += val*val;

      double av          = log_det/i;
      double av2         = var_log_det/i;
      double current_val = av;

      variance = sqrt(fabs(av2-av*av));
      
      cout << myColor::YELLOW;
      cout << setprecision(6);
      if(i == 1) cout << endl;
      cout << "\t samples = " << i << " out of " << num_samples << " log_det = " << current_val 
	   << " stderr = " << sqrt(fabs(av2-av*av))/sqrt(i) << " : " << val1 << " " << val;
      cout << endl; 
      cout << myColor::RESET;	      
    }
  cout << endl;

  log_det /= num_samples;				     
  
  return log_det;
}
