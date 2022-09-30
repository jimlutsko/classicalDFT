#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <time.h>
#include <random>

using namespace std;

#include <complex>
#include <nlopt.hpp>

#include "myColor.h"
#include "Eigenvalues.h"


static double eigenvalues_objective_func(const std::vector<double> &x, std::vector<double> &grad, void *data)
{
  Eigenvalues& eig = *((Eigenvalues*) data);
  
  eig.set_eigen_vec_from_alias(x);
  
  DFT_Vec df(x.size());

  double f = eig.calculate_gradients(df, x);
  
  df = eig.alias_gradient_from_eigen_vec_gradient(df, x);
  
  if (!grad.empty()) {
    for(long i=0;i<grad.size();i++)
      grad[i] = df.get(i);
  }
  
  return f;
}


// Note: The alias for the density maps -inf to zero density
//       and +inf to the largest meaningful density of 1/dV,
//
//   rho(x) = 0.5*rho_max*(1+tanh(2*x/rho_max))
//   drho(x)/dx = 1/cosh^2( 2*x/rho_max ) = 1-( 2*rho/rho_max -1 )^2


vector<double> Eigenvalues::get_alias_from_eigen_vec()
{
  vector<double> x(eigen_vec_.size());
  
  #ifdef USE_ALIAS_EIGEN
    double rho_max = 1/theDensity_.dV();
    for(long i=0; i<x.size(); i++) x[i] = eigen_vec_.get(i) / (1-pow( 2*theDensity_.get(i)/rho_max -1 ,2));
  #else
    for(unsigned i=0; i<x.size(); i++) x[i] = eigen_vec_.get(i);
  #endif
  
  return x;
}


void Eigenvalues::set_eigen_vec_from_alias(const vector<double> &x)
{
  if(eigen_vec_.size() != x.size()) eigen_vec_.zeros(x.size());
  
  #ifdef USE_ALIAS_EIGEN
    double rho_max = 1/theDensity_.dV();
    for(long i=0; i<x.size(); i++) eigen_vec_.set(i, x[i] * (1-pow( 2*theDensity_.get(i)/rho_max -1 ,2)) );
  #else
    for(long i=0; i<x.size(); i++) eigen_vec_.set(i,x[i]);
  #endif
}


DFT_Vec Eigenvalues::alias_gradient_from_eigen_vec_gradient(const DFT_Vec& df, const vector<double> &x)
{
  DFT_Vec df_x(x.size());
  if (x.size()==0) return df_x;
  
  #ifdef USE_ALIAS_EIGEN
    double rho_max = 1/theDensity_.dV();
    for(long i=0; i<x.size(); i++) df_x.set(i, df.get(i) * (1-pow( 2*theDensity_.get(i)/rho_max -1 ,2)) );
  #else
    df_x = df;
  #endif
  
  return df_x;
}


// Get gradient in density space
double Eigenvalues::calculate_gradients(DFT_Vec& df, const std::vector<double> &xx)
{
  matrix_.matrix_dot_v1(eigen_vec_,df,NULL);

  if(vshift_.size() == df.size())
    df.IncrementBy_Scaled_Vector(vshift_,0.5*vshift_.dotWith(eigen_vec_));
  
  df.MultBy((change_sign_ ? -1 : 1)*scale_);
    
  double v2 = eigen_vec_.dotWith(eigen_vec_);
  double f  = eigen_vec_.dotWith(df)/v2;

  df.IncrementBy_Scaled_Vector(eigen_vec_,-f);
  df.MultBy(2/v2);

  num_eval_++;

  cout << myColor::YELLOW;
  cout << '\r'; cout  <<  "\tEvaluation " << num_eval_ << " gives estimated eigenvalue " << (change_sign_ ? -1 : 1)*f/scale_ << "                 ";
  cout << myColor::RESET;
  
  if (max_num_eval_>0 && num_eval_>=max_num_eval_)
    throw runtime_error("Eigenvalues: Exceeded max number of iterations");

  return f;
}


void Eigenvalues::calculate_eigenvector(Log& theLog)
{
  //TODO: temp
  #ifdef USE_ALIAS_EIGEN
    cout << "\tUsing alias in eigenvector calculation" << endl; cout.flush();
    cout << "\tDensity at pos (1,1,1) is " << theDensity_.get(1,1,1) << endl; cout.flush();
  #else
    cout << "\tNot using alias in eigenvector calculation" << endl; cout.flush();
  #endif
  
  nlopt::opt opt("LD_LBFGS", matrix_.get_Ntot());
  opt.set_min_objective(eigenvalues_objective_func, (void*) this);  
  
  int M_old = opt.get_vector_storage();
  int M_new = 16; opt.set_vector_storage(M_new);
  
  theLog << endl;
  theLog << "\tVector storage resized from " << M_old << " to " << M_new << endl;
  
  if(eigen_vec_.size() != matrix_.get_Ntot())
  {
    eigen_vec_.zeros(matrix_.get_Ntot());
    eigen_vec_.set_random_normal();
  
    if(matrix_.is_fixed_boundary())
    {
      long pos = 0;
      do{eigen_vec_.set(pos,0.0);} while(matrix_.get_next_boundary_point(pos));
    }
  }

  vector<double> x = get_alias_from_eigen_vec();

  //opt.set_ftol_rel(tol_);
  opt.set_xtol_rel(tol_);
  return_code_ = opt.optimize(x, eigen_val_);
    
  eigen_val_ /= scale_;  
  if(change_sign_) eigen_val_ *= -1;

  set_eigen_vec_from_alias(x);
  eigen_vec_.normalise();
}
