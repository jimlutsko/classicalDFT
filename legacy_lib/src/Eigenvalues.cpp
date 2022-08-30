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


static double eigenvalues_objective_func(const std::vector<double> &xx, std::vector<double> &grad, void *data)
{
  Eigenvalues& eig = *((Eigenvalues*) data);

  eig.set_eigen_vec_(xx);
  
  DFT_Vec df(xx.size());
  
  double f = eig.calculate_gradients(df);

  if (!grad.empty()) {
    for(long i=0;i<grad.size();i++)
      grad[i] = df.get(i);
  }
  return f;
}

void Eigenvalues::set_eigen_vec_(const vector<double> &v_in)
{
  if(eigen_vec_.size() != v_in.size()) eigen_vec_.zeros(v_in.size());  

  for(long i=0;i<v_in.size();i++) eigen_vec_.set(i,v_in[i]);
}


double Eigenvalues::calculate_gradients(DFT_Vec& df)
{
  matrix_.matrix_dot_v1(eigen_vec_,df,NULL);

  df.MultBy((change_sign_ ? -1 : 1)*scale_);
    
  double x2 = eigen_vec_.dotWith(eigen_vec_);      
  double f  = eigen_vec_.dotWith(df)/x2;

  df.IncrementBy_Scaled_Vector(eigen_vec_,-f);
  df.MultBy(2/x2);

  num_eval_++;

  cout << myColor::YELLOW;
  cout << '\r'; cout  <<  "\tEvaluation " << num_eval_ << " gives estimated eigenvalue " << (change_sign_ ? -1 : 1)*f/scale_ << "                 ";
  cout << myColor::RESET;

  return f;
}

void Eigenvalues::calculate_eigenvector(Log& theLog)
{
  nlopt::opt opt("LD_LBFGS", matrix_.get_Ntot());
  opt.set_min_objective(eigenvalues_objective_func, (void*) this);  

  eigen_vec_.zeros(matrix_.get_Ntot());
  eigen_vec_.set_random_normal();

  vector<double> x(eigen_vec_.size());
  for(unsigned i=0; i<x.size();i++) x[i] = eigen_vec_.get(i);

  opt.set_ftol_rel(tol_);
  return_code_ = opt.optimize(x, eigen_val_);
    
  eigen_val_ /= scale_;  
  if(change_sign_) eigen_val_ *= -1;

  eigen_vec_.normalise();
  theLog << endl;
}
