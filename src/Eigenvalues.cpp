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

/*
#define OPTIM_USE_OPENMP
#define OPTIM_ENABLE_ARMA_WRAPPERS
#include "optim.hpp"
*/

#ifdef USE_MPI
#include <mpi.h>  
#endif

#include "myColor.h"
#include "Eigenvalues.h"


static double eigenvalues_objective_func(const std::vector<double> &xx, std::vector<double> &grad, void *data)
{
  Eigenvalues& eig = *((Eigenvalues*) data);

  if (eig.already_using_density_alias()) eig.set_eigen_vec(xx);
  else eig.set_eigen_vec_from_alias_vector(xx);
  
  DFT_Vec df(xx.size());

  double f = eig.calculate_gradients(df);
  
  if (!eig.already_using_density_alias()) eig.get_species()->convert_to_alias_deriv(df);

  if (!grad.empty()) {
    #pragma omp parallel for
    for(long i=0;i<grad.size();i++)
      grad[i] = df.get(i);
  }
  return f;
}


void Eigenvalues::set_eigen_vec(const vector<double> &v_in)
{
  if(eigen_vec_.size() != v_in.size()) eigen_vec_.zeros(v_in.size());
  for(long i=0;i<v_in.size();i++) eigen_vec_.set(i,v_in[i]);
}


void Eigenvalues::set_eigen_vec_from_alias_vector(const vector<double> &v_in)
{
  if(eigen_vec_.size() != v_in.size()) eigen_vec_.zeros(v_in.size());
  for(long i=0;i<v_in.size();i++) eigen_vec_.set(i,v_in[i]);
  species_->convert_to_density_increment(eigen_vec_);
}


double Eigenvalues::calculate_gradients(DFT_Vec& df)
{  
  matrix_dot_v(eigen_vec_,df,NULL);

  if(vshift_.size() == df.size())
    df.IncrementBy_Scaled_Vector(vshift_,0.5*vshift_.dotWith(eigen_vec_));
  
  df.MultBy((change_sign_ ? -1 : 1)*scale_);
    
  double x2 = eigen_vec_.dotWith(eigen_vec_); //TODO: include metric??
  double f  = eigen_vec_.dotWith(df)/x2;
  /*
  double x2, f;
  if (already_using_density_alias())
  {
    DFT_Vec x  = eigen_vec_; species_->convert_to_density_increment(x);  //species_->convert_to_density_increment(x);
    DFT_Vec Ax = df;         species_->convert_to_density_increment(Ax); //species_->convert_to_density_increment(Ax);
    
    x2 = x.dotWith(x);
    f  = x.dotWith(Ax)/x2;
  }
  else
  {
    x2 = eigen_vec_.dotWith(eigen_vec_);
    f  = eigen_vec_.dotWith(df)/x2;
  }
  */
  df.IncrementBy_Scaled_Vector(eigen_vec_,-f);
  df.MultBy(2/x2);
  
  // Artificial cost to keep the vectors near the unit sphere (does not change the eigenproblem)
  double vnorm2 = eigen_vec_.dotWith(eigen_vec_);
  f += (vnorm2-1)*(vnorm2-1);
  df.IncrementBy_Scaled_Vector(eigen_vec_,4*(vnorm2-1));

  num_eval_++;

  cout << myColor::YELLOW;
  cout << '\r'; cout  <<  "\tEvaluation " << num_eval_ << " gives estimated eigenvalue " << (change_sign_ ? -1 : 1)*f/scale_ << "                 ";
  cout << myColor::RESET;
  
  if (max_num_eval_>0 && num_eval_>=max_num_eval_)
    throw runtime_error("Eigenvalues: Exceeded max number of iterations");
  
  return f;
}


void Eigenvalues::matrix_dot_v(const DFT_Vec &v, DFT_Vec &result, void *param) const
{
  if (is_using_density_alias())
  {
    DFT_Vec vv; vv.set(v);
    
    species_->convert_to_density_increment(vv);
    matrix_.matrix_dot_v1(vv,result,param);
    species_->convert_to_alias_deriv(result);
    
    DFT_Vec df; df.set(species_->getDF());
    DFT_Vec d2Rhodx2; species_->get_second_derivatives_of_density_wrt_alias(d2Rhodx2);
    df.Schur(df, d2Rhodx2);
    df.Schur(df, v);
    
    result.IncrementBy(df);
  }
  else
  {
    matrix_.matrix_dot_v1(v,result,param);
  }
}


void Eigenvalues::matrix_dot_v(const vector<double> &vv, vector<double> &result, void *param)
{
  DFT_Vec v; v.zeros(vv.size());
  DFT_Vec r; r.zeros(vv.size());
  result.resize(vv.size());
  
  #pragma omp parallel for
  for (long i=0; i<vv.size(); i++) v.set(i, vv[i]);
  
  matrix_dot_v(v, r, param);
  
  #pragma omp parallel for
  for (long i=0; i<vv.size(); i++) result[i] = r.get(i);
}


/*
void Eigenvalues::calculate_eigenvector(Log& theLog)
{
  //nlopt::opt opt("LD_LBFGS", matrix_.get_Ntot());
  nlopt::opt opt("LD_MMA", matrix_.get_Ntot());
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
      
      eigen_vec_.normalise();
  }
  /*
  else // perturb --> TODO: why does that help??
  {
      DFT_Vec rnd_vec;
      rnd_vec.zeros(matrix_.get_Ntot());
      rnd_vec.set_random_normal();
      
      eigen_vec_.IncrementBy_Scaled_Vector(rnd_vec, 0.1);
      
      if(matrix_.is_fixed_boundary())
      {
        long pos = 0;
        do{eigen_vec_.set(pos,0.0);} while(matrix_.get_next_boundary_point(pos));
      }
      
      eigen_vec_.normalise();
  }
  *//*

  vector<double> x(eigen_vec_.size());
  if (already_using_density_alias()) 
  {
    for(unsigned i=0; i<x.size();i++) x[i] = eigen_vec_.get(i);
  }
  else 
  {
    DFT_Vec eigen_vec_alias; eigen_vec_alias.set(eigen_vec_);
    species_->convert_to_alias_increment(eigen_vec_alias);
    
    for(unsigned i=0; i<x.size();i++) x[i] = eigen_vec_alias.get(i);
  }

  //opt.set_ftol_rel(tol_);
  opt.set_xtol_rel(tol_);
  return_code_ = opt.optimize(x, eigen_val_);
    
  eigen_val_ /= scale_;  
  if(change_sign_) eigen_val_ *= -1;

  if (already_using_density_alias()) set_eigen_vec(x);
  else set_eigen_vec_from_alias_vector(x);

  eigen_vec_.normalise();
}*/


void Eigenvalues::calculate_eigenvector(Log& theLog)
{
  if(eigen_vec_.size() != matrix_.get_Ntot())
  {
      eigen_vec_.zeros(matrix_.get_Ntot());
      eigen_vec_.set_random_normal();
  
      if(matrix_.is_fixed_boundary())
      {
        long pos = 0;
        do{eigen_vec_.set(pos,0.0);} while(matrix_.get_next_boundary_point(pos));
      }
      
      eigen_vec_.normalise();
  }
  /*
  else // perturb --> TODO: why does that help??
  {
      DFT_Vec rnd_vec;
      rnd_vec.zeros(matrix_.get_Ntot());
      rnd_vec.set_random_normal();
      
      eigen_vec_.IncrementBy_Scaled_Vector(rnd_vec, 0.1);
      
      if(matrix_.is_fixed_boundary())
      {
        long pos = 0;
        do{eigen_vec_.set(pos,0.0);} while(matrix_.get_next_boundary_point(pos));
      }
      
      eigen_vec_.normalise();
  }
  */

  vector<double> x(eigen_vec_.size());
  if (already_using_density_alias()) 
  {
    for(unsigned i=0; i<x.size();i++) x[i] = eigen_vec_.get(i);
  }
  else 
  {
    DFT_Vec eigen_vec_alias; eigen_vec_alias.set(eigen_vec_);
    species_->convert_to_alias_increment(eigen_vec_alias);
    
    for(unsigned i=0; i<x.size();i++) x[i] = eigen_vec_alias.get(i);
  }
  
  ////////////////////////////////////
  
  // Initialize
  vector<double> v (x.size(), 0.0);
  vector<double> df(x.size(), 0.0);
  vector<double> x_old(x.size(), 0.0);
  double f = eigenvalues_objective_func(x, df, this);
  double t = 0.0;
  
  // Parameters
  double alpha_start = 0.5;
  double alpha = alpha_start;
  double dt = 0.01;
  double dt_max = 0.1;
  double dt_min = 0.0;
  double finc = 1.1;
  double falf = 1.0;
  double fdec = 0.5;
  int Npos = 0;
  int Nneg = 0;
  int Nmax = 1e6;
  int Ndelay = 10;
  int Nneg_max = 10;
  bool initialdelay = true;
  
  // Iterate
  for (int i=0; i<Nmax; i++)
  {
    double P = 0.0; 
    #pragma omp parallel for reduction(+:P)
    for (long j=0; j<v.size(); j++) P -= df[j]*v[j];
    
    #pragma omp parallel for
    for (long j=0; j<x.size(); j++) x_old[j] = x[j];
    
    if (P>0)
    {
      Npos ++; Nneg = 0;
      
      if (Npos>Ndelay)
      {
        dt = (dt*finc<dt_max)?dt*finc:dt_max;
        alpha *= falf;
      }
    }
    else // if P<=0
    {
      Nneg ++; Npos = 0;
      
      if (Nneg>Nneg_max) throw runtime_error("Eigenvalues.cpp: Cannot stop going uphill in FIRE2");
      
      if (!initialdelay || i>=Ndelay) 
      {
        if (dt*fdec>dt_min) dt *= fdec;
        alpha = alpha_start;
      }
      
      #pragma omp parallel for
      for (long j=0; j<v.size(); j++)
      {
        x[j] -= 0.5*dt*v[j];
        v[j] = 0.0;
      }
    }
    
    // Semi-implicit Euler: Inertial velocity step
    #pragma omp parallel for
    for (long j=0; j<v.size(); j++)
    {
      v[j] -= dt*df[j];
    }
    
    double vnorm = 0.0;
    double fnorm = 0.0;
    #pragma omp parallel for reduction(+:vnorm) reduction(+:fnorm)
    for (long j=0; j<v.size(); j++)
    {
      vnorm += v[j]*v[j];
      fnorm += df[j]*df[j];
    }
    vnorm = sqrt(vnorm);
    fnorm = sqrt(fnorm);
    
    // Semi-implicit Euler: velocity mixing and x-step
    #pragma omp parallel for
    for (long j=0; j<v.size(); j++)
    {
      v[j] = (1-alpha)*v[j] - alpha*df[j]*vnorm/fnorm;
      x[j] += dt*v[j];
    }
    
    double f_old = f;
    f = eigenvalues_objective_func(x, df, this);
    t += dt;
    
    double xerr  = 0.0;
    #pragma omp parallel for reduction(+:xerr)
    for (long j=0; j<v.size(); j++)
    {
      xerr  += (x_old[j]-x[j])*(x_old[j]-x[j]);
    }
    xerr  = sqrt(xerr)/x.size();
    
    vector<double> Ax; matrix_dot_v(x, Ax, NULL);
    double xdotx, xdotAx, AxdotAx; xdotx = xdotAx = AxdotAx = 0.0;
    #pragma omp parallel for reduction(+:xdotx) reduction(+:xdotAx) reduction(+:AxdotAx)
    for (long j=0; j<x.size(); j++)
    {
      xdotx   +=  x[j]* x[j];
      xdotAx  +=  x[j]*Ax[j];
      AxdotAx += Ax[j]*Ax[j];
    }
    double Axerr = 1-fabs(xdotAx)/sqrt(xdotx)/sqrt(AxdotAx);
    
    cout << endl;
    cout << "\tFire2 P  = " << P  << endl;
    cout << "\tFire2 dt = " << dt << endl;
    cout << "\tFire2 vnorm = " << vnorm << endl;
    cout << "\tFire2 fnorm = " << fnorm << endl;
    cout << "\tFire2 Axerr = " << Axerr << endl;
    cout << "\tFire2 xerr  = " << xerr  << endl;
    cout << "\tFire2 f-fold = " << f-f_old << endl;
    
    // Check if converged
    if (!initialdelay || i>=Ndelay)
    {
      //cout << endl << "convergence monitor = " << fnorm/df.size()*dt*dt << endl;
      //cout << "\rerr = " << fnorm/df.size()*dt*dt << endl;
      //if (P>0 && fnorm/df.size()*dt*dt<tol_) break;
      //if (P>0 && fabs(f-f_old)/fabs(f+f_old)<tol_) break;
      //if (P>0 && xerr<tol_) break;
      if (Axerr < tol_) break;
    }
  }
  
  // Finalize
  eigen_val_ = f;
  
  ////////////////////////////////////
    
  eigen_val_ /= scale_;  
  if(change_sign_) eigen_val_ *= -1;

  if (already_using_density_alias()) set_eigen_vec(x);
  else set_eigen_vec_from_alias_vector(x);

  eigen_vec_.normalise();
}
