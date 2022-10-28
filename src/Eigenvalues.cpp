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
  eig.set_eigen_vec(xx);
  eig.save_snapshot();
  
  DFT_Vec df(xx.size());
  double f = eig.calculate_gradients(df);

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


void Eigenvalues::save_snapshot()
{
  if(eigen_vec_.size()>1)
  {
    ofstream ofile("snapshot_eigenvector.dat");
    ofile << eigen_vec_;
    ofile.close();
  }
}


double Eigenvalues::calculate_gradients(DFT_Vec& df)
{
  matrix_dot_v(eigen_vec_,df,NULL);

  if(vshift_.size() == df.size())
    df.IncrementBy_Scaled_Vector(vshift_,0.5*vshift_.dotWith(eigen_vec_));
  
  df.MultBy((change_sign_ ? -1 : 1)*scale_);
    
  double x2 = eigen_vec_.dotWith(eigen_vec_);
  double f  = eigen_vec_.dotWith(df)/x2;
  
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
  
  cout << endl; cout << setprecision(12);
  cout << "\tIn Eigenvalues::matrix_dot_v:" << endl;
  cout << "\t  L2-Norm of v:   " << v.euclidean_norm() << endl;
  cout << "\t  L2-Norm of A*v: " << result.euclidean_norm() << endl;
  cout << "\t  L2-Norm ratio |A*v|/|v|:     " << result.euclidean_norm() / v.euclidean_norm() << endl;
  cout << "\t  Dot product vT*A*v/|A*v||v|: " << v.dotWith(result) / result.euclidean_norm() / v.euclidean_norm() << endl;
  
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

  vector<double> x(eigen_vec_.size());
  for(unsigned i=0; i<x.size(); i++) x[i] = eigen_vec_.get(i);
  
  
  ////////////////////////////////////////
  // Minimisation of xAx/xx with FIRE2
  
  // Initialize
  vector<double> v (x.size(), 0.0);
  vector<double> df(x.size(), 0.0);
  vector<double> x_old(x.size());
  for(unsigned i=0; i<x.size(); i++) x_old[i] = x[i];
  
  double f = eigenvalues_objective_func(x, df, this);
  double f_old = f;
  double t = 0.0;
  
  // Parameters
  double alpha_start = 1.0;
  double alpha_min = 0.0;
  double alpha = alpha_start;
  double dt = 0.01;
  double dt_max = 1;
  double dt_min = 0.0;
  double finc = 1.1;
  double falf = 0.9;
  double fdec = 0.1;
  int Npos = 0;
  int Nneg = 0;
  int Nmax = 1e6;
  int Ndelay = 5;
  int Nneg_max = 20;
  bool initialdelay = true;
  
  // Iterate
  for (int i=0; i<Nmax; i++)
  {
    double P = 0.0;
    #pragma omp parallel for reduction(+:P)
    for (long j=0; j<v.size(); j++) P -= df[j]*v[j];
    
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
    
    double P_normalized = P/vnorm/fnorm;
    
    if (P>=0) // && f<=f_old sometimes causes problems sending dt->0
    {
      Npos ++; Nneg = 0;
      
      if (Npos>Ndelay && P_normalized>0.999)
      {
        dt = (dt*finc<dt_max)?dt*finc:dt_max;
        alpha = (alpha*falf>alpha_min)?alpha*falf:alpha_min;
      }
      else if (Npos>Ndelay && P_normalized<=0.5)
      {
        dt = (dt/finc>dt_min)?dt/finc:dt_min;
        alpha = (alpha/falf<alpha_start)?alpha/falf:alpha_start;
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
      
      cout << "\tBacktracking..." << endl;
      #pragma omp parallel for
      for (long j=0; j<v.size(); j++)
      {
        //x[j] -= 0.5*dt*v[j];
        x[j] = x_old[j];
        v[j] = 0.0;
      }
    }
    
    #pragma omp parallel for
    for (long j=0; j<x.size(); j++) x_old[j] = x[j];
    
    // Semi-implicit Euler: Inertial velocity step
    #pragma omp parallel for
    for (long j=0; j<v.size(); j++)
    {
      v[j] -= dt*df[j];
    }
    
    vnorm = 0.0;
    fnorm = 0.0;
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
    
    f_old = f;
    f = eigenvalues_objective_func(x, df, this);
    t += dt;
    
    double xerr  = 0.0;
    #pragma omp parallel for reduction(+:xerr)
    for (long j=0; j<v.size(); j++)
    {
      xerr  += (x_old[j]-x[j])*(x_old[j]-x[j]);
    }
    xerr  = sqrt(xerr)/x.size();
    
    //TODO: compute this and check convergence only 
    //      once in a while to save a few calls??
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
    
    cout << endl; cout << setprecision(12);
    cout << "\tFire2 in Eigenvalues::calculate_eigenvector:" << endl;
    cout << "\t  P/|v||df| = " << P_normalized << endl;
    cout << "\t  dt = " << dt << endl;
    cout << "\t  alpha  = " << alpha << endl;
    cout << "\t  vnorm  = " << vnorm << endl;
    cout << "\t  fnorm  = " << fnorm << endl;
    cout << "\t  f-fold = " << f-f_old << endl;
    cout << "\t  xerr   = " << xerr  << endl;
    cout << "\t  Axerr  = " << Axerr << endl;
    
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
  
  
  ////////////////////////////////////////////////////////////////////
  // Use normalised eigenvector to compute the true eigenvalue
  // This way we cancel the (v^2-1)^2 term in the objective function 
  
  set_eigen_vec(x);
  eigen_vec_.normalise();
  
  DFT_Vec dummy; dummy.zeros(x.size());
  eigen_val_ = calculate_gradients(dummy);
  
  eigen_val_ /= scale_;  
  if(change_sign_) eigen_val_ *= -1;
}



