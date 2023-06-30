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

  if (!grad.empty()) 
  {
    #ifdef USE_OMP
    #pragma omp parallel for
    #endif
    for(long i=0;i<grad.size();i++) grad[i] = df.get(i);
  }
  
  return f;
}


static double eigenvalues_objective_func(const DFT_Vec &xx, DFT_Vec &grad, void *data)
{
  Eigenvalues& eig = *((Eigenvalues*) data);
  eig.set_eigen_vec(xx);
  eig.save_snapshot();
  
  grad.zeros(xx.size());
  return eig.calculate_gradients(grad);
}


DFT_Vec Eigenvalues::get_density_space_eigenvector() const
{
  DFT_Vec v(eigen_vec_);
  
  if (matrix_.is_dynamic()) g_dot_v(eigen_vec_, v);
  else if (is_using_density_alias()) species_->convert_to_density_increment(v);
  
  return v;
}


void Eigenvalues::set_from_density_space_eigenvector(const DFT_Vec &v)
{
  if(eigen_vec_.size() != v.size()) eigen_vec_.zeros(v.size());
  
  matrix_.matrix_dot_v1(v, eigen_vec_, NULL, true); // (d2F)v
  if(vnormal_.size() == v.size()) eigen_vec_.IncrementBy_Scaled_Vector(vnormal_, -vnormal_.dotWith(eigen_vec_));
  
  eigen_vec_.MultBy(1/norm(eigen_vec_));
}


void Eigenvalues::set_eigen_vec(const vector<double> &v_in)
{
  if(eigen_vec_.size() != v_in.size()) eigen_vec_.zeros(v_in.size());
  
  #ifdef USE_OMP
  #pragma omp parallel for
  #endif
  for(long i=0;i<v_in.size();i++) eigen_vec_.set(i,v_in[i]);
  
  if(vnormal_.size() == v_in.size()) eigen_vec_.IncrementBy_Scaled_Vector(vnormal_, -vnormal_.dotWith(eigen_vec_));
}


void Eigenvalues::g_dot_v(const DFT_Vec &v, DFT_Vec &gv) const
{
  if (gv.size() != v.size()) gv.zeros(v.size());
  
  matrix_.g_dot_x(v, gv);
  gv.MultBy(-1); //Note: minus sign because original metric g is not positive definite
}

double Eigenvalues::v_dot_w(const DFT_Vec &v, const DFT_Vec &w) const
{
  double vw = 0.0;
  
  // Here I copy twice the code that performs the actual dot calculation
  // The reason is that in case I work in a subspace orthogonal to vnormal_,
  // I have to make a copy of both input vectors, which is a computationally 
  // expensive operation that can be avoided in the more common cases where
  // we don't use this feature.
  
  if(vnormal_.size() == v.size()) 
  {
    DFT_Vec vv = v; vv.IncrementBy_Scaled_Vector(vnormal_, -vnormal_.dotWith(vv));
    DFT_Vec ww = w; ww.IncrementBy_Scaled_Vector(vnormal_, -vnormal_.dotWith(ww));
    
    /// Actual dot calculation
    if (matrix_.is_dynamic())
    {
      DFT_Vec gw(ww); g_dot_v(ww, gw);
      vw = vv.dotWith(gw);
    }
    else
    {
      vw = vv.dotWith(ww);
    }
    ///
  }
  else
  {
    /// Actual dot calculation
    if (matrix_.is_dynamic())
    {
      DFT_Vec gw(w); g_dot_v(w, gw);
      vw = v.dotWith(gw);
    }
    else
    {
      vw = v.dotWith(w);
    }
    ///
  }
  
  return vw;
}


double Eigenvalues::norm(const DFT_Vec &v) const
{
  return sqrt(v_dot_w(v,v));
}


void Eigenvalues::save_snapshot() const
{
  if (eigen_vec_.size()>1)
  {
    ofstream ofile("snapshot_eigenvector.dat");
    ofile << eigen_vec_;
    ofile.close();
    
    ofile.open("snapshot_density_eigenvector.dat");
    ofile << get_density_space_eigenvector();
    ofile.close();
  }
}


// Note: df has a lower index and we keep it that way to make subsequent 
// calculations easier. The true gradient is g^{IJ}(df)_J and is contravariant

double Eigenvalues::calculate_gradients(DFT_Vec& df)
{
  matrix_dot_v(eigen_vec_,df,NULL);
  
  if(vshift_.size() == df.size())
    df.IncrementBy_Scaled_Vector(vshift_,0.5*v_dot_w(vshift_,eigen_vec_));
  
  matrix_dot_eigen_vec_ = df;
  df.MultBy((change_sign_ ? -1 : 1)*scale_);

  double x2 = v_dot_w(eigen_vec_,eigen_vec_);
  double f  = v_dot_w(eigen_vec_,df)/x2;
  
  eigen_val_ = (change_sign_ ? -1 : 1)*f/scale_;
  
  df.IncrementBy_Scaled_Vector(eigen_vec_,-f);
  df.MultBy(2/x2);
  
  // Artificial cost to keep the vectors near the unit sphere (does not change the eigenproblem)
  // Do not scale this, otherwise the vAv/vv term for small eigenvalues is dwarfed next to this one
  f += (x2-1)*(x2-1);
  df.IncrementBy_Scaled_Vector(eigen_vec_,4*(x2-1));
  
  // Raise index of df
  //if (matrix_.is_dynamic()) g_dot_v(eigen_vec_, df);
  
  if (verbose_) cout  <<  "\tObjective function is now " << f << "                 " << endl;
  
  num_eval_++;
  
  cout << myColor::YELLOW;
  cout << '\r'; cout  <<  "\tEvaluation " << num_eval_ << " gives estimated eigenvalue " << eigen_val_ << " (res = " << calculate_residual_error(false) << ")                 ";
  cout << myColor::RESET;
  
  if (max_num_eval_>0 && num_eval_>=max_num_eval_)
    throw runtime_error("Eigenvalues: Exceeded max number of iterations");
  
  return f;
}


double Eigenvalues::calculate_residual_error(bool recompute_matrix_dot_v, bool first_order_correction_to_eigenvalue) const
{
  if (eigen_vec_.size() != matrix_.get_Ntot())
    throw runtime_error("Eigenvalues: Eigenvector not initialized yet!");
  
  DFT_Vec residual(matrix_dot_eigen_vec_);
  
  if (recompute_matrix_dot_v ||
      matrix_dot_eigen_vec_.size() != matrix_.get_Ntot())
  {
    residual.zeros(matrix_.get_Ntot());
    matrix_dot_v(eigen_vec_, residual, NULL);
  }
  
  double eigenvalue = eigen_val_;
  if (first_order_correction_to_eigenvalue)
  {
    eigenvalue = v_dot_w(eigen_vec_,residual);
    eigenvalue /= norm(eigen_vec_)*norm(eigen_vec_);
  }
  
  // Note that at this stage "residual" is still just A*v
  double alpha = v_dot_w(eigen_vec_,residual)
                 /norm(eigen_vec_)/norm(residual);
  
  residual.IncrementBy_Scaled_Vector(eigen_vec_, -eigenvalue); // now "residual" is A*v-x*v
  double r = norm(residual)/norm(eigen_vec_)/fabs(eigenvalue);
  
  double err = (r*r+1) - pow(alpha,-2);
  //if (fabs(err)>1e-8) throw runtime_error("Eigenvalues: Inconsistent residual calculation (r^2+1="+to_string(r*r+1)+" while 1/alpha^2="+to_string(pow(alpha,-2)));
  if (fabs(err)>1e-8) cerr << "Eigenvalues: Inconsistent residual calculation (r^2+1=" << r*r+1 << " while 1/alpha^2=" << pow(alpha,-2) << ")" << endl;;
  
  return r;
}


void Eigenvalues::matrix_dot_v(const DFT_Vec &v_in, DFT_Vec &result, void *param) const
{
  if (is_using_density_alias() && matrix_.is_dynamic())
    throw runtime_error("(In Eigenvalues.cpp) Incompatible use of alias with dynamics");
  
  if (is_using_density_alias() && matrix_.get_use_squared_matrix())
    throw runtime_error("(In Eigenvalues.cpp) Incompatible use of alias with squared matrix (not implemented)");
  
  DFT_Vec v = v_in;
  if(vnormal_.size() == v.size()) v.IncrementBy_Scaled_Vector(vnormal_, -vnormal_.dotWith(v));
  
  if (verbose_)
  {
    cout << endl; cout << setprecision(12);
    cout << "\tIn Eigenvalues::matrix_dot_v:" << endl;
    cout << "\t  L2-Norm of v:   " << norm(v) << endl;
  }
  
  if (is_using_density_alias())
  {
    DFT_Vec vv = v;
    species_->convert_to_density_increment(vv);
    matrix_.matrix_dot_v1(vv,result,param);
    species_->convert_to_alias_deriv(result);
    
    DFT_Vec df; df.set(species_->getDF());
    DFT_Vec d2Rhodx2; species_->get_second_derivatives_of_density_wrt_alias(d2Rhodx2);
    df.Schur(df, d2Rhodx2);
    df.Schur(df, v);
    
    result.IncrementBy(df);
  }
  else if (matrix_.is_dynamic())
  {
    // Here we first multiply by the metric THEN by the matrix of second derivatives
    // The order of the operations is deliberate: In this generalized problem using
    // an arbitrary metric, we work with the transpose of the matrix for which we want
    // to solve the eigenproblem. The actual eigenvectors of the original matrix are g*v
    
    DFT_Vec gv(v); g_dot_v(v, gv);
    matrix_.matrix_dot_v1(gv, result, param, true);
    result.MultBy(-1); //Note: minus sign to be consistent with original metric g in DDFT object, which is not positive definite
    
    if (verbose_) cout << "\t  L2-Norm of g*v: " << norm(gv) << endl;
  }
  else
  {
    matrix_.matrix_dot_v1(v,result,param);
  }
  
  if (verbose_)
  {
    cout << "\t  L2-Norm of A*v: " << norm(result) << endl;
    cout << "\t  L2-Norm ratio |A*v|/|v|:     " << norm(result) / norm(v) << endl;
    cout << "\t  Dot product vT*A*v/|A*v||v|: " << v_dot_w(v,result) / norm(result) / norm(v) << endl;
  }
}


void Eigenvalues::matrix_dot_v(const vector<double> &vv, vector<double> &result, void *param)
{
  DFT_Vec v; v.zeros(vv.size());
  DFT_Vec r; r.zeros(vv.size());
  result.resize(vv.size());
  
  #ifdef USE_OMP
  #pragma omp parallel for
  #endif
  for (long i=0; i<vv.size(); i++) v.set(i, vv[i]);
  
  matrix_dot_v(v, r, param);
  
  #ifdef USE_OMP
  #pragma omp parallel for
  #endif
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
      
      eigen_vec_.MultBy(1/norm(eigen_vec_));
  }
  
  ////////////////////////////////////////
  // Minimisation of xAx/xx with FIRE2
  
  // Initialize
  DFT_Vec v; v.zeros(eigen_vec_.size());
  DFT_Vec df; df.zeros(eigen_vec_.size());
  DFT_Vec eigen_vec_old(eigen_vec_);
  
  double f = eigenvalues_objective_func(eigen_vec_, df, this);
  double f_old = f;
  double t = 0.0;
  
  // Parameters
  double alpha_start = 1.0;
  double alpha_min = 0.0;
  double alpha = alpha_start;
  double dt = 0.01/scale_;
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
    double P = -v_dot_w(v, df);
    
    double vnorm = norm(v);
    double fnorm = norm(df);
    
    double P_normalized = P/vnorm/fnorm;
    
    if (verbose_)
    {
      cout << endl; cout << setprecision(12);
      cout << "\tFire2 in Eigenvalues::calculate_eigenvector:" << endl;
      cout << "\t  P/|v||df| = " << P_normalized << endl;
      cout << "\t  dt = " << dt << endl;
      cout << "\t  alpha  = " << alpha << endl;
      cout << "\t  vnorm  = " << vnorm << endl;
      cout << "\t  fnorm  = " << fnorm << endl;
      cout << "\t  f-fold = " << f-f_old << endl;
    }
    
    if (P>=0) // && f<=f_old sometimes causes problems sending dt->0
    {
      Npos ++; Nneg = 0;
      
      /*
      if (Npos>Ndelay) // Original FIRE2
      {
        dt = (dt*finc<dt_max)?dt*finc:dt_max;
        alpha = (alpha*falf>alpha_min)?alpha*falf:alpha_min;
      }
      */
      
      if (Npos>Ndelay && P_normalized>0.999) // Modified FIRE2
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
      
      if (verbose_) cout << "\tBacktracking..." << endl;
      
      eigen_vec_.set(eigen_vec_old);
      v.zeros();
      
      f = eigenvalues_objective_func(eigen_vec_, df, this);
    }
    
    eigen_vec_old.set(eigen_vec_);
    
    // Semi-implicit Euler: Inertial velocity step
    v.IncrementBy_Scaled_Vector(df, -dt);
    
    vnorm = norm(v);
    fnorm = norm(df);
    
    // Semi-implicit Euler: velocity mixing and x-step
    v.MultBy(1-alpha);
    v.IncrementBy_Scaled_Vector(df, -alpha*vnorm/fnorm);
    eigen_vec_.IncrementBy_Scaled_Vector(v, dt);
    
    f_old = f;
    f = eigenvalues_objective_func(eigen_vec_, df, this);
    t += dt;
    
    // Check if converged
    if (!initialdelay || i>=Ndelay)
    {
      if (calculate_residual_error(false) < tol_) break;
    }
  }
  
  ////////////////////////////////////////
  // Finalise
  
  eigen_vec_.MultBy(1/norm(eigen_vec_));
  calculate_gradients(df);
}



