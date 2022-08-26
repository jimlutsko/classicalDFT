#ifndef __LUTSKO_DFT_ARNOLDI__
#define __LUTSKO_DFT_ARNOLDI__

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <complex>

#include <boost/serialization/vector.hpp>

#include <armadillo> // for arnoldi


#include "Dynamical_Matrix.h"

class Arnoldi
{
 public:
  Arnoldi(const Dynamical_Matrix& matrix, bool verbose = false) : matrix_(matrix), verbose_(verbose){}
  ~Arnoldi(){}
  
  // k = number of eigenvalues returned (ex. 10)
  // p = additional dimensions used  (ex. 15)
  double determine_largest_eigenvalue(vector<DFT_FFT> &eigen_vector, double shift, string Filename, int k, int p, long maxSteps = 1000000, double tol = 1e-8);

  void set_verbose(bool b) { verbose_ = b;}
  void set_debug(bool b)   { debug_   = b;}
  
  void set_initial_guess(DFT_FFT vector) 
  {
    provided_guess_ = true;
    initial_guess_ = vector;
    
    initial_guess_.Real().normalise();
	initial_guess_.do_real_2_fourier();
  }
  
  bool is_success()
  {
    if (!solver_executed_) throw runtime_error("Arnoldi: Solver has yet not been executed");
    else return solver_success_;
  }
  
  double get_eigenvalue(int i)
  {
    if (!solver_executed_) throw runtime_error("Arnoldi: Solver has yet not been executed");
    else return eigenvalues_[i];
  }
  
  DFT_FFT get_eigenvector(int i)
  {
    if (!solver_executed_) throw runtime_error("Arnoldi: Solver has yet not been executed");
    else return eigenvectors_[i];
  }
  
 protected:
  void matrix_dot_v(arma::cx_vec v, arma::cx_vec& d2F, double shift) const;
  bool check_factorisation(arma::cx_mat V, arma::cx_mat H, arma::cx_vec f, double shift,  double tol) const;
  bool check_eigenvectors(arma::cx_mat eigvec, arma::cx_vec eigval, double shift, double tol) const;
  void extend_arnoldi_factorisation(arma::cx_mat &V, arma::cx_mat &H, arma::cx_vec &f, const int k, const int p, double shift, double tol) const;

 private:
  const Dynamical_Matrix& matrix_;

  bool verbose_ = false;
  bool debug_   = false;
  
  bool provided_guess_ = false;
  DFT_FFT initial_guess_;
  
  bool solver_executed_ = false;
  bool solver_success_  = false;
  
  int k_;
  vector<double> eigenvalues_;
  vector<DFT_FFT> eigenvectors_;
};

#endif // __LUTSKO_DFT_ARNOLDI_ sentinal
