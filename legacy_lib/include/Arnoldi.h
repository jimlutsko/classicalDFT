#ifndef __LUTSKO_DFT_ARNOLDI__
#define __LUTSKO_DFT_ARNOLDI__

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <complex>

#include <boost/serialization/vector.hpp>

#include "Density.h"

template <class T>
class Arnoldi
{
 public:
 Arnoldi(const Density& density, double (T::*matrix_dot_v)(double x) const) : matrix_dot_v_(matrix_dot_v), density_(density){}
  ~Arnoldi(){}

  void matrix_dot_v(arma::cx_vec v, arma::cx_vec& d2F, double shift, bool fixed_boundary, bool dynamic) const;

  //  double determine_unstable_eigenvector(vector<DFT_FFT> &eigen_vector, bool fixed_boundary, double shift, string Filename, bool dynamic, long maxSteps, double tol)const;

  double determine_unstable_eigenvector_IRArnoldi(vector<DFT_FFT> &eigen_vector, bool fixed_boundary, double shift, string Filename, bool dynamic, int k, int p, long maxSteps, double tol) const;
  //  void compute_and_sort_eigenvectors(arma::cx_mat H, arma::cx_mat &eigvec, arma::cx_vec &eigval);  
  bool check_factorisation(arma::cx_mat V, arma::cx_mat H, arma::cx_vec f, double shift, bool fixed_boundary, bool dynamic, double tol) const;
  bool check_eigenvectors(arma::cx_mat eigvec, arma::cx_vec eigval, double shift, bool fixed_boundary, bool dynamic, double tol) const;
  void extend_arnoldi_factorisation(arma::cx_mat &V, arma::cx_mat &H, arma::cx_vec &f, const int k, const int p, double shift, bool fixed_boundary, bool dynamic, double tol) const;

  //double determine_unstable_eigenvector_IRArnoldi(vector<DFT_FFT> &eigen_vector, bool fixed_boundary, double shift, string Filename, bool dynamic, int k, int p, long maxSteps, double tol) const;
 private:
    double (T::*matrix_dot_v_)(double x) const;
    Density &density_;
};
void save_Arnoldi_matrices(arma::cx_mat V, arma::cx_mat H);
void compute_and_sort_eigenvectors(arma::cx_mat H, arma::cx_mat &eigvec, arma::cx_vec &eigval);
#endif // __LUTSKO_DFT_ARNOLDI_ sentinal
