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

class Arnoldi
{
 public:
  Arnoldi(const DFT_Matrix& matrix, bool fixed_boundary = false) : matrix_(matrix), fixed_boundary_(fixed_boundary){}
  ~Arnoldi(){}

  void matrix_dot_v(arma::cx_vec v, arma::cx_vec& d2F, double shift) const;

  double determine_unstable_eigenvector_IRArnoldi(vector<DFT_FFT> &eigen_vector, double shift, string Filename, int k, int p, long maxSteps, double tol) const;

  bool check_factorisation(arma::cx_mat V, arma::cx_mat H, arma::cx_vec f, double shift,  double tol) const;
  bool check_eigenvectors(arma::cx_mat eigvec, arma::cx_vec eigval, double shift, double tol) const;
  void extend_arnoldi_factorisation(arma::cx_mat &V, arma::cx_mat &H, arma::cx_vec &f, const int k, const int p, double shift, double tol) const;

 private:
  const DFT_Matrix& matrix_;
  bool fixed_boundary_ = false;
};
#endif // __LUTSKO_DFT_ARNOLDI_ sentinal
