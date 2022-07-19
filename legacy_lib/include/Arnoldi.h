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
  Arnoldi(const Dynamical_Matrix& matrix) : matrix_(matrix){}
  ~Arnoldi(){}

  double determine_unstable_eigenvector(vector<DFT_FFT> &eigen_vector, double shift, string Filename, int k, int p, long maxSteps, double tol) const;
  
protected:
  void matrix_dot_v(arma::cx_vec v, arma::cx_vec& d2F, double shift) const;
  bool check_factorisation(arma::cx_mat V, arma::cx_mat H, arma::cx_vec f, double shift,  double tol) const;
  bool check_eigenvectors(arma::cx_mat eigvec, arma::cx_vec eigval, double shift, double tol) const;
  void extend_arnoldi_factorisation(arma::cx_mat &V, arma::cx_mat &H, arma::cx_vec &f, const int k, const int p, double shift, double tol) const;

 private:
  const Dynamical_Matrix& matrix_;
};
#endif // __LUTSKO_DFT_ARNOLDI_ sentinal
