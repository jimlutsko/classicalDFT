#ifndef __LUTSKO_DFT_LOG_DET__
#define __LUTSKO_DFT_LOG_DET__

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <complex>

#include <boost/serialization/vector.hpp>

#include "Dynamical_Matrix.h"

// This is the method of Han, Maliouto, Avron and Shin from https://arxiv.org/pdf/1606.00942.pdf
// It is similar to the above method, but Chebychev interpolation is used since this is supposed to be more accurate.
// If do_condition = true, there is a zero eigenvalue associated with vector v0 = (1/sqrt(Ntot),...,1/sqrt(Ntot)) and so we use A->A+lam_max v0v0

class Log_Det
{
 public:
  Log_Det(const Dynamical_Matrix& matrix, int order, double lam_max, double lam_min, bool verbose = false);
  ~Log_Det(){}

  double calculate_log_det(long seed, int num_samples, bool has_zero_eigenvalue, double& variance);  

  void set_verbose(bool b) { verbose_ = b;}
  void set_debug(bool b)   { debug_   = b;}
  
protected:
  void get_coefficients(int order);
  void matrix_dot_v1(const DFT_Vec &v, DFT_Vec &result) const;
  
 private:
  const Dynamical_Matrix& matrix_;
  vector<double> c_;

  double lam_max_ = 0;
  double lam_min_ = 0;
  double a_       = 0;
  double b_       = 0;
  double scale_   = 0;
  
  bool verbose_ = false;
  bool debug_   = false;
};
#endif // __LUTSKO_DFT_LOG_DET_ sentinal
