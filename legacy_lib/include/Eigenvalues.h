#ifndef __LUTSKO_DFT_EIGENVALUES__
#define __LUTSKO_DFT_EIGENVALUES__

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <complex>

#include <boost/serialization/vector.hpp>

#include "Dynamical_Matrix.h"
#include "Log.h"

class Eigenvalues
{
 public:
  Eigenvalues(const Dynamical_Matrix& matrix, bool verbose = false) : matrix_(matrix), verbose_(verbose){}
  ~Eigenvalues(){}

  void set_scale(double d)     {scale_ = d;}
  void set_tolerence(double d) {tol_ = d;}  
  void set_change_sign(bool b) {change_sign_ = b;}

  double get_scale()       const {return scale_;}
  double get_tolerence()   const {return tol_;}
  double get_eigenvalue()  const {return eigen_val_;}
  bool   get_change_sign() const {return change_sign_;}
  int    get_num_eval()    const {return num_eval_;}
  int    get_return_code() const {return return_code_;}
  
  void   calculate_eigenvector(Log &theLog);
  
  void   calculate_eigen_value();
  void   set_eigen_vec_(const vector<double> &v_in);  
  double calculate_gradients(DFT_Vec& df);
  
 private:
  const Dynamical_Matrix& matrix_;
  
  double scale_     = 1e6;
  double tol_       = 1e-6;
  bool change_sign_ = false;  
  int num_eval_     = 0;
  int return_code_  = 0; // nlopt return code
  bool verbose_     = false;

  DFT_Vec eigen_vec_;
  double  eigen_val_ = 0;
};

#endif //  __LUTSKO_DFT_EIGENVALUES__
