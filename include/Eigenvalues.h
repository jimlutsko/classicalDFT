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
#include "Species.h"
#include "Log.h"

class Eigenvalues
{
 public:
 
  Eigenvalues(const Dynamical_Matrix& matrix, Species* species, bool verbose = false) 
   : matrix_(matrix), species_(species), verbose_(verbose) {}
  
  ~Eigenvalues(){}

  void set_scale(double d)     {scale_ = d;}
  void set_verbose(bool v)     {verbose_ = v;}
  void set_tolerence(double d) {tol_ = d;}  
  void set_change_sign(bool b) {change_sign_ = b;}
  void set_vshift(DFT_Vec& v)  {vshift_ = v;}
  void set_max_num_eval(int i) {max_num_eval_ = i;}
  void set_use_density_alias(bool v) { using_density_alias_  = v;}
  
  double get_scale()        const {return scale_;}
  double get_tolerence()    const {return tol_;}
  double get_eigenvalue()   const {return eigen_val_;}
  bool   get_change_sign()  const {return change_sign_;}
  int    get_num_eval()     const {return num_eval_;}
  int    get_max_num_eval() const {return max_num_eval_;}
  int    get_return_code()  const {return return_code_;}
  const  DFT_Vec& get_eigen_vec() const {return eigen_vec_;}
  double get_eigen_val()          const { return eigen_val_;}
  Species *get_species()    const { return species_;}
  const Dynamical_Matrix& get_matrix() const {return matrix_;}
  
  DFT_Vec get_density_space_eigenvector() const;
  void set_from_density_space_eigenvector(const DFT_Vec &v_in); 
  
  void   reset_num_eval() {num_eval_ = 0;}
  void   calculate_eigenvector(Log &theLog);
  void   set_eigen_vec(const vector<double> &v_in);
  void   set_eigen_vec(const DFT_Vec &v_in) {eigen_vec_ = v_in;}
  void   save_snapshot() const;
  
  void g_dot_v(const DFT_Vec &v, DFT_Vec &gv) const;
  double v_dot_w(const DFT_Vec &v, const DFT_Vec &w) const;
  double norm(const DFT_Vec &v) const;
  
  double calculate_gradients(DFT_Vec& df);
  double calculate_residual_error(bool recompute_matrix_dot_v=true) const;
  
  void matrix_dot_v(const DFT_Vec &v, DFT_Vec &result, void *param) const;
  void matrix_dot_v(const vector<double> &vv, vector<double> &result, void *param);
  
  bool is_using_density_alias() const {return using_density_alias_;}
  //bool is_using_density_alias() const {return matrix_.is_using_density_alias();}
  
  void clear() { eigen_vec_.zeros(1); matrix_dot_eigen_vec_.zeros(1); }
  
 private:
  const Dynamical_Matrix& matrix_;
  Species* species_;
  DFT_Vec vshift_;
  
  double scale_     = 1;
  double tol_       = 1e-6;
  bool change_sign_ = false;  
  int num_eval_     = 0;
  int max_num_eval_ = -1;
  int return_code_  = 0; // nlopt return code
  bool verbose_     = false;
  bool using_density_alias_ = false;

  DFT_Vec eigen_vec_;
  DFT_Vec matrix_dot_eigen_vec_;
  double  eigen_val_ = 0;
};

#endif //  __LUTSKO_DFT_EIGENVALUES__
