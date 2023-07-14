#ifndef __LUTSKO_DYNAMICAL_MATRIX
#define __LUTSKO_DYNAMICAL_MATRIX

#include "DFT_LinAlg.h"

class Dynamical_Matrix
{
 public:
  Dynamical_Matrix(){}
  ~Dynamical_Matrix(){}
  
    // Dynamical_Matrix interface
  virtual unsigned get_dimension(int direction) const =0;
  virtual long     get_Nboundary()              const =0;
  virtual bool     get_next_boundary_point(int &ix, int &iy, int &iz) const =0;
  virtual bool     get_next_boundary_point(long &pos) const =0;
  virtual long     boundary_pos_2_pos(int p)    const =0;
  virtual bool     is_boundary_point(long p) const =0;
  virtual bool     is_fixed_boundary() const =0;
  virtual bool     is_dynamic() const = 0;

  // The first does A*v : the second will do A*A*v if use_squared_matrix_= true,
  // the third is a simpler interface for the second  
  virtual void matrix_dot_v_intern(const vector<DFT_FFT> &v, vector<DFT_Vec> &result, void *param, bool only_d2F) const =0;
  void matrix_dot_v(const vector<DFT_FFT> &v, vector<DFT_Vec> &result, void *param, bool only_d2F=false) const;
  void matrix_dot_v1(const DFT_Vec &v, DFT_Vec &result, void *param = NULL, bool only_d2F=false) const;
  
  virtual void get_matrix_diag(DFT_Vec &diag) const = 0;
  virtual void get_matrix_diag_nonhermetian(DFT_Vec &diag) const = 0;
  virtual void get_metric_diag(DFT_Vec &diag) const = 0;
  
  virtual void g_dot_x(const DFT_Vec& x, DFT_Vec& gx) const =0;

  void set_boundary_points_to_zero(DFT_Vec &v) const;
    
  long get_Ntot() const { return get_dimension(0)*get_dimension(1)*get_dimension(2);}

  void set_verbose(bool v)                 { dm_verbose_ = v;}
  void set_use_squared_matrix(bool v)      { use_squared_matrix_  = v;}
  void set_dm_in_alias_coordinates(bool v) { dm_in_alias_coordinates_ = v;}
  
  bool get_use_squared_matrix() const { return use_squared_matrix_;}
  bool is_matrix_in_alias_coordinates() const { return dm_in_alias_coordinates_;}

protected:

  bool dm_verbose_              = false;
  bool use_squared_matrix_      = false;
  bool dm_in_alias_coordinates_ = false; // Not compatible with DDFT dynamics (not implemented)
};

#endif
