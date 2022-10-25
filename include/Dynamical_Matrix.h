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
  virtual void matrix_dot_v_intern(const vector<DFT_FFT> &v, vector<DFT_Vec> &result, void *param) const =0;
          void matrix_dot_v(const vector<DFT_FFT> &v, vector<DFT_Vec> &result, void *param) const;
          void matrix_dot_v1(const DFT_Vec &v, DFT_Vec &result, void *param) const;

  void set_boundary_points_to_zero(DFT_Vec &v) const;
    
  long get_Ntot() const { return get_dimension(0)*get_dimension(1)*get_dimension(2);}

  void set_verbose(bool v)            { dm_verbose_ = v;}
  void set_use_squared_matrix(bool v) { use_squared_matrix_  = v;}

protected:

  bool dm_verbose_         = false;
  bool use_squared_matrix_ = false;
};

#endif
