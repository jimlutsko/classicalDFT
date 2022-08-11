#ifndef __LUTSKO_DYNAMICAL_MATRIX
#define __LUTSKO_DYNAMICAL_MATRIX

#include "DFT_LinAlg.h"

class Dynamical_Matrix
{
 public:
  Dynamical_Matrix(){}
  ~Dynamical_Matrix(){}

  double log_det_1(double lam_max, int num_samples, int order);
  double log_det_2(double lam_max, double lam_min, int num_samples, int order, bool do_condition = false);
  double log_fabs_det_2(double lam_max, double lam_min, int num_samples, int order);
  
    // Dynamical_Matrix interface
  virtual unsigned get_dimension(int direction) const =0;
  virtual long     get_Nboundary()              const =0;
  virtual bool     get_next_boundary_point(int &ix, int &iy, int &iz) const =0;
  virtual bool     get_next_boundary_point(long &pos) const =0;
  virtual long     boundary_pos_2_pos(int p)    const =0;
  virtual bool     is_boundary_point(long p) const =0;
  virtual bool     is_fixed_boundary() const =0;

  virtual void     matrix_dot_v(const vector<DFT_FFT> &v, vector<DFT_Vec> &result, void *param) const =0;

  void matrix_dot_v(const DFT_Vec &v, DFT_Vec &result, void *param) const;
  long get_Ntot() const { return get_dimension(0)*get_dimension(1)*get_dimension(2);}
};

#endif
