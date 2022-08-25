


#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <time.h>
#include <random>

#ifdef USE_OMP
#include <omp.h>
#endif

using namespace std;

#include "Dynamical_Matrix.h"
#include "myColor.h"

 
void Dynamical_Matrix::set_boundary_points_to_zero(DFT_Vec &v) const
{
  long pos = 0;
  
  do{
    v.set(pos,0.0);
  } while(get_next_boundary_point(pos));
  
}


void Dynamical_Matrix::matrix_dot_v1(const DFT_Vec &v, DFT_Vec &result, void *param) const
{
  vector<DFT_FFT> vwork(1);
  vwork[0].initialize(get_dimension(0), get_dimension(1), get_dimension(2));
  vwork[0].Real().set(v);

  if(is_fixed_boundary()) set_boundary_points_to_zero(vwork[0].Real());
  
  vwork[0].do_real_2_fourier();

  vector<DFT_Vec> rwork(1, result.size());

  matrix_dot_v(vwork,rwork,param);
  result.set(rwork[0]);

  if(is_fixed_boundary()) set_boundary_points_to_zero(result);
}
