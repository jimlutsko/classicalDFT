#ifdef USE_SLEPC

#include <slepceps.h>
#include <slepcsys.h>
#include <petscerror.h>

#include "Dynamical_Matrix.h"
#include "Log.h"

class DFT_Petsc
{
 public:
  DFT_Petsc(Dynamical_Matrix &dm, Log &theLog, int argc, char**argv);
  ~DFT_Petsc();

  PetscErrorCode Solve_g_x_equals_b(DFT_Vec & dft_x, DFT_Vec &dft_b, stringstream &ret);
  virtual PetscErrorCode write_version_info(ostream &os);

  void set_abs_tol(double x) { abs_tol_ = x;}
  void set_rel_tol(double x) { rel_tol_ = x;}
  void set_max_its(int x)    { max_its_ = x;}
  
  void update_message(string s);
  
  PetscErrorCode MatMultA1(Vec &x,Vec &y);
  PetscErrorCode MatGetDiagonalA1(Vec &diag);
  PetscErrorCode MatMultTransposeA1(Vec &x,Vec &y);

  PetscErrorCode MatMultG1(Vec &x,Vec &y);
  PetscErrorCode MatGetDiagonalG1(Vec &diag);
  
  void DFT_to_SLEPC(const DFT_Vec &vin,  Vec& v) const;
  void SLEPC_to_DFT(const  Vec &v, DFT_Vec &vout) const;

 protected:
  Dynamical_Matrix &dm_;
  Log &theLog_;
  long Ntot_      = 0;
  long Ndynamic_  = 0;

  double abs_tol_ = 1e-5;
  double rel_tol_ = 1e-7;
  int    max_its_ = 100000;

  Mat A_;
};


class DFT_Slepc : public DFT_Petsc
{
 public:
  DFT_Slepc(Dynamical_Matrix &dm, Log &theLog, int argc, char**argv);
  ~DFT_Slepc();

  virtual PetscErrorCode write_version_info(ostream &os);
  
  int run_eigenproblem(int argc, char**argv);
  
  void set_input_vec(string  filename);
  void set_input_vec(const DFT_Vec &v);

  void set_eps_tol(double v) { eps_tol_ = v;}
  void set_num_eig(int n) { num_eigenvalues_ = n;}
  void set_two_sided(bool two_sided) { two_sided_ = two_sided;}
  void set_smallest(bool s) { smallest_ = s;}
  
  PetscErrorCode MatMultA1(Vec &x,Vec &y);
  PetscErrorCode MatGetDiagonalA1(Vec &diag);
  PetscErrorCode MatMultTransposeA1(Vec &x,Vec &y);

  PetscErrorCode check();
  PetscErrorCode check(DFT_Vec &eigen, double kr);
  PetscErrorCode DisplaySolution(EPS &eps);
  
  PetscErrorCode init_from_dft_vector();
  PetscErrorCode init_with_random_vector();
  
  PetscErrorCode write_output_vectors_dft();

  PetscErrorCode get_eigenvalue(double &eigenvalue, double &eigenvalue_imag, int which = 0) const;
  PetscErrorCode get_eigenvector(DFT_Vec &eigenvec, DFT_Vec &eigenvector_imag, int which = 0) const;
  PetscErrorCode get_eigenvector_left(DFT_Vec &eigenvec, DFT_Vec &eigenvector_imag, int which = 0) const;
  PetscErrorCode get_number_converged(int &nconv) const;
  
  PetscErrorCode compute_error(double &err, int which = 0) const;
  
  double  get_eigenvalue (int which = 0) const {double  kr,ki; get_eigenvalue(kr,ki,which);  return kr;}
  DFT_Vec get_eigenvector(int which = 0) const {DFT_Vec vr(Ntot_),vi(Ntot_); get_eigenvector(vr,vi,which); return vr;}
  double  get_res_error(int which = 0) const {double err; compute_error(err, which); return err;}
  
 protected:
  DFT_Vec input_vector_dft_;

  double eps_tol_      = 1e-4; // relative tolerance
  int num_eigenvalues_ = 1;
  bool two_sided_      = false;
  bool smallest_       = true;
  
  EPS eps_;
};

#endif //USE_SLEPC
