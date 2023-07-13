#ifdef USE_SLEPC

#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>
#include <time.h>

using namespace std;

#ifdef USE_OMP
#include <omp.h>
#endif

#include <slepceps.h>
#include <slepcsys.h>
#include <petscerror.h>


#include "DFT_Slepc.h"

#include "myColor.h"

static char help[] = "Standard symmetric eigenproblem corresponding to the Laplacian operator in 1 dimension.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions = matrix dimension.\n\n";

PetscErrorCode DFT_Petsc::write_version_info(ostream &os)
{
  PetscFunctionBeginUser;
  char version[128];  
  //  PetscCall(PetscGetVersion(version,sizeof(version)));
  PetscGetVersion(version,sizeof(version));  
  os << version << endl;
  PetscFunctionReturn(PETSC_SUCCESS);
}

double bnorm = 0.0;
PetscErrorCode MyKSPMonitor(KSP ksp, PetscInt n, PetscReal rnorm, void *obj)
{
   PetscFunctionBeginUser;

   // In principle, this should be used if there are multiple processes in order to avoid jumbled output.
   //   PetscCall(PetscPrintf(PETSC_COMM_WORLD, "iteration %" PetscInt_FMT " KSP Residual norm %e \n", n, (double)rnorm));

   stringstream s;
   s << "\t\tSolution of g*x = dcluster: iteration " << n << " residual = " << (double) rnorm << " rel = " << rnorm/bnorm << "                         ";
   
   //   DFT_Petsc& pet = (*(DFT_Petsc*) obj);   
   //   pet.update_message(s.str());
   cout << '\r'; cout << myColor::YELLOW << s.str() << myColor::RESET;  cout.flush(); 

   PetscFunctionReturn(PETSC_SUCCESS);
 }

void DFT_Petsc::update_message(string s) 
{
  theLog_ << '\r';
  theLog_ << myColor::YELLOW << s << myColor::RESET;
  theLog_.flush();
}

void DFT_Petsc::DFT_to_SLEPC(const DFT_Vec &vin,  Vec &v) const
{
  PetscFunctionBegin;
  PetscScalar *pv;  
  PetscCallVoid(VecGetArray(v,&pv));    
  
  long p1 = 0;
  for(long p=0;p<vin.size();p++)
    if(!dm_.is_boundary_point(p))
      pv[p1++] = vin.get(p);

  PetscCallVoid(VecRestoreArray(v,&pv));
  PetscFunctionReturnVoid();
}

void DFT_Petsc::SLEPC_to_DFT(const  Vec &v, DFT_Vec &vout) const
{
  PetscFunctionBegin;  
  const PetscScalar *pv;
  PetscCallVoid(VecGetArrayRead(v,&pv));
  
  long p1 = 0;
  for(long p=0;p<vout.size();p++)
    vout.set(p,(dm_.is_boundary_point(p) ? 0.0 : pv[p1++]));

  PetscCallVoid(VecRestoreArrayRead(v,&pv));
  PetscFunctionReturnVoid();
}

PetscErrorCode MatMultA(Mat A,Vec x,Vec y)
{
  PetscFunctionBegin;
  
  void* ctx = NULL;
  PetscErrorCode ierr = MatShellGetContext(A,&ctx);CHKERRQ(ierr);

  DFT_Petsc *obj = (DFT_Petsc*) ctx;  
  PetscCall(obj->MatMultA1(x,y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DFT_Petsc::MatMultA1(Vec &x,Vec &y)
{
  PetscFunctionBegin;  

  DFT_Vec vx(Ntot_);
  DFT_Vec vy(Ntot_);

  SLEPC_to_DFT(x,vx);
  
  dm_.matrix_dot_v1(vx,vy,NULL);
  if(dm_.is_dynamic()) vy.MultBy(-1);

  DFT_to_SLEPC(vy,y);

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatMultTransposeA(Mat A,Vec x,Vec y)
{
  PetscFunctionBegin;
  
  void* ctx = NULL;
  PetscCall(MatShellGetContext(A,&ctx));
  DFT_Petsc *obj = (DFT_Petsc*) ctx;
  PetscCall(obj->MatMultTransposeA1(x,y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DFT_Petsc::MatMultTransposeA1(Vec &x,Vec &y)
{
  PetscFunctionBegin;
  
  DFT_Vec vx(Ntot_);
  DFT_Vec gx(Ntot_);
  DFT_Vec vy(Ntot_);

  SLEPC_to_DFT(x,vx);

  dm_.g_dot_x(vx,gx);  
  dm_.matrix_dot_v1(gx,vy,NULL,/*only_d2F=*/ true);
  if(dm_.is_dynamic()) vy.MultBy(-1);

  DFT_to_SLEPC(vy,y);

  PetscFunctionReturn(PETSC_SUCCESS);  
}

PetscErrorCode MatGetDiagonalA(Mat A,Vec diag)
{
  PetscFunctionBegin;
  
  cout << "MatGetDiagonal_A" << endl;

  void* ctx            = NULL;
  PetscCall(MatShellGetContext(A,&ctx));
  DFT_Petsc *obj = (DFT_Petsc*) ctx;
  PetscCall(obj->MatGetDiagonalA1(diag));

  PetscFunctionReturn(PETSC_SUCCESS);  
}

PetscErrorCode DFT_Petsc::MatGetDiagonalA1(Vec &diag)
{
  PetscFunctionBegin;
  
  DFT_Vec dft_diag(Ntot_); dft_diag.zeros();
  
  ifstream in("diag_out.dat");
  if(in.good())
    {
      in >> dft_diag;
      in.close();
    } else {
    dm_.get_matrix_diag(dft_diag);
    if(dm_.is_dynamic()) dft_diag.MultBy(-1);
    
    ofstream of("diag_out.dat");
    of << dft_diag;
    of.close();
  }
  
  DFT_to_SLEPC(dft_diag,diag);

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode MatMultG(Mat A,Vec x,Vec y)
{
  PetscFunctionBegin;
  
  void* ctx = NULL;
  PetscCall(MatShellGetContext(A,&ctx));

  DFT_Petsc *obj = (DFT_Petsc*) ctx;  
  PetscCall(obj->MatMultG1(x,y));

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode DFT_Petsc::MatMultG1(Vec &x,Vec &y)
{
  PetscFunctionBegin;  

  DFT_Vec vx(Ntot_);
  DFT_Vec vy(Ntot_);

  SLEPC_to_DFT(x,vx);
  
  dm_.g_dot_x(vx,vy);
  //  if(dm_.is_dynamic()) vy.MultBy(-1);

  DFT_to_SLEPC(vy,y);
 
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode MatGetDiagonalG(Mat A,Vec diag)
{
  PetscFunctionBegin;
  
  void* ctx            = NULL;
  PetscCall(MatShellGetContext(A,&ctx));
  DFT_Petsc *obj = (DFT_Petsc*) ctx;
  PetscCall(obj->MatGetDiagonalG1(diag));

  PetscFunctionReturn(PETSC_SUCCESS);  
}

PetscErrorCode DFT_Petsc::MatGetDiagonalG1(Vec &diag)
{
  PetscFunctionBegin;
  
  DFT_Vec dft_diag(Ntot_); dft_diag.zeros();  
  dm_.get_metric_diag(dft_diag);
  //  if(dm_.is_dynamic()) dft_diag.MultBy(-1);
      
  DFT_to_SLEPC(dft_diag,diag);

  PetscFunctionReturn(PETSC_SUCCESS);
}


DFT_Petsc::DFT_Petsc(Dynamical_Matrix &dm, Log &theLog, int argc, char**argv) : dm_(dm), theLog_(theLog)
{
  Ntot_ = dm_.get_Ntot();
  Ndynamic_ = Ntot_  - dm_.get_Nboundary();

  
  PetscCallVoid(PetscInitialize(&argc, &argv, (char *)0, help));
  
  
  // Make the basic objects
  // Create the operator matrix that defines the eigensystem, Ax=kx
  PetscCallVoid(MatCreateShell(PETSC_COMM_WORLD,Ndynamic_,Ndynamic_,Ndynamic_,Ndynamic_,this,&A_)); 
  PetscCallVoid(MatShellSetOperation(A_,MATOP_MULT,(void(*)(void))MatMultA));  
  PetscCallVoid(MatShellSetOperation(A_,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonalA));  
  PetscCallVoid(MatShellSetOperation(A_,MATOP_MULT_TRANSPOSE,(void(*)(void))MatMultTransposeA));
  
}

DFT_Petsc::~DFT_Petsc()
{
  PetscFunctionBegin;
  PetscCallVoid(MatDestroy(&A_));
  PetscFunctionReturnVoid();  
}  


PetscErrorCode DFT_Petsc::Solve_g_x_equals_b(DFT_Vec & dft_x, DFT_Vec &dft_b, stringstream &ret)
{
  PetscFunctionBegin;
  
  Vec x,b;
  PetscCall(VecCreate(PETSC_COMM_SELF, &x));
  PetscCall(PetscObjectSetName((PetscObject)x, "Solution"));
  PetscCall(VecSetSizes(x, Ndynamic_, Ndynamic_));

  PetscCall(VecSetFromOptions(x));
  
  PetscInt ss;
  PetscCall(VecGetSize(x,&ss));  
  PetscCall(VecDuplicate(x, &b));

  DFT_to_SLEPC(dft_b, b);  

  bnorm = dft_b.euclidean_norm();

  PetscCall(MatShellSetOperation(A_,MATOP_MULT,(void(*)(void))MatMultG));  
  PetscCall(MatShellSetOperation(A_,MATOP_GET_DIAGONAL,(void(*)(void))MatGetDiagonalG));  

  KSP ksp;    
  PetscCall(KSPCreate(PETSC_COMM_SELF, &ksp));
  PetscCall(KSPSetOperators(ksp, A_, A_));

  // Avoid using options
  
  PetscCall(KSPSetType(ksp, KSPBCGS));
  //  KSPSetType(ksp, KSPQMRCGS);
  PetscCall(KSPMonitorSet(ksp, MyKSPMonitor, NULL, 0));
  PetscCall(KSPSetTolerances(ksp, rel_tol_, abs_tol_, PETSC_DEFAULT, max_its_));

  PC pc;
  PetscCall(KSPGetPC(ksp, &pc));
  PetscCall(PCSetType(pc, PCJACOBI));

  //PetscCall(KSPSetFromOptions(ksp));

  ///////// SOLVE

  PetscCall(KSPSolve(ksp, b, x));

  SLEPC_to_DFT(x,dft_x);

  PetscInt its;
  PetscCall(KSPGetIterationNumber(ksp, &its));
  
  const char *strreason;
  PetscCall(KSPGetConvergedReasonString(ksp, &strreason));

  ret << "Num iterations = " << its << " converged because " << strreason;
  

  //  PetscCall(KSPView(ksp, PETSC_VIEWER_STDOUT_SELF));

  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&b));
  PetscCall(KSPDestroy(&ksp));

  PetscFunctionReturn(PETSC_SUCCESS);
  
}

#endif //USE_SLEPC
