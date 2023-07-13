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


PetscErrorCode DFT_Slepc::write_version_info(ostream &os)
{
  PetscFunctionBeginUser;

  DFT_Petsc::write_version_info(os);

  char version[128];  
  PetscCall(SlepcGetVersion(version,sizeof(version)));  
  os << version << endl;  
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode DFT_Slepc::check()
{
  PetscFunctionBegin;
  
  DFT_Vec eigenr; eigenr.zeros(Ntot_);
  DFT_Vec eigeni; eigeni.zeros(Ntot_);
  double kr,ki;
  get_eigenvalue(kr,ki);
  get_eigenvector(eigenr, eigeni);
  check(eigenr, kr);

  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode DFT_Slepc::check(DFT_Vec &eigen, double kr)
{
  PetscFunctionBegin;
  
  theLog_ << endl;
  theLog_ << "CHECK whether we understand the results (by recalculating the eigenvalue, residuals, etc. by hand)" << endl << endl;
  
    // my check
  DFT_Vec Ax(Ntot_);
  dm_.matrix_dot_v1(eigen,Ax,NULL);  
  if(dm_.is_dynamic()) Ax.MultBy(-1);
  theLog_ << "xAx = " << eigen.dotWith(Ax) << endl;

  Ax.IncrementBy_Scaled_Vector(eigen,-kr);
  theLog_ << "||Ax-kx|| = " << Ax.euclidean_norm() << endl;
  theLog_ << "||Ax-kx||/||kx|| = " << Ax.euclidean_norm()/eigen.euclidean_norm()/abs(kr) << endl << endl;

  PetscFunctionReturn(PETSC_SUCCESS);
}


DFT_Slepc::DFT_Slepc(Dynamical_Matrix &dm, Log &theLog, int argc, char**argv) :
  DFT_Petsc(dm, theLog, argc, argv)
{
  PetscFunctionBegin;

  PetscCallVoid(SlepcInitialize(&argc, &argv,(char*)0,help));

  // Create eigensolver context
  PetscCallVoid(EPSCreate(PETSC_COMM_WORLD,&eps_));

  // Attach a monitor
  PetscViewerAndFormat *vf;
  PetscCallVoid(PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,&vf));
  PetscCallVoid(EPSMonitorSet(eps_,(PetscErrorCode (*)(EPS,PetscInt,PetscInt,PetscScalar*,PetscScalar*,PetscReal*,PetscInt,void*))EPSMonitorFirst,vf,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy));

  PetscFunctionReturnVoid();
}


DFT_Slepc::~DFT_Slepc()
{
  PetscFunctionBegin;
  PetscCallVoid(EPSDestroy(&eps_));
  
  //Free work space
  PetscCallVoid(SlepcFinalize());
  PetscFunctionReturnVoid();  
}  


int DFT_Slepc::run_eigenproblem(int argc, char **argv)
{
  PetscFunctionBegin;
  
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"\nDFT Eigenproblem, n=%ld two_sided = %d\n\n",Ndynamic_,two_sided_));
  
  // Define problem and set operators: Hermetian for static problem ... and we use JD
  if(dm_.is_dynamic() == false)
    {
      PetscCall(EPSSetOperators(eps_,A_,NULL));
      PetscCall(EPSSetProblemType(eps_,EPS_HEP));
      PetscCall(EPSSetType(eps_, EPSJD));

      ST st;
      PetscCall(EPSGetST(eps_, &st));
      KSP ksp;
      PetscCall(STGetKSP(st,&ksp));
      PetscCall(KSPSetType(ksp,KSPGMRES));
      // Note that changing maxits to a larger value (up to about 40) sometimes helps
      PetscCall(KSPSetTolerances(ksp, /*rtol =*/min(0.1*eps_tol_,0.99), /*abstol=*/PETSC_DEFAULT, /*dtol=*/PETSC_DEFAULT, /*maxits = */ 10));
    
      PC pc;
      PetscCall(KSPGetPC(ksp, &pc));
      PetscCall(PCSetType(pc, PCJACOBI));
  } else {
    // Non=hermetian with Krylov-Schur
    PetscCall(EPSSetOperators(eps_,A_,NULL));
    PetscCall(EPSSetProblemType(eps_,EPS_NHEP));
    PetscCall(EPSSetType(eps_, EPSKRYLOVSCHUR));
    
    ST st;
    PetscCall(EPSGetST(eps_, &st));
    KSP ksp;
    PetscCall(STGetKSP(st,&ksp));
    PetscCall(KSPSetType(ksp,KSPGMRES)); 
    // Somewhere, it says the tolerance here should be lower than that demanded for the eigenvalues
    PetscCall(KSPSetTolerances(ksp, /*rtol =*/min(0.1*eps_tol_,0.99), /*abstol=*/PETSC_DEFAULT, /*dtol=*/PETSC_DEFAULT, /*maxits = */PETSC_DEFAULT));

    PC pc;
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(PCSetType(pc, PCJACOBI));
    
    // This appears to be very important ...
    PetscCall(EPSSetBalance(eps_,EPS_BALANCE_TWOSIDE,/*PetscInt its=*/PETSC_DEFAULT,/*PetscReal cutoff=*/ PETSC_DEFAULT));    

  }

  // Should get this working someday so I can eliminate the messy argc/argv stuff
  //PetscCall(EPSMonitorSet(eps, (PetscErrorCode (*)(EPS, PetscInt, PetscInt, PetscScalar*, PetscScalar*, PetscReal*, PetscInt, void*)) EPSMonitorFirst, NULL, NULL));

  //Define problem (smallest real eigevalue) and convergence criterion
  PetscCall(EPSSetDimensions(eps_,num_eigenvalues_,PETSC_DEFAULT,PETSC_DEFAULT));      
  if(smallest_) PetscCall(EPSSetWhichEigenpairs(eps_,EPS_SMALLEST_REAL));
  else PetscCall(EPSSetWhichEigenpairs(eps_,EPS_LARGEST_REAL));
  PetscCall(EPSSetConvergenceTest(eps_,EPS_CONV_REL));
  PetscCall(EPSSetTolerances(eps_, eps_tol_, PETSC_DEFAULT));
  PetscCall(EPSSetTwoSided(eps_,(two_sided_ ? PETSC_TRUE : PETSC_FALSE)));
  
  // Set solver parameters at runtime
  // This is good for playing but I suppress it for now.
  PetscCall(EPSSetFromOptions(eps_));

  PetscCall(EPSSetBalance(eps_,EPS_BALANCE_TWOSIDE,/*PetscInt its=*/PETSC_DEFAULT,/*PetscReal cutoff=*/ PETSC_DEFAULT));    
  
  // Initial guess
  //if(input_vector_dft_.empty() == false && input_vector_dft_.find_first_not_of(' ') != std::string::npos)
  //  read_input_vec();
  if(input_vector_dft_.size()>1) init_from_dft_vector();
  else init_with_random_vector();

  // Solve
  PetscCall(EPSSolve(eps_));

  PetscCall(DisplaySolution(eps_));
  
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode DFT_Slepc::compute_error(double &err, int which) const
{
  PetscFunctionBegin;
  EPSComputeError(eps_, which, EPS_ERROR_RELATIVE, &err);
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode DFT_Slepc::get_number_converged(int &nconv) const
{
  PetscFunctionBegin;  
  EPSGetConverged(eps_,&nconv);
  PetscFunctionReturn(PETSC_SUCCESS); 
}


PetscErrorCode DFT_Slepc::DisplaySolution(EPS &eps)
{
  PetscFunctionBegin;
  
  // Optional: Get some information from the solver and display it
  EPSType      type;
  PetscInt     nev;
  
  PetscCall(EPSGetType(eps_,&type));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type));
  PetscCall(EPSGetDimensions(eps_,&nev,NULL,NULL));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %d\n",nev));
  
  // Display solution and clean up
  // show detailed info unless -terse option is given by user
  
  PetscBool terse;  

  PetscCall(PetscOptionsHasName(NULL,NULL,"-terse",&terse));
  if (terse) {
    PetscCall(EPSErrorView(eps_,EPS_ERROR_RELATIVE,NULL));
  } else {
    PetscCall(PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_ASCII_INFO_DETAIL));
    PetscCall(EPSConvergedReasonView(eps_,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(EPSErrorView(eps_,EPS_ERROR_RELATIVE,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(EPSErrorView(eps_,EPS_ERROR_ABSOLUTE,PETSC_VIEWER_STDOUT_WORLD));
    PetscCall(PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}


void DFT_Slepc::set_input_vec(string filename)
{
  theLog_ << "reading DFT input vector from " << filename;
  ifstream in(filename.c_str());
  in >> input_vector_dft_;
  in.close();
}

void DFT_Slepc::set_input_vec(const DFT_Vec &v)
{
  input_vector_dft_.set(v);
}

PetscErrorCode DFT_Slepc::init_from_dft_vector()
{
  PetscFunctionBegin;
  
  Vec v0;
  PetscCall(MatCreateVecs(A_,&v0,NULL));
  DFT_to_SLEPC(input_vector_dft_,v0);
  PetscCall(EPSSetInitialSpace(eps_,1,&v0));

  theLog_  << endl;

  DFT_Vec Ax(Ntot_);
  dm_.matrix_dot_v1(input_vector_dft_,Ax,NULL,/*only_d2F=*/ false);
  if(dm_.is_dynamic()) Ax.MultBy(-1);
  theLog_ << "xAx = " << Ax.dotWith(input_vector_dft_) << endl;
  
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode DFT_Slepc::init_with_random_vector()
{
  PetscFunctionBegin;
  
  PetscRandom rctx;
  PetscRandomCreate(PETSC_COMM_WORLD,&rctx);
  Vec v0;
  PetscCall(MatCreateVecs(A_,&v0,NULL));
  VecSetRandom(v0,rctx);
  PetscCall(EPSSetInitialSpace(eps_,1,&v0));
  PetscRandomDestroy(&rctx);
  
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode DFT_Slepc::get_eigenvalue(double &eigenvalue, double &eigenvalue_imag, int which) const
{
  PetscFunctionBegin;
  PetscCall(EPSGetEigenvalue(eps_,which, &eigenvalue, &eigenvalue_imag));
  PetscFunctionReturn(PETSC_SUCCESS);  
}


PetscErrorCode DFT_Slepc::get_eigenvector(DFT_Vec &eigenvec, DFT_Vec &eigenvec_imag, int which) const
{
  PetscFunctionBegin;
  
  Vec xr,xi;
  PetscCall(MatCreateVecs(A_,NULL,&xr));
  PetscCall(MatCreateVecs(A_,NULL,&xi));
  
  PetscCall(EPSGetEigenvector(eps_,which, xr, xi));

  SLEPC_to_DFT(xr,eigenvec);
  SLEPC_to_DFT(xi,eigenvec_imag);
  
  PetscFunctionReturn(PETSC_SUCCESS);  
}


PetscErrorCode DFT_Slepc::get_eigenvector_left(DFT_Vec &eigenvec, DFT_Vec &eigenvec_imag, int which) const
{
  PetscFunctionBegin;
  
  Vec xr,xi;
  PetscCall(MatCreateVecs(A_,NULL,&xr));
  PetscCall(MatCreateVecs(A_,NULL,&xi));
  
  PetscCall(EPSGetLeftEigenvector(eps_,which, xr, xi));

  SLEPC_to_DFT(xr,eigenvec);
  SLEPC_to_DFT(xi,eigenvec_imag);
  
  PetscFunctionReturn(PETSC_SUCCESS);  
}


PetscErrorCode DFT_Slepc::write_output_vectors_dft()
{  
  PetscFunctionBegin;

  ofstream of("eigenvalues.dat");
  
  PetscInt nconv = 0;
  PetscCall(EPSGetConverged(eps_, &nconv));  
  of << nconv << endl;

  double kr,ki;
  for (PetscInt i=0;i<nconv;i++)
    {
      get_eigenvalue(kr, ki, i);
      of << kr << " " << ki << endl;
      theLog_  << kr << " " << ki << endl;
    }  

  of.close();
  
  DFT_Vec eigenr(Ntot_);
  DFT_Vec eigeni(Ntot_);  
  for (PetscInt i=0;i<nconv;i++)
    {            
      get_eigenvector(eigenr, eigeni, i);

      stringstream ss;
      ss << "eigenvector_" << i << ".dat";
      ofstream of(ss.str().c_str());      
      of << eigenr << eigeni;
      of.close();
    }

  if(two_sided_)
    for (PetscInt i=0;i<nconv;i++)
      {            
	get_eigenvector_left(eigenr, eigeni, i);
	
	stringstream ss;
	ss << "eigenvector_left_" << i << ".dat";
	ofstream of(ss.str().c_str());      
	of << eigenr << eigeni;
	of.close();
      }


  theLog_ << "Wrote " << nconv << " eigenvals/eigenvecs " << endl;	  	  
  
  PetscFunctionReturn(PETSC_SUCCESS);
}

#endif //USE_SLEPC
