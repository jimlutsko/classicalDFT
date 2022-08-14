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

#include <complex>

#include "myColor.h"

#include "Arnoldi.h"

static void save_Arnoldi_matrices(arma::cx_mat V, arma::cx_mat H);
static void compute_and_sort_eigenvectors(arma::cx_mat H, arma::cx_mat &eigvec, arma::cx_vec &eigval);

void Arnoldi::matrix_dot_v(arma::cx_vec v, arma::cx_vec& d2F, double shift) const
{
	const int species = 0;
	//	const Density& density = density_; //dft_->getDensity(species);
	const int Nx = matrix_.get_dimension(0);
	const int Ny = matrix_.get_dimension(1);
	const int Nz = matrix_.get_dimension(2);
	const long Ntot = matrix_.get_Ntot();
		
	//real part
	vector<DFT_FFT> dft_v(1); dft_v[0].initialize(Nx,Ny,Nz);
	for(long i=0; i<Ntot; i++)
	  {
	    if(matrix_.is_fixed_boundary() && matrix_.is_boundary_point(i)) dft_v[0].Real().set(i,0.0);
	    else dft_v[0].Real().set(i,v[i].real());
	  }
	    
	dft_v[0].do_real_2_fourier();
	
	//imag part
	vector<DFT_FFT> dft_w(1); dft_w[0].initialize(Nx,Ny,Nz);
	for(long i=0; i<Ntot; i++)
	  {
	    if(matrix_.is_fixed_boundary() && matrix_.is_boundary_point(i)) dft_w[0].Real().set(i,0.0);
	    else dft_w[0].Real().set(i,v[i].imag());
	  }
	dft_w[0].do_real_2_fourier();
	
	//real part
	vector<DFT_Vec> dft_d2F(1); dft_d2F[0].zeros(Ntot);
	matrix_.matrix_dot_v(dft_v, dft_d2F, NULL); // JFL , fixed_boundary, dynamic);
	
	//imag part
	vector<DFT_Vec> dft_d2F_imag(1); dft_d2F_imag[0].zeros(Ntot);
	matrix_.matrix_dot_v(dft_w, dft_d2F_imag, NULL); //JFL , fixed_boundary, dynamic);
	
	d2F.zeros(Ntot);
	for (long i=0; i<Ntot; i++) d2F[i] = dft_d2F[0].get(i) + std::complex<double>(0,1)*dft_d2F_imag[0].get(i) + shift*v[i];

	// This is the responsibility of the matrix_.matrix_dot_v function ...
	if (matrix_.is_fixed_boundary())
	  for (long i=0; i<matrix_.get_Nboundary(); i++)
	    d2F[matrix_.boundary_pos_2_pos(i)] = 0.0;
}



double Arnoldi::determine_largest_eigenvalue(vector<DFT_FFT> &eigen_vector, double shift, string Filename, int k, int p, long maxSteps, double tol) const
{
  if(verbose_) cout << endl;
  if(verbose_) cout << myColor::YELLOW;
  if(verbose_) cout << "\tFixed boundary = " << matrix_.is_fixed_boundary() << ", MaxIterations = " << maxSteps << ", tolerence = " << tol << endl;
  if(verbose_) cout << myColor::RESET;
  if(verbose_) cout << endl;
	
  int sysres = system("zip -r -q arnoldi_backup.zip arnoldi/ eigenvectors/");
  sysres = system("rm -r arnoldi");
  sysres = system("mkdir arnoldi");
  sysres = system("rm -r eigenvectors");
  sysres = system("mkdir eigenvectors");
	
  ofstream ofile_iter("arnoldi/iterations.dat");
  ofile_iter << "# Largest eigenvalues from Implicitely Restarted Arnoldi method" << endl;
  ofile_iter << "# " << endl;
  ofile_iter << "# These are Ritz values (real and imaginary components) and Ritz" << endl;
  ofile_iter << "# estimates for the errors. In other words these are Rayleigh quotients" << endl;
  ofile_iter << "# using the best approximations of the associated eigenvectors and the" << endl;
  ofile_iter << "# error is equivalent to the residual |Av-xv|." << endl;
  ofile_iter << "# " << endl;
  ofile_iter << "#" << setw(7) << "iter*p" << setw(16) <<  "real"  << setw(16) <<  "imag"  << setw(16) << "error" << endl;
	
  const int species = 0;
  const int Nx = matrix_.get_dimension(0);
  const int Ny = matrix_.get_dimension(1);
  const int Nz = matrix_.get_dimension(2);
  const long Ntot = Nx*Ny*Nz;	
	
  // Pass negative tolerance to tell the algorithm it must keep iterating
  bool dont_stop_iterating = false;
  if (tol<0) {dont_stop_iterating = true; tol = abs(tol);}
	
  arma::cx_mat Vk,Hk;
  arma::cx_vec fk;
	
  extend_arnoldi_factorisation(Vk, Hk, fk, 0, k, shift, tol);
	
  arma::cx_vec Hk_eigval;
  arma::cx_mat Hk_eigvec;
	
  compute_and_sort_eigenvectors(Hk, Hk_eigvec, Hk_eigval);
	
  // Iterate
  int iter = 0;
  bool converged = false;
  double eigen_value_old = 0.0;
	
  arma::cx_vec eigval;
  arma::cx_mat eigvec;
	
  while (!converged || dont_stop_iterating)
    {
      // p-step extension of Arnoldi factorisation
		
      arma::cx_mat Vp = Vk;
      arma::cx_mat Hp = Hk;
      arma::cx_vec fp = fk;
		
      extend_arnoldi_factorisation(Vp, Hp, fp, k, p, shift, tol);
		
      arma::cx_vec Hp_eigval;
      arma::cx_mat Hp_eigvec;
		
      compute_and_sort_eigenvectors(Hp, Hp_eigvec, Hp_eigval);
		
      // Sort eigenvalues for extended H-matrix
		
      arma::cx_vec Hk_eigval(k);
      arma::cx_vec Hp_eigval_discard(p);
		
      for (int i=0; i<k; i++) Hk_eigval[i] = Hp_eigval[i];
      for (int i=k; i<k+p; i++) Hp_eigval_discard[i] = Hp_eigval[i];
		
      // Filter out unwanted eigenvalues
		
      arma::cx_vec q(k+p); q.zeros(); q(k+p-1) = 1;
		
      for (int i=0; i<p; i++)
	{
	  arma::cx_mat Hp_shifted = Hp;
	  for (int j=0; j<k+p; j++) Hp_shifted(j,j) -= Hp_eigval_discard[i];
			
	  arma::cx_mat Q,R;
	  arma::qr(Q,R,Hp_shifted);
			
	  Hp = Q.t()*Hp*Q;
	  Vp = Vp*Q;
			
	  Hp_eigvec = Q.t()*Hp_eigvec; // Missing information in Radke's thesis
			
	  q = Q.t()*q; // Not the original line: misktake in Radke's thesis
	}
		
      // Update Arnoldi factorisation
      long nelem = Vk.n_rows;
      for (long i=0; i<nelem; i++) for (int j=0; j<k; j++) Vk(i,j) = Vp(i,j);
      for (int i=0; i<k; i++) for (int j=0; j<k; j++) Hk(i,j) = Hp(i,j);
      for (int i=0; i<k; i++) for (int j=0; j<k; j++) Hk_eigvec(j,i) = Hp_eigvec(j,i);
		
      arma::cx_vec vv(nelem); for (long i=0; i<nelem; i++) vv[i] = Vp(i,k);
      fk = vv*Hp(k,k-1) + fp*q(k-1);
		
      save_Arnoldi_matrices(Vk,Hk);
		
      if (!check_factorisation(Vk,Hk,fk,shift,tol))
	throw runtime_error("(IRA) Arnoldi factorisation does not check out");
		
      // Update approximations for eigenvalues/vectors
		
      eigval = Hk_eigval;
      eigvec = arma::cx_mat(nelem,k);
		
      for (int i=0; i<k; i++)
	{
	  arma::cx_vec wi(k);
	  for (int j=0; j<k; j++) wi[j] = Hk_eigvec(j,i);
			
	  arma::cx_vec vi = Vk*wi;
	  for (long j=0; j<nelem; j++) eigvec(j,i) = vi[j];
	}
		
      // TODO: They are not exactly the same, why??
      // Note: For testing only -- The norm of the residuals are the same as the Ritz estimates
      check_eigenvectors(eigvec,eigval,shift,tol);
		
      // Convergence checks
		
      iter++;
      //JFL
      //		if (iter>=maxSteps) throw runtime_error("IRArnoldi method: Exceeded max iterations");
      //		converged = true;
      if (iter>=maxSteps) 	return eigval[0].real()-shift;
      
      for (int i=0; i<k; i++)
	{
	  double ritz_estimate = abs(Hk_eigvec(k-1,i))*norm(fk);
			
	  // Not converged if (large ritz estimate) and 
	  // (eigenvalue not degenerated with the largest discarded value)
	  // The second condition is necessary because if the last eigenvalue
	  // is degenerate with another that is discarded, we never complete
	  // the calculation of the corresponding vector subspace since we
	  // never record all necessary Arnoldi vectors to contruct the basis.
			
	  if ( ritz_estimate > tol*abs(eigval[i]) &&
	       abs(eigval[i]-Hp_eigval_discard[0]) > tol) 
	    {
	      converged = false;
	    }
			
	  ofile_iter << scientific << setprecision(6);
	  ofile_iter << setw(8) << iter*p << setw(16) << eigval[i].real()-shift << setw(16) << eigval[i].imag() << setw(16) << ritz_estimate << endl;
	}
      ofile_iter << endl;
		
      // Save leading eigenvectors
      // JFL: This code needs to be replaced by something that only uses native (Armadillo) functions
		
      eigen_vector[species].initialize(Nx,Ny,Nz);
		
      for (int i=k-1; i>=0; i--)
	{
	  for(long j=0; j<Ntot; j++)
	    eigen_vector[species].Real().set(j,eigvec(j,i).real());
			
	  eigen_vector[species].Real().normalise();
	  eigen_vector[species].do_real_2_fourier();
			
	  ofstream of("eigenvectors/density_eigenvector_"+to_string(i)+".dat");
	  of << eigen_vector[species].Real();
	  of.close();
	}
	       
      // Report in terminal
		
      cout << myColor::YELLOW;
      cout << setprecision(6);
      cout << '\r'; cout << "\t" << "iteration = " << iter << " shift = " << shift << " eigen_value = " << setw(12) << eigval[0].real()-shift;
      cout << myColor::RESET;
		
      ofstream debug("debug.dat", (iter == 0 ? ios::trunc : ios::app));
      debug << iter << " " << eigval[0].real()-shift << " " << fabs(eigval[0].real()-shift - eigen_value_old)/fabs(eigval[0].real()-shift) << " " << fabs(eigval[0].real()-shift - eigen_value_old) << " " << abs(Hk_eigvec(k-1,0))*norm(fk) << endl;
      debug.close();
		
      eigen_value_old = eigval[0].real()-shift;
    }
  cout << endl;
	
  // This is not already initialized ???
  eigen_vector[species].initialize(Nx,Ny,Nz);		
  for(long j=0; j<Ntot; j++)
    eigen_vector[species].Real().set(j,eigvec(j,0).real()); 	    
  eigen_vector[species].Real().normalise();
  eigen_vector[species].do_real_2_fourier();
	
  if(verbose_) cout << endl;
	
  return eigval[0].real()-shift;
}

void compute_and_sort_eigenvectors(arma::cx_mat H, arma::cx_mat &eigvec, arma::cx_vec &eigval)
{
	if (!H.is_square()) throw runtime_error("(DDFT_IF.cpp) Attempting to compute eigenvectors of non-square matrix");
	int k = H.n_cols;
	
	eigvec = arma::cx_mat(k,k);
	eigval = arma::cx_vec(k);
	
	arma::cx_mat eigvec_raw;
	arma::cx_vec eigval_raw;
	arma::eig_gen(eigval_raw, eigvec_raw, H);
	
	arma::vec eigval_real(k);
	for (int i=0; i<k; i++) eigval_real[i] = abs(eigval_raw[i]); //eigval_raw[j].real();
	arma::uvec indices = sort_index(eigval_real, "descend");
	
	for (int i=0; i<k; i++) eigval[i] = eigval_raw[indices[i]];
	for (int i=0; i<k; i++) for (int j=0; j<k; j++) eigvec(j,i) = eigvec_raw(j,indices[i]);
}



bool Arnoldi::check_factorisation(arma::cx_mat V, arma::cx_mat H, arma::cx_vec f, double shift, double tol) const
{
	if (V.n_cols!=H.n_rows || !H.is_square()) throw runtime_error("Inconsistent dimensions in Arnoldi factorisation");
	int k = V.n_cols;
	long Ntot = V.n_rows;
	
	arma::cx_mat R(Ntot, k); //= A*V - V*H - fk*ek;
	for (int j=0; j<k; j++)
	{
		arma::cx_vec v(Ntot); for (long i=0; i<Ntot; i++) v[i] = V(i,j);
		arma::cx_vec w(Ntot); matrix_dot_v(v, w, shift);
		
		for (long i=0; i<Ntot; i++)
		{
			R(i,j) = w[i];
			for (int l=0; l<k; l++) R(i,j) -= V(i,l)*H(l,j);
			if (j==k-1) R(i,j) -= f[i];
		}
	}
	
	double error = arma::norm(R,"fro");
	
	ofstream ofile_check("arnoldi/check_factorisation.dat", ios::app);
	
	ofile_check << scientific << setprecision(2);
	ofile_check << "Norm of Residual:     " << error << endl;
	ofile_check << "While tol*sqrt(N*k):  " << tol*sqrt(Ntot*k) << endl;
	ofile_check << endl;
	
	return (error<tol*sqrt(Ntot*k));
}


bool Arnoldi::check_eigenvectors(arma::cx_mat eigvec, arma::cx_vec eigval, double shift, double tol) const
{
	const long Ntot = matrix_.get_Ntot();
	
	if (eigvec.n_cols!=eigval.n_rows) throw runtime_error("Inconsistent dimensions in eigvec/eigval");
	int k = eigval.n_rows;
	
	if (matrix_.is_fixed_boundary())
	  for (long i=0; i<matrix_.get_Nboundary(); i++)
	    for (int j=0; j<k; j++)
	      {
		if (abs(eigvec(matrix_.boundary_pos_2_pos(i),j))>tol) 
		  throw runtime_error("Eigenvector has non-zero values on boundaries");
	      }
	
	ofstream ofile_check("arnoldi/check_eigenvectors.dat", ios::app);
	ofile_check << scientific << setprecision(2);
	
	double err_max = 0.0;
	
	for (int j=0; j<k; j++)
	{
		arma::cx_vec v(Ntot); for (long i=0; i<Ntot; i++) v[i] = eigvec(i,j);
		arma::cx_vec w(Ntot); matrix_dot_v(v, w, shift);
		for (long i=0; i<Ntot; i++) w[i] -= eigval[j]*v[i];
		
		double err = arma::norm(w,2);
		if (err>err_max) err_max = err;
		ofile_check << "Norm of Residual:   " << err << endl;
	}
	
	ofile_check << "While tol*sqrt(N):  " << tol*sqrt(Ntot) << endl;
	ofile_check << endl;
	
	return (err_max<tol);
}


void save_Arnoldi_matrices(arma::cx_mat V, arma::cx_mat H)
{
	if (V.n_cols!=H.n_rows || !H.is_square()) throw runtime_error("Inconsistent dimensions in Arnoldi factorisation");
	int k = V.n_cols;
	long Ntot = V.n_rows;
	
	/*
	ofstream ofile_Vmatrix("arnoldi/Vmatrix.dat");
	
	ofile_Vmatrix << fixed << setprecision(6);
	for (long i=0; i<Ntot; i++) 
	{
		for (int j=0; j<k; j++) ofile_Vmatrix << setw(12) << V(i,j).real();
		ofile_Vmatrix << endl;
	}
	
	ofile_Vmatrix << endl;
	for (long i=0; i<Ntot; i++) 
	{
		for (int j=0; j<k; j++) ofile_Vmatrix << setw(12) << V(i,j).imag();
		ofile_Vmatrix << endl;
	}
	*/
	
	ofstream ofile_Hmatrix("arnoldi/Hmatrix.dat");
	
	ofile_Hmatrix << fixed << setprecision(6);
	for (int i=0; i<k; i++) 
	{
		for (int j=0; j<k; j++) ofile_Hmatrix << setw(12) << H(i,j).real();
		ofile_Hmatrix << endl;
	}
	
	ofile_Hmatrix << endl;
	for (int i=0; i<k; i++) 
	{
		for (int j=0; j<k; j++) ofile_Hmatrix << setw(12) << H(i,j).imag();
		ofile_Hmatrix << endl;
	}
}


void Arnoldi::extend_arnoldi_factorisation(arma::cx_mat &V, arma::cx_mat &H, arma::cx_vec &f, const int k, const int p, double shift, double tol) const
{
	const int species = 0;
	const long Ntot = matrix_.get_Ntot();
	
	// Prepare extended V,H matrices
	
	arma::cx_mat V_new(Ntot, k+p); V_new.zeros();
	arma::cx_mat H_new(k+p,  k+p); H_new.zeros();
	
	for (long i=0; i<Ntot; i++) for (int j=0; j<k; j++) V_new(i,j) = V(i,j);
	for (int i=0; i<k; i++) for (int j=0; j<k; j++) H_new(i,j) = H(i,j);
	
	V = V_new;
	H = H_new;
	
	// Initial arnoldi vector
	
	arma::cx_vec v(Ntot);
	
	if (k==0) // random vector
	{
		random_device r;
		int seed = 1; //r();
		
		mt19937 rng(seed);
		uniform_real_distribution<double> urd;
		
		for (long i=0; i<Ntot; i++) v[i] = 2*urd(rng)-1;
		
		if (matrix_.is_fixed_boundary())
		  for (long i=0; i<matrix_.get_Nboundary(); i++)
		    v[matrix_.boundary_pos_2_pos(i)] = 0.0;
		
		v /= norm(v);
		
		ofstream file_guess("arnoldi/guess.dat");
		file_guess << fixed << setprecision(6);
		for (long i=0; i<Ntot; i++) file_guess << setw(10) << real(v[i]); file_guess << endl;
		file_guess << endl;
		for (long i=0; i<Ntot; i++) file_guess << setw(10) << imag(v[i]); file_guess << endl;
		
		arma::cx_vec w; matrix_dot_v(v, w, shift);
		
		for (long j=0; j<Ntot; j++) V(j,0) = v[j];
		H(0,0) = cdot(v,w);
		v = w - v*H(0,0);
	}
	else // use the remainder of the previous factorization
	{
		v = f;
	}
	
	// p-steps extension of V,H matrices
	
	for (int i=(k>0?k-1:0); i<k+p-1; i++)
	{
		double beta = norm(v); 
		v /= beta;
		
		for (long j=0; j<Ntot; j++) V(j,i+1) = v[j];
		for (long j=0; j<=i;   j++) H(i+1,j) = 0; H(i+1,i) = beta;
		
		// new vector
		matrix_dot_v(v, v, shift);
		
		// orthogonalise
		arma::cx_vec h = V.t()*v;
		v -= V*h;
		
		// refine orthogonalisation
		for (int j=0; j<3; j++)
		{
			arma::cx_vec s = V.t()*v;
			v -= V*s;
			h += s;
		}
		
		
		for (int j=0; j<=i+1; j++) H(j,i+1) = h[j];
	}
	
	f = v;
	
	save_Arnoldi_matrices(V,H);
	
	if (!check_factorisation(V,H,f,shift,tol)) 
		throw runtime_error("Arnoldi factorisation does not check out");
}

