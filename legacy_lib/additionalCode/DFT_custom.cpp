#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>
#include <vector>

using namespace std;

#include "DFT_LinAlg.h"

#include <armadillo>

#define VEC (*(static_cast<vector<double>*>(data_)))
#define v_VEC (*(static_cast<vector<double>*>(v.data_)))
#define v1_VEC (*(static_cast<vector<double>*>(v1.data_)))
#define v2_VEC (*(static_cast<vector<double>*>(v2.data_)))

#define DATA (*(static_cast<arma::vec*>(data_)))
#define v_DATA (*(static_cast<arma::vec*>(v.data_)))
#define v1_DATA (*(static_cast<arma::vec*>(v1.data_)))
#define v2_DATA (*(static_cast<arma::vec*>(v2.data_)))

DFT_Vec::DFT_Vec(unsigned N) : data_(new vector<double>(N,0.0)){} // TODO: remains slow
DFT_Vec::DFT_Vec(const DFT_Vec& v) {VEC = v_VEC;}
DFT_Vec::DFT_Vec(): data_(new vector<double>(1,0.0)) {}

DFT_Vec::~DFT_Vec() { if(data_) delete static_cast<vector<double>*>(data_);}

void   DFT_Vec::set(unsigned pos, double val) { VEC[pos] = val;}
double DFT_Vec::get(unsigned pos) const { return VEC[pos];}

void DFT_Vec::set(const DFT_Vec& v) {VEC = v_VEC;} // TODO: remains slow
void DFT_Vec::set(const DFT_Vec& v1, const DFT_Vec& v2, double scale)
{
	long N = VEC.size();
	//for (long i=0;i<N;i++) VEC[i] = v1_VEC[i] + v2_VEC[i] * scale;
	
	long Nq = CACHESIZE/4/sizeof(double);
	int Q = N/Nq;
	
	for (int q=0; q<Q; q++)
	{
		#pragma omp parallel for
		for (long i=Nq*q; i<Nq*(q+1); i++) VEC[i] = v1_VEC[i] + v2_VEC[i] * scale;
	}
	
	if (N%Nq!=0)
	{
		#pragma omp parallel for
		for (long i=Nq*Q; i<N; i++) VEC[i] = v1_VEC[i] + v2_VEC[i] * scale;
	}
}
  
void DFT_Vec::set(const double *x, unsigned N)
{
	VEC.resize(N);
	//for (long i=0; i<N; i++) VEC[i] = x[i];
	
	long Nq = CACHESIZE/4/sizeof(double);
	int Q = N/Nq;
	
	for (int q=0; q<Q; q++)
	{
		#pragma omp parallel for
		for (long i=Nq*q; i<Nq*(q+1); i++) VEC[i] = x[i];
	}
	
	if (N%Nq!=0)
	{
		#pragma omp parallel for
		for (long i=Nq*Q; i<N; i++) VEC[i] = x[i];
	}
}
  
void DFT_Vec::resize(long N) {VEC.resize(N);}
void DFT_Vec::zeros(long N)
{
	VEC.resize(N);
	//for (long i=0; i<N; i++) VEC[i] = 0.0;
	
	long Nq = CACHESIZE/4/sizeof(double);
	int Q = N/Nq;
	
	for (int q=0; q<Q; q++)
	{
		#pragma omp parallel for
		for (long i=Nq*q; i<Nq*(q+1); i++) VEC[i] = 0.0;
	}
	
	if (N%Nq!=0)
	{
		#pragma omp parallel for
		for (long i=Nq*Q; i<N; i++) VEC[i] = 0.0;
	}
}
void DFT_Vec::zeros()
{
	long N = VEC.size();
	//for (long i=0; i<N; i++) VEC[i] = 0.0;
	
	long Nq = CACHESIZE/4/sizeof(double);
	int Q = N/Nq;
	
	for (int q=0; q<Q; q++)
	{
		#pragma omp parallel for
		for (long i=Nq*q; i<Nq*(q+1); i++) VEC[i] = 0.0;
	}
	
	if (N%Nq!=0)
	{
		#pragma omp parallel for
		for (long i=Nq*Q; i<N; i++) VEC[i] = 0.0;
	}
}

double DFT_Vec::inf_norm() const {throw runtime_error("Inf norm not implemented in DFT_custom.cpp");}
double DFT_Vec::euclidean_norm() const 
{
	// TODO: proper summation
	
	double sum = 0.0;
	
	long N = VEC.size();
	long Nq = CACHESIZE/4/sizeof(double);
	int Q = N/Nq;
	
	//for (long i=0; i<N; i++) sum += VEC[i]*VEC[i];
	
	for (int q=0; q<Q; q++)
	{
		#pragma omp parallel for
		for (long i=Nq*q; i<Nq*(q+1); i++) sum += VEC[i];
	}
	
	if (N%Nq!=0)
	{
		#pragma omp parallel for
		for (long i=Nq*Q; i<N; i++) sum += VEC[i];
	}
	
	return sum;
}

void DFT_Vec::normalise() { MultBy(1.0/euclidean_norm());}
  
double DFT_Vec::min() const
{
	long N = VEC.size();
	long Nq = CACHESIZE/4/sizeof(double);
	int Q = N/Nq;
	
	double min_value = VEC[0];
	
	//for (long i=1; i<N; i++) if (VEC[i]<min_value) min_value = VEC[i];
	
	for (int q=0; q<Q; q++)
	{
		#pragma omp parallel for
		for (long i=Nq*q; i<Nq*(q+1); i++) if (VEC[i]<min_value) min_value = VEC[i];
	}
	
	if (N%Nq!=0)
	{
		#pragma omp parallel for
		for (long i=Nq*Q; i<N; i++) if (VEC[i]<min_value) min_value = VEC[i];
	}
	
	return min_value;
}

double DFT_Vec::max() const
{
	long N = VEC.size();
	long Nq = CACHESIZE/4/sizeof(double);
	int Q = N/Nq;
	
	double max_value = VEC[0];
	
	//for (long i=1; i<N; i++) if (VEC[i]>max_value) max_value = VEC[i];
	
	for (int q=0; q<Q; q++)
	{
		#pragma omp parallel for
		for (long i=Nq*q; i<Nq*(q+1); i++) if (VEC[i]>max_value) max_value = VEC[i];
	}
	
	if (N%Nq!=0)
	{
		#pragma omp parallel for
		for (long i=Nq*Q; i<N; i++) if (VEC[i]>max_value) max_value = VEC[i];
	}
	
	return max_value;
}

double *DFT_Vec::memptr() {return VEC.data();} // TODO: is this correct??

unsigned DFT_Vec::size() const { return VEC.size();}

double DFT_Vec::dotWith(const DFT_Vec &v) const
{
	// TODO: proper summation
	
	double sum = 0.0;
	
	long N = VEC.size();
	long Nq = CACHESIZE/4/sizeof(double);
	int Q = N/Nq;
	
	//for (long i=0; i<N; i++) sum += VEC[i]*v_VEC[i];
	
	for (int q=0; q<Q; q++)
	{
		#pragma omp parallel for
		for (long i=Nq*q; i<Nq*(q+1); i++) sum += VEC[i] * v_VEC[i];
	}
	
	if (N%Nq!=0)
	{
		#pragma omp parallel for
		for (long i=Nq*Q; i<N; i++) sum += VEC[i] * v_VEC[i];
	}
	
	return sum;
}

double DFT_Vec::accu() const
{
	// TODO: is this a usual sum?
	// TODO: do a proper sum with sorting instead?
	
	double sum = 0.0;
	
	long N = VEC.size();
	long Nq = CACHESIZE/4/sizeof(double);
	int Q = N/Nq;
	
	//for (long i=0; i<N; i++) sum += VEC[i];
	
	for (int q=0; q<Q; q++)
	{
		#pragma omp parallel for
		for (long i=Nq*q; i<Nq*(q+1); i++) sum += VEC[i];
	}
	
	if (N%Nq!=0)
	{
		#pragma omp parallel for
		for (long i=Nq*Q; i<N; i++) sum += VEC[i];
	}
	
	return sum;
}

void DFT_Vec::MultBy(double val)
{
	//long N = VEC.size();
	//for (long i=0; i<N; i++) VEC[i] *= val;
	
	long N = VEC.size();
	long Nq = CACHESIZE/4/sizeof(double);
	int Q = N/Nq;
	
	for (int q=0; q<Q; q++)
	{
		#pragma omp parallel for
		for (long i=Nq*q; i<Nq*(q+1); i++) VEC[i] *= val;
	}
	
	if (N%Nq!=0)
	{
		#pragma omp parallel for
		for (long i=Nq*Q; i<N; i++) VEC[i] *= val;
	}
}
void DFT_Vec::IncrementBy(const DFT_Vec& v)
{
	//long N = VEC.size();
	//for (long i=0; i<N; i++) VEC[i] += v_VEC[i];
	
	long N = VEC.size();
	long Nq = CACHESIZE/4/sizeof(double);
	int Q = N/Nq;
	
	for (int q=0; q<Q; q++)
	{
		#pragma omp parallel for
		for (long i=Nq*q; i<Nq*(q+1); i++) VEC[i] += v_VEC[i];
	}
	
	if (N%Nq!=0)
	{
		#pragma omp parallel for
		for (long i=Nq*Q; i<N; i++) VEC[i] += v_VEC[i];
	}
}
void DFT_Vec::DecrementBy(const DFT_Vec& v)
{
	//long N = VEC.size();
	//for (long i=0; i<N; i++) VEC[i] -= v_VEC[i];
	
	long N = VEC.size();
	long Nq = CACHESIZE/4/sizeof(double);
	int Q = N/Nq;
	
	for (int q=0; q<Q; q++)
	{
		#pragma omp parallel for
		for (long i=Nq*q; i<Nq*(q+1); i++) VEC[i] -= v_VEC[i];
	}
	
	if (N%Nq!=0)
	{
		#pragma omp parallel for
		for (long i=Nq*Q; i<N; i++) VEC[i] -= v_VEC[i];
	}
}  
void DFT_Vec::ShiftBy(double shift)
{
	//long N = VEC.size();
	//for (long i=0; i<N; i++) VEC[i] += -shift; // TODO: Why +-? Is it really more efficient?
	
	long N = VEC.size();
	long Nq = CACHESIZE/4/sizeof(double);
	int Q = N/Nq;
	
	for (int q=0; q<Q; q++)
	{
		#pragma omp parallel for
		for (long i=Nq*q; i<Nq*(q+1); i++) VEC[i] += -shift;
	}
	
	if (N%Nq!=0)
	{
		#pragma omp parallel for
		for (long i=Nq*Q; i<N; i++) VEC[i] += -shift;
	}
}

void DFT_Vec::IncrementBy(unsigned pos, double val) {VEC[pos] += val;}

void DFT_Vec::IncrementBy_Scaled_Vector(const DFT_Vec& v,double scale)
{
	//long N = VEC.size();
	//for (long i=0; i<N; i++) VEC[i] += v_VEC[i]*scale;
	
	long N = VEC.size();
	long Nq = CACHESIZE/4/sizeof(double);
	int Q = N/Nq;
	
	for (int q=0; q<Q; q++)
	{
		#pragma omp parallel for
		for (long i=Nq*q; i<Nq*(q+1); i++) VEC[i] += v_VEC[i]*scale;
	}
	
	if (N%Nq!=0)
	{
		#pragma omp parallel for
		for (long i=Nq*Q; i<N; i++) VEC[i] += v_VEC[i]*scale;
	}
}

void DFT_Vec::Schur(const DFT_Vec &v1, const DFT_Vec &v2)
{
	//int N = VEC.size();
	//for (int i=0; i<N; i++) VEC[i] = v1_VEC[i] * v2_VEC[i];
	
	long N = VEC.size();
	long Nq = CACHESIZE/4/sizeof(double);
	int Q = N/Nq;
	
	for (int q=0; q<Q; q++)
	{
		#pragma omp parallel for
		for (long i=Nq*q; i<Nq*(q+1); i++) VEC[i] = v1_VEC[i] * v2_VEC[i];
	}
	
	if (N%Nq!=0)
	{
		#pragma omp parallel for
		for (long i=Nq*Q; i<N; i++) VEC[i] = v1_VEC[i] * v2_VEC[i];
	}
}

template<class Archive> void DFT_Vec::save(Archive & ar, const unsigned int version) const
{
  unsigned N = VEC.size();
  boost::serialization::binary_object buf_wrap(VEC.data(), N*sizeof(double));
  ar & N;
  ar & buf_wrap;
}
  
template<class Archive>  void DFT_Vec::load(Archive & ar, const unsigned int version)
{
  unsigned N = 0;
  ar & N;
  VEC.resize(N);    
  boost::serialization::binary_object buf_wrap(VEC.data(), N*sizeof(double));
  ar & buf_wrap;
}
  

// These are legacy functions that should be removed at some point. 
void DFT_Vec::save(ofstream &of) const {DATA.save(of);}
void DFT_Vec::load(ifstream &in) {DATA.load(in);}





////////////////////////////
// Complex vector

#define cVEC (*(static_cast<vector<complex<double>>*>(data_)))
#define v_cVEC (*(static_cast<vector<complex<double>>*>(v.data_)))
#define v1_cVEC (*(static_cast<vector<complex<double>>*>(v1.data_)))
#define v2_cVEC (*(static_cast<vector<complex<double>>*>(v2.data_)))

#define cDATA (*(static_cast<arma::cx_vec*>(data_)))
#define v_cDATA (*(static_cast<arma::cx_vec*>(v.data_)))
#define v1_cDATA (*(static_cast<arma::cx_vec*>(v1.data_)))
#define v2_cDATA (*(static_cast<arma::cx_vec*>(v2.data_)))

DFT_Vec_Complex::DFT_Vec_Complex(unsigned N) : data_(new vector<complex<double>>(N,0.0)){}
DFT_Vec_Complex::DFT_Vec_Complex(const DFT_Vec_Complex& v) {cVEC = v_cVEC;}
DFT_Vec_Complex::DFT_Vec_Complex() : data_(new vector<complex<double>>(1,0.0)) {}

DFT_Vec_Complex::~DFT_Vec_Complex() { if(data_) delete static_cast<vector<complex<double>>*>(data_);}
  
void DFT_Vec_Complex::set(const DFT_Vec_Complex& v) {cVEC = v_cVEC;}
void DFT_Vec_Complex::set(unsigned pos, complex<double> val) { cVEC[pos] = val;}
complex<double> DFT_Vec_Complex::get(unsigned pos) const { return cVEC[pos];}

void DFT_Vec_Complex::MultBy(double val)
{
	long N = cVEC.size();
	//for (long i=0; i<N; i++) cVEC[i] *= val;
	
	long Nq = CACHESIZE/4/sizeof(double);
	int Q = N/Nq;
	
	for (int q=0; q<Q; q++)
	{
		#pragma omp parallel for
		for (long i=Nq*q; i<Nq*(q+1); i++) cVEC[i] *= val;
	}
	
	if (N%Nq!=0)
	{
		#pragma omp parallel for
		for (long i=Nq*Q; i<N; i++) cVEC[i] *= val;
	}
}

void DFT_Vec_Complex::Schur(const DFT_Vec_Complex &v1, const DFT_Vec_Complex &v2, bool bUseConj)
{
	if (bUseConj)
	{
		//long N = cVEC.size();
		//for (long i=0; i<N; i++) cVEC[i] = v1_cVEC[i] * conj(v2_cVEC[i]);
		
		long N = cVEC.size();
		long Nq = CACHESIZE/4/sizeof(double);
		int Q = N/Nq;
		
		for (int q=0; q<Q; q++)
		{
			#pragma omp parallel for
			for (long i=Nq*q; i<Nq*(q+1); i++) cVEC[i] = v1_cVEC[i] * conj(v2_cVEC[i]);
		}
		
		if (N%Nq!=0)
		{
			#pragma omp parallel for
			for (long i=Nq*Q; i<N; i++) cVEC[i] = v1_cVEC[i] * conj(v2_cVEC[i]);
		}
	}
	else 
	{
		//long N = cVEC.size();
		//for (long i=0; i<N; i++) cVEC[i] = v1_cVEC[i] * v2_cVEC[i];
		
		long N = cVEC.size();
		long Nq = CACHESIZE/4/sizeof(double);
		int Q = N/Nq;
		
		for (int q=0; q<Q; q++)
		{
			#pragma omp parallel for
			for (long i=Nq*q; i<Nq*(q+1); i++) cVEC[i] = v1_cVEC[i] * v2_cVEC[i];
		}
		
		if (N%Nq!=0)
		{
			#pragma omp parallel for
			for (long i=Nq*Q; i<N; i++) cVEC[i] = v1_cVEC[i] * v2_cVEC[i];
		}
	}
}

void DFT_Vec_Complex::incrementSchur(const DFT_Vec_Complex &v1, const DFT_Vec_Complex &v2, bool bUseConj)
{
	if (bUseConj)
	{
		//long N = cVEC.size();
		//for (long i=0; i<N; i++) cVEC[i] += v1_cVEC[i] * conj(v2_cVEC[i]);
		
		long N = cVEC.size();
		long Nq = CACHESIZE/4/sizeof(double);
		int Q = N/Nq;
		
		for (int q=0; q<Q; q++)
		{
			#pragma omp parallel for
			for (long i=Nq*q; i<Nq*(q+1); i++) cVEC[i] += v1_cVEC[i] * conj(v2_cVEC[i]);
		}
		
		if (N%Nq!=0)
		{
			#pragma omp parallel for
			for (long i=Nq*Q; i<N; i++) cVEC[i] += v1_cVEC[i] * conj(v2_cVEC[i]);
		}
	}
	else 
	{
		//long N = cVEC.size();
		//for (long i=0; i<N; i++) cVEC[i] += v1_cVEC[i] * v2_cVEC[i];
		
		long N = cVEC.size();
		long Nq = CACHESIZE/4/sizeof(double);
		int Q = N/Nq;
		
		for (int q=0; q<Q; q++)
		{
			#pragma omp parallel for
			for (long i=Nq*q; i<Nq*(q+1); i++) cVEC[i] += v1_cVEC[i] * v2_cVEC[i];
		}
		
		if (N%Nq!=0)
		{
			#pragma omp parallel for
			for (long i=Nq*Q; i<N; i++) cVEC[i] += v1_cVEC[i] * v2_cVEC[i];
		}
	}
}

void DFT_Vec_Complex::resize(long N) {cVEC.resize(N);}
void DFT_Vec_Complex::zeros(long N)
{
	cVEC.resize(N);
	//for (long i=0; i<N; i++) cVEC[i] = 0.0;
	
	long Nq = CACHESIZE/4/sizeof(double);
	int Q = N/Nq;
	
	for (int q=0; q<Q; q++)
	{
		#pragma omp parallel for
		for (long i=Nq*q; i<Nq*(q+1); i++) cVEC[i] = 0.0;
	}
	
	if (N%Nq!=0)
	{
		#pragma omp parallel for
		for (long i=Nq*Q; i<N; i++) cVEC[i] = 0.0;
	}
}
void DFT_Vec_Complex::zeros()
{
	long N = cVEC.size();
	//for (long i=0; i<N; i++) cVEC[i] = 0.0;
	
	long Nq = CACHESIZE/4/sizeof(double);
	int Q = N/Nq;
	
	for (int q=0; q<Q; q++)
	{
		#pragma omp parallel for
		for (long i=Nq*q; i<Nq*(q+1); i++) cVEC[i] = 0.0;
	}
	
	if (N%Nq!=0)
	{
		#pragma omp parallel for
		for (long i=Nq*Q; i<N; i++) cVEC[i] = 0.0;
	}
}

complex<double> DFT_Vec_Complex::max() const
{
	long N = cVEC.size();
	complex<double> max_value = cVEC[0];
	
	//for (long i=1; i<N; i++) if (abs(cVEC[i])>abs(max_value)) max_value = cVEC[i];
	
	long Nq = CACHESIZE/4/sizeof(double);
	int Q = N/Nq;
	
	for (int q=0; q<Q; q++)
	{
		#pragma omp parallel for
		for (long i=Nq*q; i<Nq*(q+1); i++) if (abs(cVEC[i])>abs(max_value)) max_value = cVEC[i];
	}
	
	if (N%Nq!=0)
	{
		#pragma omp parallel for
		for (long i=Nq*Q; i<N; i++) if (abs(cVEC[i])>abs(max_value)) max_value = cVEC[i];
	}
	
	return max_value;
}

complex<double> *DFT_Vec_Complex::memptr() { return cVEC.data();}
 
unsigned DFT_Vec_Complex::size() const { return cVEC.size();}


template<class Archive> void DFT_Vec_Complex::save(Archive & ar, const unsigned int version) const
{
  unsigned N = size();
  boost::serialization::binary_object buf_wrap(cVEC.data(), N*sizeof(complex<double>));
  ar & N;
  ar & buf_wrap;
}

template<class Archive>  void DFT_Vec_Complex::load(Archive & ar, const unsigned int version)
{
  unsigned N  = 0;
  ar & N;
  resize(N);    
  boost::serialization::binary_object buf_wrap(cVEC.data(), N*sizeof(complex<double>));
  ar & buf_wrap;
}
