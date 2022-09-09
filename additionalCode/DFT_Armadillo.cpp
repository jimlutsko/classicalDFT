#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>

using namespace std;

#include "DFT_LinAlg.h"

#include <armadillo>

#define DATA (*(static_cast<arma::vec*>(data_)))
#define v_DATA (*(static_cast<arma::vec*>(v.data_)))
#define v1_DATA (*(static_cast<arma::vec*>(v1.data_)))
#define v2_DATA (*(static_cast<arma::vec*>(v2.data_)))

DFT_Vec::DFT_Vec(unsigned N) : data_(new arma::vec(N)){}
DFT_Vec::DFT_Vec(const DFT_Vec& v) : data_(new arma::vec(v.size())) {DATA = v_DATA;}
DFT_Vec::DFT_Vec(): data_(new arma::vec(1)) {}

DFT_Vec& DFT_Vec::operator= (const DFT_Vec& v){DATA = v_DATA; return *this;}

DFT_Vec::~DFT_Vec() { if(data_) delete static_cast<arma::vec*>(data_);}

void   DFT_Vec::set(unsigned pos, double val) { DATA[pos] = val;}
double DFT_Vec::get(unsigned pos) const { return DATA[pos];}

void DFT_Vec::set(const DFT_Vec& v) { DATA = v_DATA;} 
void DFT_Vec::set(const DFT_Vec& v1, const DFT_Vec& v2, double scale) { DATA = v1_DATA+v2_DATA*scale;}
void DFT_Vec::set(double d) { DATA.fill(d);}  
void DFT_Vec::set(const double *x, unsigned n) { DATA.set_size(n); memcpy(DATA.memptr(),x,sizeof(double)*n);}
void DFT_Vec::set_random_normal(){arma::arma_rng::set_seed_random(); DATA.randn();}
void DFT_Vec::set_random() {arma::arma_rng::set_seed_random(); DATA.randu(); DATA = 2*DATA-1;}

void DFT_Vec::resize(long N) {DATA.resize(N);}
void DFT_Vec::zeros(long N)  {DATA.zeros(N);}
void DFT_Vec::zeros()  {DATA.zeros();}

double DFT_Vec::inf_norm() const { return arma::norm(DATA,"inf");}
double DFT_Vec::euclidean_norm() const { return arma::norm(DATA);}
void DFT_Vec::normalise() { MultBy(1.0/euclidean_norm());} // better data_ = arma::normalise(data_) ???
  
double DFT_Vec::min() const { return DATA.min();} 
double DFT_Vec::max() const { return DATA.max();} 

double *DFT_Vec::memptr() { return DATA.memptr();}

unsigned DFT_Vec::size() const { return DATA.size();}

double DFT_Vec::dotWith(const DFT_Vec &v) const { return arma::dot(v_DATA,DATA);}

double DFT_Vec::accu() const { return arma::accu(DATA);}

void DFT_Vec::MultBy(double val)            { DATA *= val;}
void DFT_Vec::IncrementBy(const DFT_Vec& v) { DATA += v_DATA;}
void DFT_Vec::DecrementBy(const DFT_Vec& v) { DATA -= v_DATA;}  
void DFT_Vec::add(double shift)         { DATA += shift;}

void DFT_Vec::IncrementBy(unsigned pos, double val) { DATA[pos] += val;}

void DFT_Vec::IncrementBy_Scaled_Vector(const DFT_Vec& v,double scale) {DATA += v_DATA*scale;}

void DFT_Vec::Schur(const DFT_Vec &v1, const DFT_Vec &v2) { DATA = v1_DATA%v2_DATA;}

template<class Archive> void DFT_Vec::save(Archive & ar, const unsigned int version) const
{
  unsigned N = DATA.size();
  boost::serialization::binary_object buf_wrap(DATA.memptr(), N*sizeof(double));
  ar & N;
  ar & buf_wrap;
}
  
template<class Archive>  void DFT_Vec::load(Archive & ar, const unsigned int version)
{
  unsigned N  = 0;
  ar & N;
  DATA.resize(N);    
  boost::serialization::binary_object buf_wrap(DATA.memptr(), N*sizeof(double));
  ar & buf_wrap;
}
  

// These are legacy functions that should be removed at some point. 
void DFT_Vec::save(ofstream &of) const {DATA.save(of);}
void DFT_Vec::load(ifstream &in) {DATA.load(in);}


////////////////////////////
// Complex vector

#define cDATA (*(static_cast<arma::cx_vec*>(data_)))
#define v_cDATA (*(static_cast<arma::cx_vec*>(v.data_)))
#define v1_cDATA (*(static_cast<arma::cx_vec*>(v1.data_)))
#define v2_cDATA (*(static_cast<arma::cx_vec*>(v2.data_)))

DFT_Vec_Complex::DFT_Vec_Complex(unsigned N) : data_(new arma::cx_vec(N)){}
DFT_Vec_Complex::DFT_Vec_Complex(const DFT_Vec_Complex& v) : data_(new arma::cx_vec(v.size())) {cDATA = v_cDATA;}
DFT_Vec_Complex::DFT_Vec_Complex() : data_(new arma::cx_vec(1)) {}

DFT_Vec_Complex::~DFT_Vec_Complex() { if(data_) delete static_cast<arma::cx_vec*>(data_);}
  
void DFT_Vec_Complex::set(const DFT_Vec_Complex& v) { cDATA = v_cDATA;}
void DFT_Vec_Complex::set(unsigned pos, complex<double> val) { cDATA[pos] = val;}
complex<double> DFT_Vec_Complex::get(unsigned pos) const { return cDATA[pos];}

void DFT_Vec_Complex::MultBy(double val)  { cDATA *= val;}

void DFT_Vec_Complex::Schur(const DFT_Vec_Complex &v1, const DFT_Vec_Complex &v2, bool bUseConj)
{
 
 if(bUseConj) cDATA = v1_cDATA%conj(v2_cDATA);
 else cDATA = v1_cDATA%v2_cDATA;
}

void DFT_Vec_Complex::incrementSchur(const DFT_Vec_Complex &v1, const DFT_Vec_Complex &v2, bool bUseConj)
{
 if(bUseConj) cDATA += v1_cDATA%conj(v2_cDATA);
 else cDATA += v1_cDATA%v2_cDATA;
}


void DFT_Vec_Complex::resize(long N) {cDATA.resize(N);}
void DFT_Vec_Complex::zeros(long N)  {cDATA.zeros(N);}
void DFT_Vec_Complex::zeros()        {cDATA.zeros();}

complex<double> DFT_Vec_Complex::max() const { return cDATA.max();}
complex<double> DFT_Vec_Complex::min() const { return cDATA.min();}

complex<double> *DFT_Vec_Complex::memptr() { return cDATA.memptr();}
 
unsigned DFT_Vec_Complex::size() const { return cDATA.size();}


template<class Archive> void DFT_Vec_Complex::save(Archive & ar, const unsigned int version) const
{
  unsigned N = size();
  boost::serialization::binary_object buf_wrap(cDATA.memptr(), N*sizeof(complex<double>));
  ar & N;
  ar & buf_wrap;
}

template<class Archive>  void DFT_Vec_Complex::load(Archive & ar, const unsigned int version)
{
  unsigned N  = 0;
  ar & N;
  resize(N);    
  boost::serialization::binary_object buf_wrap(cDATA.memptr(), N*sizeof(complex<double>));
  ar & buf_wrap;
}
