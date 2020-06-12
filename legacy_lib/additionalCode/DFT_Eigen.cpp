#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include <stdexcept>

using namespace std;

#include "DFT_LinAlg.h"

#include <Eigen/Dense>

using namespace Eigen;

#define convert(x) (*(static_cast<Eigen::VectorXd*>(data_)))

#define DATA convert(data_)
#define v_DATA (*(static_cast<Eigen::VectorXd*>(v.data_)))
#define v1_DATA (*(static_cast<Eigen::VectorXd*>(v1.data_)))
#define v2_DATA (*(static_cast<Eigen::VectorXd*>(v2.data_)))


DFT_Vec::DFT_Vec(unsigned N) : data_(new Eigen::VectorXd(N)){}
DFT_Vec::DFT_Vec(const DFT_Vec& v) { DATA = v_DATA;}
DFT_Vec::DFT_Vec(): data_(new Eigen::VectorXd(1)) {}

DFT_Vec::~DFT_Vec() { if(data_) delete static_cast<Eigen::VectorXd*>(data_);}

void   DFT_Vec::set(unsigned pos, double val) { DATA[pos] = val;}
double DFT_Vec::get(unsigned pos) const { return DATA[pos];}

void DFT_Vec::set(const DFT_Vec& v) { DATA = v_DATA;} 
void DFT_Vec::set(const DFT_Vec& v1, const DFT_Vec& v2, double scale) { DATA = v1_DATA+v2_DATA*scale;}

void DFT_Vec::setFromAlias(const DFT_Vec &v) { DATA = v_DATA.cwiseProduct(v_DATA); DATA = (DATA.array() + 1e-20).matrix();} 
void DFT_Vec::setAliasFromValues(const DFT_Vec &v)
{
  for(long i=0;i<v.size();i++)
    set(i, sqrt(std::max(0.0, v.get(i)-1e-20)));          
}
void DFT_Vec::alias_Jacobian(const DFT_Vec &v) { DATA = DATA.cwiseProduct(2*v_DATA);}

void DFT_Vec::set(const double *x, unsigned n) { DATA.resize(n); memcpy(DATA.data(),x,sizeof(double)*n);} 
  
void DFT_Vec::resize(long N) {DATA.resize(N);}
void DFT_Vec::zeros(long N)  {DATA = Eigen::VectorXd::Zero(N);}
void DFT_Vec::zeros()  {DATA = Eigen::VectorXd::Zero(DATA.size());}

double DFT_Vec::inf_norm() const { return DATA.lpNorm<Infinity>();}
double DFT_Vec::euclidean_norm() const { return DATA.norm();}
void DFT_Vec::normalise() { MultBy(1.0/euclidean_norm());} 
  
double DFT_Vec::min() const { return DATA.minCoeff();} 
double DFT_Vec::max() const { return DATA.maxCoeff();} 

double *DFT_Vec::memptr() { return DATA.data();}

unsigned DFT_Vec::size() const { return DATA.size();}

double DFT_Vec::dotWith(const DFT_Vec &v) const { return DATA.dot(v_DATA);}

double DFT_Vec::accu() const { return DATA.sum();}

void DFT_Vec::save(ofstream &of) const
{
  of << DATA.size() << endl;
  for(long i=0;i<DATA.size();i++)
    of << DATA[i] << " ";
  of << endl;
}
void DFT_Vec::load(ifstream &in)
{
  unsigned N;
  in >> N;
  zeros(N);
  for(long i=0;i<DATA.size();i++)
    in >> DATA[i];
}
    
void DFT_Vec::MultBy(double val)            { DATA *= val;}
void DFT_Vec::IncrementBy(const DFT_Vec& v) { DATA += v_DATA;}
void DFT_Vec::DecrementBy(const DFT_Vec& v) { DATA -= v_DATA;}  
void DFT_Vec::ShiftBy(double shift)         { DATA = DATA.array() -shift;}

void DFT_Vec::IncrementBy(unsigned pos, double val) { DATA[pos] += val;}

void DFT_Vec::IncrementBy_Scaled_Vector(const DFT_Vec& v,double scale) {DATA += v_DATA*scale;}

void DFT_Vec::Schur(const DFT_Vec &v1, const DFT_Vec &v2) { DATA = v1_DATA.cwiseProduct(v2_DATA);}
  
template<class Archive> void DFT_Vec::save(Archive & ar, const unsigned int version) const
{
  unsigned N = DATA.size();
  boost::serialization::binary_object buf_wrap(DATA.data(), N*sizeof(double));
  ar & N;
  ar & buf_wrap;
}
  
template<class Archive>  void DFT_Vec::load(Archive & ar, const unsigned int version)
{
  unsigned N  = 0;
  ar & N;
  DATA.resize(N);    
  boost::serialization::binary_object buf_wrap(DATA.data(), N*sizeof(double));
  ar & buf_wrap;
}
  
////////////////////////////
// Complex vector

#define cDATA (*(static_cast<Eigen::VectorXcd*>(data_)))
#define v_cDATA (*(static_cast<Eigen::VectorXcd*>(v.data_)))
#define v1_cDATA (*(static_cast<Eigen::VectorXcd*>(v1.data_)))
#define v2_cDATA (*(static_cast<Eigen::VectorXcd*>(v2.data_)))

DFT_Vec_Complex::DFT_Vec_Complex(unsigned N) : data_(new Eigen::VectorXcd(N)){}
DFT_Vec_Complex::DFT_Vec_Complex(const DFT_Vec_Complex& v) { cDATA = v_cDATA;}
DFT_Vec_Complex::DFT_Vec_Complex() : data_(new Eigen::VectorXcd(1)) {}

DFT_Vec_Complex::~DFT_Vec_Complex() { if(data_) delete static_cast<Eigen::VectorXcd*>(data_);}


void DFT_Vec_Complex::set(const DFT_Vec_Complex& v) { cDATA = v_cDATA;}
void DFT_Vec_Complex::set(unsigned pos, complex<double> val) { cDATA[pos] = val;}
complex<double> DFT_Vec_Complex::get(unsigned pos) const { return cDATA[pos];}

void DFT_Vec_Complex::MultBy(double val)  { cDATA *= val;}
  
void DFT_Vec_Complex::Schur(const DFT_Vec_Complex &v1, const DFT_Vec_Complex &v2, bool bUseConj)
{
if(bUseConj) cDATA = v1_cDATA.cwiseProduct(v2_cDATA.conjugate());
else cDATA         = v1_cDATA.cwiseProduct(v2_cDATA);
}
void DFT_Vec_Complex::incrementSchur(const DFT_Vec_Complex &v1, const DFT_Vec_Complex &v2, bool bUseConj)
{
 if(bUseConj) cDATA += v1_cDATA.cwiseProduct(v2_cDATA.conjugate());
 else         cDATA += v1_cDATA.cwiseProduct(v2_cDATA);
}

void DFT_Vec_Complex::resize(long N) {cDATA.resize(N);}
void DFT_Vec_Complex::zeros(long N)  {cDATA = Eigen::VectorXcd::Zero(N);}
void DFT_Vec_Complex::zeros()        {cDATA = Eigen::VectorXcd::Zero(cDATA.size());}

//complex<double> DFT_Vec_Complex::max() const { return cDATA.maxCoeff();}
  
complex<double> *DFT_Vec_Complex::memptr() { return cDATA.data();}

unsigned DFT_Vec_Complex::size() const { return cDATA.size();}

template<class Archive> void DFT_Vec_Complex::save(Archive & ar, const unsigned int version) const
{
  unsigned N = size();
  boost::serialization::binary_object buf_wrap(DATA.data(), N*sizeof(complex<double>));
  ar & N;
  ar & buf_wrap;
}

template<class Archive>  void DFT_Vec_Complex::load(Archive & ar, const unsigned int version)
{
  unsigned N  = 0;
  ar & N;
  resize(N);    
  boost::serialization::binary_object buf_wrap(DATA.data(), N*sizeof(complex<double>));
  ar & buf_wrap;
}
