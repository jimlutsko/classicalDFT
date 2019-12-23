#ifndef __DFT_LINALG_EIGEN_LUTSKO__
#define __DFT_LINALG_EIGEN_LUTSKO__


#include <Eigen/Dense>

using namespace Eigen;

#include "FMT_FFTW.h"



/**
  *  @brief UTILITY: A wrapper for linear algebra packages. The idea is to be able to easily change libraries without changing any other code. Probably overkill and could be eliminated.
  */  
class DFT_Vec
{
 public:
 DFT_Vec(unsigned N) : data_(N){}
  DFT_Vec(const DFT_Vec& c) { data_ = c.data_;}
  DFT_Vec() {}
  
  void   set(unsigned pos, double val) { data_[pos] = val;}
  double get(unsigned pos) const { return data_[pos];}

  void set(const DFT_Vec& x) { data_ = x.data_;} 
  void set(const DFT_Vec& v1, const DFT_Vec& v2, double scale) { data_ = v1.data_+v2.data_*scale;}
  
  void setFromAlias(const DFT_Vec &x) { data_ = x.data_.cwiseProduct(x.data_); data_ = (data_.array() + 1e-20).matrix();} //{ data_ = 1e-10+exp(x.data_);}
  void setAliasFromValues(const DFT_Vec &x)
  {
    for(long i=0;i<x.size();i++)
      set(i, sqrt(std::max(0.0, x.get(i)-1e-20)));          
    //      set(i, log(std::max(1e-15, x.get(i)-1e-10)));    
  }
  void alias_Jacobian(const DFT_Vec &x) { data_ = data_.cwiseProduct(2*x.data_);}

  
  void set(const double *x, unsigned n) { data_.resize(n); memcpy(data_.data(),x,sizeof(double)*n);} 
  
  void resize(long N) {data_.resize(N);}
  void zeros(long N)  {data_ = Eigen::VectorXd::Zero(N);}
  void zeros()  {data_ = Eigen::VectorXd::Zero(data_.size());}

  double inf_norm() const { return data_.lpNorm<Infinity>();}
  double euclidean_norm() const { return data_.norm();}
  void normalise() { MultBy(1.0/euclidean_norm());} // better data_ = arma::normalise(data_) ???
  
  double min() const { return data_.minCoeff();} 
  double max() const { return data_.maxCoeff();} 

  double *memptr() { return data_.data();}

  unsigned size() const { return data_.size();}
  
  double dotWith(const DFT_Vec &v) const { return data_.dot(v.data_);}

  double accu() const { return data_.sum();}

  void save(ofstream &of) const
  {
    of << data_.size() << endl;
    for(long i=0;i<data_.size();i++)
      of << data_[i] << " ";
    of << endl;
  }
  void load(ifstream &in)
  {
    unsigned N;
    in >> N;
    zeros(N);
    for(long i=0;i<data_.size();i++)
      in >> data_[i];
  }
    

  void MultBy(double val)            { data_ *= val;}
  void IncrementBy(const DFT_Vec& v) { data_ += v.data_;}
  void DecrementBy(const DFT_Vec& v) { data_ -= v.data_;}  
  void ShiftBy(double shift)         { data_ = data_.array() -shift;}

  void IncrementBy(unsigned pos, double val) { data_[pos] += val;}
    
  void IncrementBy_Scaled_Vector(const DFT_Vec& v,double scale) {data_ += v.data_*scale;}

  void Schur(const DFT_Vec &v1, const DFT_Vec &v2) { data_ = v1.data_.cwiseProduct(v2.data_);}
  
 protected:
  Eigen::VectorXd data_;
};


/**
  *  @brief UTILITY: A wrapper for linear algebra packages. This one handles complex values.
  */  
class DFT_Vec_Complex
{
 public:
 DFT_Vec_Complex(unsigned N) : data_(N){}
  DFT_Vec_Complex(const DFT_Vec_Complex& c) { data_ = c.data_;}
  DFT_Vec_Complex() {}

  void   set(const DFT_Vec_Complex& c) { data_ = c.data_;}
  void   set(unsigned pos, complex<double> val) { data_[pos] = val;}
  complex<double> get(unsigned pos) const { return data_[pos];}

  void MultBy(double val)  { data_ *= val;}
  
  void Schur(const DFT_Vec_Complex &v1, const DFT_Vec_Complex &v2, bool bUseConj=false) { if(bUseConj) data_ = v1.data_.cwiseProduct(v2.data_.conjugate()); else data_ = v1.data_.cwiseProduct(v2.data_);}
  void incrementSchur(const DFT_Vec_Complex &v1, const DFT_Vec_Complex &v2) { data_ += v1.data_.cwiseProduct(v2.data_);}
 
  void resize(long N) {data_.resize(N);}
  void zeros(long N)  {data_ = Eigen::VectorXcd::Zero(N);}
  void zeros()        {data_ = Eigen::VectorXcd::Zero(data_.size());}

  //complex<double> max() const { return data_.maxCoeff();}
  
  complex<double> *memptr() { return data_.data();}
 
  unsigned size() const { return data_.size();}
 
 protected:
  Eigen::VectorXcd data_;
};

/**
  *  @brief UTILITY: This class encapsulates a basic FFT. It is basically an interface to fftw library while holding both a real and a complex vector.
  */  
class DFT_FFT
{
 public:
 DFT_FFT(unsigned Nx, unsigned Ny, unsigned Nz) : RealSpace_(Nx*Ny*Nz), FourierSpace_(Nx*Ny*((Nz/2)+1)), four_2_real_(NULL), real_2_four_(NULL)
    {
      four_2_real_ = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, reinterpret_cast<fftw_complex*>(FourierSpace_.memptr()), RealSpace_.memptr(), FMT_FFTW);
      real_2_four_ = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, RealSpace_.memptr(),reinterpret_cast<fftw_complex*>(FourierSpace_.memptr()),  FMT_FFTW);
    };

 DFT_FFT() : four_2_real_(NULL), real_2_four_(NULL){};
  
  ~DFT_FFT()
    {
      if(four_2_real_) fftw_destroy_plan(four_2_real_);
      if(real_2_four_) fftw_destroy_plan(real_2_four_);
    }

  void zeros() {RealSpace_.zeros(); FourierSpace_.zeros();}
  
  void initialize(unsigned Nx, unsigned Ny, unsigned Nz)
  {
    RealSpace_.zeros(Nx*Ny*Nz);
    FourierSpace_.zeros(Nx*Ny*((Nz/2)+1));
    four_2_real_ = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, reinterpret_cast<fftw_complex*>(FourierSpace_.memptr()), RealSpace_.memptr(), FMT_FFTW);
    real_2_four_ = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, RealSpace_.memptr(),reinterpret_cast<fftw_complex*>(FourierSpace_.memptr()),  FMT_FFTW);
  };
  
  DFT_Vec &Real() { return RealSpace_;}
  DFT_Vec_Complex &Four() { return FourierSpace_;}

  const DFT_Vec &cReal() const { return RealSpace_;}
  const DFT_Vec_Complex &cFour() const { return FourierSpace_;}

  void do_real_2_fourier() {fftw_execute(real_2_four_);}
  void do_fourier_2_real() {fftw_execute(four_2_real_);}

 protected:
  DFT_Vec RealSpace_;
  DFT_Vec_Complex FourierSpace_;
  fftw_plan four_2_real_;   ///< fftw3 "plan" to perform fft on density
  fftw_plan real_2_four_;   ///< fftw3 "plan" to perform fft on density
};


#endif // __DFT_LINALG_EIGEN_LUTSKO__
