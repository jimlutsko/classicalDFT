#ifndef __DFT_LINALG_ATLAS_LUTSKO__
#define __DFT_LINALG_ATLAS_LUTSKO__

//#include <functional>
#include <numeric>
#include <complex>

using namespace std;

extern "C"
{
       #include <cblas.h>
}

#include "FMT_FFTW.h"

class DFT_Vec
{
 public:
  //Constructor: reserves size of array but does not initialize elements
 DFT_Vec(unsigned N) : data_(N){}

  // Copy constructor
  DFT_Vec(const DFT_Vec& c) { data_.resize(c.size()); cblas_dcopy(c.size(),&c.data_[0],1,&data_[0],1);}

  // Do nothing constructor
  DFT_Vec() {}

  /**
   *  @brief set data_[pos] = val
   *  @param pos unsigned int
   *  @param val double - the value assigned
   *  @return nothing
   */
  void   set(unsigned pos, double val) { data_[pos] = val;}

  /**
   *  @brief returns data_[pos]
   *  @param pos unsigned int
   *  @return data_[pos]
   */
  double get(unsigned pos) const { return data_[pos];}

  /**
   *  @brief copies a vector.
   *  @detailed for all i, set data_[i] = x[i]. Resizes data_ if necessary.
   *  @param pos unsigned int
   *  @param x DFT_Vec - the vector copied
   *  @return nothing
   */
  void set(const DFT_Vec& x) { data_.resize(x.size()); cblas_dcopy(x.size(),&x.data_[0],1,&data_[0],1);}

  /**
   *  @brief for all i, sets data_[i] = v1[i] + scale*v2[i]. Resiwes data_ if necessary.
   *  @detailed barfs if either &v1 = this or &v2 = this
   *  @param v1 DFT_Vec 
   *  @param v2 DFT_Vec 
   *  @param scale double
   *  @return nothing
   */  
  void set(const DFT_Vec& v1, const DFT_Vec& v2, double scale) 
  {
    if(data_.size() != v1.size()) data_.resize(v1.size());

    if(&v1 == this || &v2 == this)throw std::runtime_error("Bad arguments in DFT_Vec: set");

    cblas_dcopy(v1.size(),&v1.data_[0],1,&data_[0],1);
    cblas_daxpy(v2.size(),scale,&v2.data_[0],1,&data_[0],1);
  } 

  void setFromAmplitude(double v, const DFT_Vec &x) // problem!
  { 
    data_.resize(x.size());
    for(int i=0;i<x.size();i++) {double z = x.get(i); data_[i] = v+z*z;}
  }
  
  void addTo(DFT_Vec &v) { cblas_daxpy(v.size(), 1.0 ,&v.data_[0],1,&data_[0],1);}

  void addTo(unsigned pos, double val) { data_[pos] += val;}

  void addTo(double val) // problem!
  { for(int i=0;i<data_.size(); i++) data_[i] += val;}
  
  void resize(long N) {data_.resize(N);}
  void zeros(long N)  {resize(N); std::fill(data_.begin(), data_.end(), 0.0);} 

  double inf_norm() const { unsigned location = cblas_idamax(data_.size(),&data_[0],1); return fabs(data_[location]);}
  //  double min() const 
  //  double max() const { return arma::max(data_);}

  double *memptr() { return &data_[0];}

  unsigned size() const { return data_.size();}
  
  double dotWith(const DFT_Vec &v) const { return cblas_ddot(data_.size(), &data_[0], 1, &v.data_[0],1);} 

  void multBy(double val) { cblas_dscal(data_.size(), val, &data_[0],1);} 
  void scaleBy(long pos, double val) { data_[pos] /= val;} //cblas_dscal(data_.size(), 1.0/val, &data_[0],1);} 

  //  double accu() const { return arma::accu(data_);}
  double accu_abs_vals() const { return cblas_dasum(data_.size(),&data_[0],1);} 

  double accu() const { return accumulate(begin(data_), end(data_), 0.0);} 

  void Schur(const DFT_Vec &v1, const DFT_Vec &v2) // problem!
  { for(int i=0;i<v1.size();i++) data_[i] = v1.data_[i]*v2.data_[i];}

  void save(ofstream &of) const { long N = data_.size(); of.write((char*) &N, sizeof(long)); of.write((char*) &data_[0], N*sizeof(double));} 
  void load(ifstream &in) { long N = 0; in.read((char*) &N, sizeof(long)); data_.resize(N); in.read((char*) &data_[0], N*sizeof(double));} 

  void Increment_Shift_And_Scale(const DFT_Vec& v,double scale, double shift) // problem!!
  {
    cblas_daxpy(v.size(),scale,&v.data_[0],1,&data_[0],1);
    addTo(-shift*scale);
  }

  void Increment_And_Scale(const DFT_Vec& v,double scale)
  {
    cblas_daxpy(v.size(),scale,&v.data_[0],1,&data_[0],1); 
  }

  void IncrementBy(const DFT_Vec& v) { cblas_daxpy(v.size(),1.0,&v.data_[0],1,&data_[0],1);} 

 protected:
  vector<double> data_;
};



class DFT_Vec_Complex
{
 public:
 DFT_Vec_Complex(unsigned N) : data_(N){}
  DFT_Vec_Complex(const DFT_Vec_Complex& c) { data_.resize(c.size()); cblas_zcopy(c.size(),&c.data_[0],1,&data_[0],1);} 
  DFT_Vec_Complex() {}

  void           set(unsigned pos, complex<double> val) { data_[pos] = val;}
  complex<double> get(unsigned pos) const { return data_[pos];}

  void scaleBy(double val) { complex<double> s = 1.0/val; cblas_zscal(data_.size(), (const void*) &s, &data_[0],1);} 
  void multBy(double val)  { complex<double> s = val; cblas_zscal(data_.size(), (const void*) &s, &data_[0],1);} 
  
  void Schur(const DFT_Vec_Complex &v1, const DFT_Vec_Complex &v2, bool bUseConj=false) // problem!!
  {
    for(int i=0;i<v1.size();i++)
      data_[i] = v1.data_[i]*(bUseConj ? conj(v2.data_[i]) : v2.data_[i]);
  }
  void incrementSchur(const DFT_Vec_Complex &v1, const DFT_Vec_Complex &v2) // problem!!
  { 
    for(int i=0;i<v1.size();i++)
      data_[i] += v1.data_[i]*v2.data_[i];
  }
 
  void resize(long N) {data_.resize(N);}
  void zeros(long N)  {resize(N); std::fill(data_.begin(), data_.end(), 0.0);} 
  void zeros()        {std::fill(data_.begin(), data_.end(), 0.0);} 
 
  complex<double> *memptr() { return &data_[0];} 
 
  unsigned size() const { return data_.size();}
  
 protected:
  vector< complex<double> > data_;
};


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


#endif // __DFT_LINALG_ATLAS_LUTSKO__
