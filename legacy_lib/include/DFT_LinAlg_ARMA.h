#ifndef __DFT_LINALG_ARMA_LUTSKO__
#define __DFT_LINALG_ARMA_LUTSKO__

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/binary_object.hpp>


#include <armadillo>

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
  
  void setFromAlias(const DFT_Vec &x) { data_ = 1e-20+x.data_%x.data_;} 
  void setAliasFromValues(const DFT_Vec &x)
  {
    for(long i=0;i<x.size();i++)
      set(i, sqrt(std::max(0.0, x.get(i)-1e-20)));          
  }
  void alias_Jacobian(const DFT_Vec &x) { data_ %= 2*x.data_;}

  
  void set(const double *x, unsigned n) { data_.set_size(n); memcpy(data_.memptr(),x,sizeof(double)*n);} 
  
  void resize(long N) {data_.resize(N);}
  void zeros(long N)  {data_.zeros(N);}
  void zeros()  {data_.zeros();}

  double inf_norm() const { return arma::norm(data_,"inf");}
  double euclidean_norm() const { return arma::norm(data_);}
  void normalise() { MultBy(1.0/euclidean_norm());} // better data_ = arma::normalise(data_) ???
  
  double min() const { return data_.min();} 
  double max() const { return data_.max();} 

  double *memptr() { return data_.memptr();}

  unsigned size() const { return data_.size();}
  
  double dotWith(const DFT_Vec &v) const { return arma::dot(v.data_,data_);}

  double accu() const { return arma::accu(data_);}

  void MultBy(double val)            { data_ *= val;}
  void IncrementBy(const DFT_Vec& v) { data_ += v.data_;}
  void DecrementBy(const DFT_Vec& v) { data_ -= v.data_;}  
  void ShiftBy(double shift)         { data_ += -shift;}

  void IncrementBy(unsigned pos, double val) { data_[pos] += val;}
    
  void IncrementBy_Scaled_Vector(const DFT_Vec& v,double scale) {data_ += v.data_*scale;}

  void Schur(const DFT_Vec &v1, const DFT_Vec &v2) { data_ = v1.data_%v2.data_;}

  // Armadillo doesn't provide real streaming so we have to construct it.
  friend ostream &operator<<(ostream &of, const DFT_Vec &v)
  {
    unsigned N = v.data_.size();;
    of.write((char*) &N, sizeof(unsigned));
    of.write((char *)(v.data_.memptr()), N*sizeof(double));
    return of;
  }
  friend istream &operator>>(istream  &in, DFT_Vec &v )
  {
    unsigned N = 0;
    in.read((char*) &N, sizeof(unsigned));
    v.resize(N);
    in.read((char *) (v.data_.memptr()), N*sizeof(double));
    return in;
  }    


  template<class Archive> void save(Archive & ar, const unsigned int version) const
  {
    unsigned N = data_.size();
    boost::serialization::binary_object buf_wrap(data_.memptr(), N*sizeof(double));
    ar & N;
    ar & buf_wrap;
  }
  
  template<class Archive>  void load(Archive & ar, const unsigned int version)
  {
    unsigned N  = 0;
    ar & N;
    data_.resize(N);    
    boost::serialization::binary_object buf_wrap(data_.memptr(), N*sizeof(double));
    ar & buf_wrap;
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER()  


  
  // These are legacy functions that should be removed at some point. 
  void save(ofstream &of) const {data_.save(of);}
  void load(ifstream &in) {data_.load(in);}
  
 protected:
  arma::vec data_;
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
  
  void Schur(const DFT_Vec_Complex &v1, const DFT_Vec_Complex &v2, bool bUseConj=false) { if(bUseConj) data_ = v1.data_%conj(v2.data_); else data_ = v1.data_%v2.data_;}
  void incrementSchur(const DFT_Vec_Complex &v1, const DFT_Vec_Complex &v2, bool bUseConj=false) { if(bUseConj) data_ += v1.data_%conj(v2.data_); else data_ += v1.data_%v2.data_;}
 
  void resize(long N) {data_.resize(N);}
  void zeros(long N)  {data_.zeros(N);}
  void zeros()        {data_.zeros();}

  complex<double> max() const { return data_.max();}
  
  complex<double> *memptr() { return data_.memptr();}
 
  unsigned size() const { return data_.size();}


  friend ostream &operator<<(ostream &of, const DFT_Vec_Complex &v)
  {    
    unsigned N = v.data_.size();;
    of.write((char *) &N, sizeof(unsigned));
    of.write((char *)(v.data_.memptr()), N*sizeof(complex<double>));
    return of;
  }
  friend istream &operator>>(istream  &in, DFT_Vec_Complex &v )
  {
    unsigned N = 0;
    in.read((char *) &N, sizeof(unsigned));
    v.resize(N);
    in.read((char *)(v.data_.memptr()), N*sizeof(complex<double>));
    return in;
  }
  
  template<class Archive> void save(Archive & ar, const unsigned int version) const
  {
    unsigned N = data_.size();
    boost::serialization::binary_object buf_wrap(data_.memptr(), N*sizeof(complex<double>));
    ar & N;
    ar & buf_wrap;
  }
  
  template<class Archive>  void load(Archive & ar, const unsigned int version)
  {
    unsigned N  = 0;
    ar & N;
    data_.resize(N);    
    boost::serialization::binary_object buf_wrap(data_.memptr(), N*sizeof(complex<double>));
    ar & buf_wrap;
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER()  
  
protected:
  arma::cx_vec data_;
};

/**
  *  @brief UTILITY: This class encapsulates a basic FFT. It is basically an interface to fftw library while holding both a real and a complex vector.
  */  
class DFT_FFT
{
 public:
  DFT_FFT(unsigned Nx, unsigned Ny, unsigned Nz) : RealSpace_(Nx*Ny*Nz), FourierSpace_(Nx*Ny*((Nz/2)+1)), four_2_real_(NULL), real_2_four_(NULL), Nx_(Nx), Ny_(Ny), Nz_(Nz)
    {
      four_2_real_ = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, reinterpret_cast<fftw_complex*>(FourierSpace_.memptr()), RealSpace_.memptr(), FMT_FFTW);
      real_2_four_ = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, RealSpace_.memptr(),reinterpret_cast<fftw_complex*>(FourierSpace_.memptr()),  FMT_FFTW);
    };

  DFT_FFT(const DFT_FFT &Other) : RealSpace_(Other.RealSpace_), FourierSpace_(Other.FourierSpace_), four_2_real_(NULL), real_2_four_(NULL), Nx_(Other.Nx_), Ny_(Other.Ny_), Nz_(Other.Nz_)
    {
      four_2_real_ = fftw_plan_dft_c2r_3d(Nx_, Ny_, Nz_, reinterpret_cast<fftw_complex*>(FourierSpace_.memptr()), RealSpace_.memptr(), FMT_FFTW);
      real_2_four_ = fftw_plan_dft_r2c_3d(Nx_, Ny_, Nz_, RealSpace_.memptr(),reinterpret_cast<fftw_complex*>(FourierSpace_.memptr()),  FMT_FFTW);
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
    Nx_ = Nx;
    Ny_ = Ny;
    Nz_ = Nz;
  };
  
  DFT_Vec &Real() { return RealSpace_;}
  DFT_Vec_Complex &Four() { return FourierSpace_;}

  const DFT_Vec &cReal() const { return RealSpace_;}
  const DFT_Vec_Complex &cFour() const { return FourierSpace_;}

  void do_real_2_fourier() {fftw_execute(real_2_four_);}
  void do_fourier_2_real() {fftw_execute(four_2_real_);}

  void MultBy(double val) { RealSpace_.MultBy(val);  FourierSpace_.MultBy(val);}

  friend ostream &operator<<(ostream &of, const DFT_FFT &v)
  {    
    of << v.RealSpace_ << v.FourierSpace_ << v.Nx_ << " " << v.Ny_ << " " << v.Nz_; 
    return of;
  }  
  friend istream &operator>>(istream  &in, DFT_FFT &v )
  {
    in >> v.RealSpace_  >> v.FourierSpace_ >> v.Nx_ >> v.Ny_ >> v.Nz_;    

    // This is really unpleasant in terms of wasted memory:
    // fftw messes up the data so we have to restore it after making the plans
    // copy the data:
    DFT_Vec r(v.RealSpace_);
    DFT_Vec_Complex f(v.FourierSpace_);
    
    v.four_2_real_ = fftw_plan_dft_c2r_3d(v.Nx_, v.Ny_, v.Nz_, reinterpret_cast<fftw_complex*>(v.FourierSpace_.memptr()), v.RealSpace_.memptr(), FMT_FFTW);
    v.real_2_four_ = fftw_plan_dft_r2c_3d(v.Nx_, v.Ny_, v.Nz_, v.RealSpace_.memptr(),reinterpret_cast<fftw_complex*>(v.FourierSpace_.memptr()),  FMT_FFTW);

    // restore data:
    // I could just exchange the memptr's but i hate to count on library internals
    
    v.RealSpace_.set(r);
    v.FourierSpace_.set(f);

    return in;
  }

  friend class boost::serialization::access;  
  template<class Archive> void save(Archive & ar, const unsigned int version) const
  {
    ar & RealSpace_;
    ar & FourierSpace_;
    ar & Nx_;
    ar & Ny_;
    ar & Nz_;
  }
  template<class Archive>  void load(Archive & ar, const unsigned int version)
  {
    ar & RealSpace_;
    ar & FourierSpace_;
    ar & Nx_;
    ar & Ny_;
    ar & Nz_;    

    
    // I do not know how to stream the fftw plans, so I reconstruct them.
    // BUT, this is really unpleasant in terms of wasted memory:
    // fftw messes up the data so we have to restore it after making the plans
    // copy the data:
    DFT_Vec r(RealSpace_);
    DFT_Vec_Complex f(FourierSpace_);
    
    four_2_real_ = fftw_plan_dft_c2r_3d(Nx_, Ny_, Nz_, reinterpret_cast<fftw_complex*>(FourierSpace_.memptr()), RealSpace_.memptr(), FMT_FFTW);
    real_2_four_ = fftw_plan_dft_r2c_3d(Nx_, Ny_, Nz_, RealSpace_.memptr(),reinterpret_cast<fftw_complex*>(FourierSpace_.memptr()),  FMT_FFTW);

    // restore data:
    // I could just exchange the memptr's but i hate to count on library internals
    
    RealSpace_.set(r);
    FourierSpace_.set(f);
    
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER()  
  
 protected:
  DFT_Vec RealSpace_;
  DFT_Vec_Complex FourierSpace_;
  fftw_plan four_2_real_;   ///< fftw3 "plan" to perform fft on density
  fftw_plan real_2_four_;   ///< fftw3 "plan" to perform fft on density
  unsigned Nx_; ///< spatial dimension: only needed for reading and writing plans.
  unsigned Ny_; ///< spatial dimension: only needed for reading and writing plans.
  unsigned Nz_; ///< spatial dimension: only needed for reading and writing plans.
};


#endif // __DFT_LINALG_ARMA_LUTSKO__
