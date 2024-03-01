#ifndef __DFT_LINALG_LUTSKO__
#define __DFT_LINALG_LUTSKO__

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/binary_object.hpp>

#include "FMT_FFTW.h"
#include <complex>

/**
  *  @brief UTILITY: A wrapper for linear algebra packages. The idea is to be able to easily change libraries without changing any other code. Probably overkill and could be eliminated.
  */  
class DFT_Vec
{
 public:
  DFT_Vec(unsigned N);
  DFT_Vec(const DFT_Vec& c);
  DFT_Vec();
  ~DFT_Vec();

  static const char *get_library_name();
  
  DFT_Vec& operator= (const DFT_Vec& c){set(c); return *this;}
  
  void   set(unsigned pos, double val);
  double get(unsigned pos) const;

  void set(const DFT_Vec& x);
  void set(const DFT_Vec& v1, const DFT_Vec& v2, double scale); // = v1+v2*scale
  void set(const double *x, unsigned n); // copy x (length n) into data
  void set(double d); // set all values to constant
  void set_random_normal(); // set all values to random, gaussian distributed numbers with mean zero and variance one.
  void set_random(); // uniform random elements from (-1,1)
  
  void setFromAlias(const DFT_Vec &x);
  void setAliasFromValues(const DFT_Vec &x);
  void alias_Jacobian(const DFT_Vec &x);

  void resize(long N);
  void zeros(long N);
  void zeros();

  double inf_norm() const;
  double euclidean_norm() const;
  void normalise();
  
  double min() const;
  double max() const;

  double *memptr();
  const void* get_raw_data() const { return data_;} // ugly! eventually need a better solution
  unsigned size() const;
  
  double dotWith(const DFT_Vec &v) const;

  double accu() const;

  void MultBy(double val);
  void IncrementBy(const DFT_Vec& v);
  void DecrementBy(const DFT_Vec& v);
  void add(double shift);
  void IncrementBy(unsigned pos, double val);

  void operator*=(double a) { MultBy(a);}
  void operator/=(double a) { MultBy(1.0/a);}
  void operator+=(double d) { add(d);}
  void operator+=(const DFT_Vec &v) { IncrementBy(v);}
  void operator-=(const DFT_Vec &v) { DecrementBy(v);}

  void IncrementBy_Scaled_Vector(const DFT_Vec& v,double scale);

  void Schur(const DFT_Vec &v1, const DFT_Vec &v2);

  // There may be value in having library-specific implementations  of these functions
  // but for now, I do not see it. 
  friend ostream &operator<<(ostream &of, const DFT_Vec &v)
    {
      unsigned N = v.size();
      of.write((char*) &N, sizeof(unsigned));
      of.write((const char *)(const_cast<DFT_Vec&>(v).memptr()), N*sizeof(double));
      return of;
    }
  friend istream &operator>>(istream  &in, DFT_Vec &v )
    {
      unsigned N = 0;
      in.read((char*) &N, sizeof(unsigned));
      v.resize(N);
      in.read((char *) (v.memptr()), N*sizeof(double));
      return in;
    }    

  // These need to be inline ...
  template<class Archive> void save(Archive & ar, unsigned int version) const
  {
    unsigned N = size();
    boost::serialization::binary_object buf_wrap(((DFT_Vec*) this)->memptr(), N*sizeof(double));
    ar & N;
    ar & buf_wrap;
  }    
  template<class Archive>  void load(Archive & ar, unsigned int version)
  {
    unsigned N  = 0;
    ar & N;
    resize(N);    
    boost::serialization::binary_object buf_wrap(memptr(), N*sizeof(double));
    ar & buf_wrap;
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER()

  // These are legacy functions that should be removed at some point. 
  void save(ofstream &of) const;
  void load(ifstream &in);
  
 protected:
  void *data_;
};

double operator*(const DFT_Vec &v1, const DFT_Vec &v2);

/**
  *  @brief UTILITY: A wrapper for linear algebra packages. This one handles complex values.
  */  
class DFT_Vec_Complex
{
 public:
  DFT_Vec_Complex(unsigned N);
  DFT_Vec_Complex(const DFT_Vec_Complex& c);
  DFT_Vec_Complex();
  ~DFT_Vec_Complex();

  DFT_Vec_Complex& operator= (const DFT_Vec_Complex& c) {set(c); return *this;}  
  
  void   set(const DFT_Vec_Complex& c);
  void   set(unsigned pos, complex<double> val);
  complex<double> get(unsigned pos) const;

  void MultBy(double val);
  
  void Schur(const DFT_Vec_Complex &v1, const DFT_Vec_Complex &v2, bool bUseConj=false);
  void incrementSchur(const DFT_Vec_Complex &v1, const DFT_Vec_Complex &v2, bool bUseConj=false);

  void resize(long N);
  void zeros(long N);
  void zeros();

  complex<double> max() const;
  complex<double> min() const;
  
  complex<double> *memptr();
 
  unsigned size() const;


  friend ostream& operator<<(ostream &of, const DFT_Vec_Complex &v)
    {    
      unsigned N = v.size();;
      of.write((char *) &N, sizeof(unsigned));
      of.write((char *)(const_cast<DFT_Vec_Complex&>(v).memptr()), N*sizeof(complex<double>));
      return of;
    }
  friend istream& operator>>(istream  &in, DFT_Vec_Complex &v )
    {
      unsigned N = 0;
      in.read((char *) &N, sizeof(unsigned));
      v.resize(N);
      in.read((char *)(v.memptr()), N*sizeof(complex<double>));
      return in;
    }
  template<class Archive> void save(Archive & ar, const unsigned int version) const
  {
    unsigned N = size();
    boost::serialization::binary_object buf_wrap(((DFT_Vec_Complex*) this)->memptr(), N*sizeof(complex<double>));
    ar & N;
    ar & buf_wrap;
  }
  template<class Archive>  void load(Archive & ar, const unsigned int version)
  {
    unsigned N  = 0;
    ar & N;
    resize(N);    
    boost::serialization::binary_object buf_wrap(memptr(), N*sizeof(complex<double>));
    ar & buf_wrap;
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER()  
  
  protected:
  void* data_;
};

/**
  *  This class encapsulates a basic FFT. It is basically an interface to fftw library while holding both a real and a complex vector.
  *  Note that the flag is_dirty_ is meant to track whether or not the two components are mutually consistent. Since we do not have the means to tell when
  *  the real or fourier space vectors is changed (i.e. dirty) I just set it to true whenever they are accessed via non-constant references. 
  *  This is overly cautious but at least lets us avoid disastors. 
  */  
class DFT_FFT
{
 public:
  DFT_FFT(unsigned Nx, unsigned Ny, unsigned Nz) : RealSpace_(Nx*Ny*Nz), FourierSpace_(Nx*Ny*((Nz/2)+1)), four_2_real_(NULL), real_2_four_(NULL), Nx_(Nx), Ny_(Ny), Nz_(Nz)
    {
      zeros();
      four_2_real_ = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, reinterpret_cast<fftw_complex*>(FourierSpace_.memptr()), RealSpace_.memptr(), FMT_FFTW);
      real_2_four_ = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, RealSpace_.memptr(),reinterpret_cast<fftw_complex*>(FourierSpace_.memptr()),  FMT_FFTW);
    };

  DFT_FFT(const DFT_FFT &Other) : RealSpace_(Other.RealSpace_), FourierSpace_(Other.FourierSpace_), four_2_real_(NULL), real_2_four_(NULL), Nx_(Other.Nx_), Ny_(Other.Ny_), Nz_(Other.Nz_)
    {
      four_2_real_ = fftw_plan_dft_c2r_3d(Nx_, Ny_, Nz_, reinterpret_cast<fftw_complex*>(FourierSpace_.memptr()), RealSpace_.memptr(), FMT_FFTW);
      real_2_four_ = fftw_plan_dft_r2c_3d(Nx_, Ny_, Nz_, RealSpace_.memptr(),reinterpret_cast<fftw_complex*>(FourierSpace_.memptr()),  FMT_FFTW);
      is_dirty_ = Other.is_dirty_;
    };  

 DFT_FFT() : four_2_real_(NULL), real_2_four_(NULL){};

  DFT_FFT& operator= (const DFT_FFT& c)
  {
    if(four_2_real_) fftw_destroy_plan(four_2_real_);
    if(real_2_four_) fftw_destroy_plan(real_2_four_);
    
    RealSpace_    = c.RealSpace_;
    FourierSpace_ = c.FourierSpace_;
    
    Nx_ = c.Nx_;
    Ny_ = c.Ny_;
    Nz_ = c.Nz_;
    is_dirty_ = true; // mutable - so default to worse case

    // I don't know how to copy these ... so recreate
    four_2_real_ = fftw_plan_dft_c2r_3d(Nx_, Ny_, Nz_, reinterpret_cast<fftw_complex*>(FourierSpace_.memptr()), RealSpace_.memptr(), FMT_FFTW);
    real_2_four_ = fftw_plan_dft_r2c_3d(Nx_, Ny_, Nz_, RealSpace_.memptr(),reinterpret_cast<fftw_complex*>(FourierSpace_.memptr()),  FMT_FFTW);
    return *this;
  }
  
  ~DFT_FFT()
    {
      if(four_2_real_) fftw_destroy_plan(four_2_real_);
      if(real_2_four_) fftw_destroy_plan(real_2_four_);
    }

  void zeros() {RealSpace_.zeros(); FourierSpace_.zeros(); is_dirty_ = false;}


  bool get_is_dirty() const {return is_dirty_;}
  
  void initialize(unsigned Nx, unsigned Ny, unsigned Nz)
  {
    if(four_2_real_) fftw_destroy_plan(four_2_real_);
    if(real_2_four_) fftw_destroy_plan(real_2_four_);
    
    RealSpace_.zeros(Nx*Ny*Nz);
    FourierSpace_.zeros(Nx*Ny*((Nz/2)+1));
    four_2_real_ = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, reinterpret_cast<fftw_complex*>(FourierSpace_.memptr()), RealSpace_.memptr(), FMT_FFTW);
    real_2_four_ = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, RealSpace_.memptr(),reinterpret_cast<fftw_complex*>(FourierSpace_.memptr()),  FMT_FFTW);
    Nx_ = Nx;
    Ny_ = Ny;
    Nz_ = Nz;
    is_dirty_ = false;
  };

  // Assume the worse (is_dirty) with these non-constant accessors
  DFT_Vec &Real()         { is_dirty_ = true; return RealSpace_;}
  DFT_Vec_Complex &Four() { is_dirty_ = true; return FourierSpace_;}

  const DFT_Vec &cReal()         const { return RealSpace_;}
  const DFT_Vec_Complex &cFour() const { return FourierSpace_;}

  void set(double x)                  { RealSpace_.set(x); is_dirty_ = true;}
  void set(const DFT_Vec &x)          { RealSpace_.set(x); is_dirty_ = true;}
  void set(unsigned pos, double val)  { RealSpace_.set(pos,val); is_dirty_ = true;}    

  void do_real_2_fourier() {if(is_dirty_) fftw_execute(real_2_four_); is_dirty_ = false;}
  void do_fourier_2_real() {if(is_dirty_) fftw_execute(four_2_real_); is_dirty_ = false;}

  void MultBy(double val)       { RealSpace_.MultBy(val);  FourierSpace_.MultBy(val);} // doesn't change is_dirty_
  void IncrementBy(DFT_Vec & v) { RealSpace_.IncrementBy(v); is_dirty_ = true;}
  void DecrementBy(DFT_Vec & v) { RealSpace_.DecrementBy(v); is_dirty_ = true;}
  void IncrementBy(unsigned pos, double val) { RealSpace_.IncrementBy(pos,val); is_dirty_ = true;}     
  
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
    DFT_Vec         r(v.RealSpace_);
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
    DFT_Vec         r(RealSpace_);
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
  mutable bool is_dirty_ = true;
};

#endif // __DFT_LINALG_LUTSKO__

