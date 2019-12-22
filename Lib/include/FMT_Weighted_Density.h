#ifndef __LUTSKO_FMT_Weighted_Density3__
#define __LUTSKO_FMT_Weighted_Density3__


#include "DFT_LinAlg.h"


/**
  *  @brief UTILITY: abstracts the FMT weighted densities. It has vectors to hold the weight functions and the resulting fundamental measures and knows how to do convolutions via FFT.
  *         It also carries around a field dPhi_ which is d Phi/ d weighted_density_(i,j,k) where weighted_density_(i,j,k) is this weighted density. 
  */  

class FMT_Weighted_Density
{
 public:
  FMT_Weighted_Density(){ }

  FMT_Weighted_Density(const FMT_Weighted_Density &) = delete;
  
  ~FMT_Weighted_Density(){}
    
  void initialize(long Nx, long Ny, long Nz)
  {
    weighted_density_.initialize(Nx,Ny,Nz);
    weight_.initialize(Nx,Ny,Nz);
    dPhi_.initialize(Nx,Ny,Nz);
  }

  // This function does an fft of the weights from real to fourier space
  // NOTE: I throw in a factor of 1/Ntot because these are only ever used to 
  // perform convolutions and then transform back to real space and that process
  // requires a multiplication by 1/Ntot due to fftw3's normalization convention.
  void transformWeights()
  {
    weight_.do_real_2_fourier();
    weight_.Four().MultBy(1.0/weight_.Real().size());
  }

  // Tricky point: here, we have to use the conjugates because the original expression has the form
  // e[i] = sum_{j} w(i+j) rho(j)
  void convolute(const DFT_Vec_Complex& density) 
  {
    weighted_density_.Four().Schur(density,weight_.Four(),true);

    weighted_density_.do_fourier_2_real();
  }

  void extractDensity(DFT_FFT& density)
  {
    DFT_Vec cp = density.Real();
    
    density.do_real_2_fourier();
    density.do_fourier_2_real();

    density.Real().MultBy(1.0/density.Real().size());
    
    double dmax = 0.0;
    double imax = -1;
    for(long i=0;i<cp.size();i++)
      {
	double d = fabs(cp.get(i)-density.Real().get(i));
	if(d > dmax) { dmax = d; imax = i;}
      }

    cout << "imax = " << imax << " dmax = " << dmax << endl;
    cout << setprecision(20) << cp.get(imax) << endl;
    cout << setprecision(20) << density.Real().get(imax) << endl;

    for(long i=0;i<10;i++)
      cout << cp.get(i) << " " << density.Real().get(i) << " " << cp.get(i)-density.Real().get(i) << endl;



    
    cout << (1.0/weight_.Real().size())*density.Four().get(0) << " " << (1.0/weight_.Real().size())*density.Four().get(1000) << endl;
    
    weighted_density_.do_real_2_fourier();

    //    cout << "Here: " << density.Four().get(1000)*(weight_.Four().get(1000)) << " " << weighted_density_.Four().get(1000) << endl;
    cout << density.Four().get(1000).real() << " " << density.Four().get(1000).imag() << endl;
    cout << weight_.Four().get(1000).real() << " " << weight_.Four().get(1000).imag() << endl;
    cout << "Here: " << density.Four().get(1000)*weight_.Four().get(1000) << " " << weighted_density_.Four().get(1000) << endl;
    //    cout << "Here: " << density.Four().get(1000)*std::conj(weight_.Four().get(1000)) << " " << weighted_density_.Four().get(1000) << endl;

    for(long i=0;i<density.Four().size(); i++)
      density.Four().set(i, (1.0/weight_.Real().size())*(1.0/weight_.Real().size())*weighted_density_.Four().get(i)/conj(weight_.Four().get(i)));

    cout << density.Four().get(0) << " " << density.Four().get(1000) << endl;
    
    density.do_fourier_2_real();
    
  }
  
  void add_to_dPhi(DFT_Vec_Complex& dPhi)
  {
    dPhi_.do_real_2_fourier();
    dPhi.incrementSchur(dPhi_.Four(), weight_.Four());
  }
    
  void   setWeight(long pos, double x) {weight_.Real().set(pos,x);} 
  double getWeight(long pos) const {return weight_.cReal().get(pos);}
  void   addToWeight(long pos, double x) {weight_.Real().IncrementBy(pos,x);}
  void   Set_dPhi(long pos, double x) {dPhi_.Real().set(pos,x);} 
 
  double real(long i) const { return weighted_density_.cReal().get(i); }
 
  const DFT_Vec_Complex &wk()   const {return weight_.cFour();}
  const DFT_Vec         &Real() const {return weighted_density_.cReal();}
  const DFT_Vec_Complex &Four() const {return weighted_density_.cFour();}
 
  void setWk(long pos, double x, double y) {weight_.Four().set(pos, complex<double>(x,y));} 

  void dump(ofstream &of){weight_.cReal().save(of);}
  void load(ifstream &in){weight_.Real().load(in);}
  
 protected:    
  DFT_FFT weighted_density_;
  DFT_FFT weight_;
  DFT_FFT dPhi_;
};


#endif //  __LUTSKO_FMT_Weighted_Density3__ sentinal
