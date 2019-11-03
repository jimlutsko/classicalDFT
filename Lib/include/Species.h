#ifndef __LUTSKO__SPECIES__
#define __LUTSKO__SPECIES__


#include "FMT_Weighted_Density.h"
#include "Potential1.h"
#include "Density.h"

class Species
{
 public:
 Species(Density &density, double mu = 0) : density_(density), dF_(density.Ntot()), mu_(mu), fixedMass_(-1) { seq_num_ = SequenceNumber_++;}
  ~Species(){}

  int getSequenceNumber() const { return seq_num_;}

  void doFFT(){ density_.doFFT();}
  
  void setFixedMass(double m) { fixedMass_ = m; if(m > 0.0) mu_ = 0.0;}

  double getChemPotential() const {return mu_;}
  void setChemPotential(double m) {mu_ = m;}

  const Lattice& getLattice() const { return density_;}
  const Density& getDensity() const { return density_;}

  void doDisplay(string &title, string &file) const { density_.doDisplay(title,file, seq_num_);}
  void set_density_from_amplitude(DFT_Vec &x) {density_.set_density_from_amplitude(x);}

  void zeroForce() {dF_.zeros();}
  void addToForce(long i, double v) {dF_.IncrementBy(i,v);}
  void addToForce(const DFT_Vec &f) {dF_.IncrementBy(f);}
  void setForce(const DFT_Vec &f) {dF_.set(f);}
  void multForce(double f) {dF_.MultBy(f);}  
  
  double get_convergence_monitor() const { return dF_.inf_norm()/density_.dV();}
  
  DFT_Vec &getDF() {return dF_;}
  void setDF(DFT_Vec &df) {return dF_.set(df);}
  
  double externalField(bool bCalcForces)
  {
    double dV = density_.dV();
    double Fx = density_.getExternalFieldEnergy()*dV - density_.getNumberAtoms()*mu_;
    if(bCalcForces)
      {
	dF_.IncrementBy_Scaled_Vector(density_.getExternalField(),dV);
	dF_.ShiftBy(mu_*dV);
      }
    return Fx;
  }

  virtual double getHSD() const { return 0.0;}

  void beginForceCalculation()
  {
    if(fixedMass_ > 0.0)
      {
	density_.scale_to(fixedMass_/density_.getNumberAtoms());
	mu_ = 0.0;
      }
  }

  void endForceCalculation()
  {
    if(fixedMass_ > 0.0)
      {
	mu_ = 0.0;
      
	for(long p=0;p<density_.Ntot();p++)
	  mu_ += dF_.get(p)*density_.getDensity(p);
	mu_ *= density_.dV()/fixedMass_;
	  
	for(long p=0;p<density_.Ntot();p++)
	  dF_.set(p, dF_.get(p)-mu_);
      }
  }
  
 private:
  static int SequenceNumber_;
  
 protected:
  Density &density_;
  DFT_Vec dF_;
  double mu_;
  int seq_num_;
  double fixedMass_; // if this is > 0, then the mass of this species is being held fixed at this value.
};

  /**
   *  @brief Species Class: hard sphere diameters, etc.
   *
   *    This class holds an 11-dimentional array called d_ . This in turn holds the weighted densities as
   *             d_ = {eta(N),s(N),V1(N),...,VD(N), T11(N), T12(N), ... T1D(N), T22(N), ..., TDD(N)}
   *             where for D=3 there are 1+1+3+(3+2+1) = 11 entries.
   *             Note that each of these, e.g. eta(N), is an N-dimensional vector holding the value of the weighted density at each point.
   */

class FMT_Species : public Species
{
public:
  /**
   *   @brief  Default  constructor for FMT_Species 
   *  
   *   @param  hsd is the hard-sphere diameter
   *   @param  lattice describes the mesh
   *   @param  pointsFile contains the points for spherical integration
   *   @return nothing 
   */    
  FMT_Species(Density& density, double hsd, string &pointsFile);

  FMT_Species(const FMT_Species &) = delete;
  
  ~FMT_Species(){}

  /**
   *   @brief  Accessor for hard sphere diameter
   *  
   *   @param  none
   *   @return hsd_
   */      
  virtual double getHSD() const { return hsd_;}

  /**
   *   @brief  get value of Eta at pos
   *  
   *   @param  pos is the mesh position
   *   @return value of Eta at pos
   */    
  double getEta(long pos)           const { return d_[EI()].real(pos);}

  /**
   *   @brief  get value of S at pos
   *  
   *   @param  pos is the mesh position
   *   @return value of S at pos
   */      
  double getS(long pos)             const { return d_[SI()].real(pos);}

  /**
   *   @brief  get value of component j of V at pos
   *  
   *   @param  j is index of V
   *   @param  pos is the mesh position
   *   @return value of V(j) at pos
   */      
  double getV(int j, long pos)      const { return d_[VI(j)].real(pos);}

  /**
   *   @brief  get value of component j,k of T at pos
   *  
   *   @param  j is first index of T
   *   @param  k is second index of T
   *   @param  pos is the mesh position
   *   @return value of T(j,k) at pos
   */      
  double getT(int j,int k,long pos) const { return d_[TI(j,k)].real(pos);}
  
  const DFT_Vec_Complex& getWEK() const { return d_[EI()].wk();}

  /**
   *   @brief  set value of dPhi_dEta at pos
   *  
   *   @param  pos is the mesh position
   *   @return none
   */        
  void Set_dPhi_Eta(long i,double val)             {d_[EI()].Set_dPhi(i,val);}

  /**
   *   @brief  set value of dPhi_dS at pos
   *  
   *   @param  pos is the mesh position
   *   @return none
   */        
  void Set_dPhi_S(long i,double val)               {d_[SI()].Set_dPhi(i,val);}

  /**
   *   @brief  set value of dPhi_dV_j at pos
   *  
   *   @param  j is index of V
   *   @param  pos is the mesh position
   *   @return none
   */        
  void Set_dPhi_V(int j, long i,double val)        {d_[VI(j)].Set_dPhi(i,val);}

  /**
   *   @brief  set value of dPhi_dT(j,k) at pos
   *  
   *   @param  j is first index of T
   *   @param  k is second index of T
   *   @param  pos is the mesh position
   *   @return none
   */        
  void Set_dPhi_T(int j, int k, long i,double val) {d_[TI(j,k)].Set_dPhi(i,val);}

  /**
   *   @brief This does the convolution of the density and the weight for each weighted density after which it converts back to real space 
   *          ( so this computes the weighted densities n(r) = int w(r-r')rho(r')dr'). The results are all stored in parts of FMT_Weighted_Density
   *  
   *   @return none
   */        
  void convoluteDensities(bool needsTensor)
  {
    // reference to Fourier-space array of density
    density_.doFFT();
    const DFT_Vec_Complex &rho_k = density_.getDK();

    int imax = (needsTensor ? d_.size() : 5);

    for(int i=0;i<imax;i++)
      d_[i].convolute(rho_k);      
  }

  /**
   *   @brief Loop over the weighted densities and ask each one to add its contribution to dPhi
   *          In other words:   SUM_{a} d PHI/d n_{a}
   *  
   *   @return none
   */  
  void Accumulate_dPhi(DFT_Vec_Complex& dPhi, bool needsTensor)
  {
    int imax = (needsTensor ? d_.size() : 5);

    for(int i=0;i<imax;i++)      
      d_[i].add_to_dPhi(dPhi);
  }

  // Used in DFT_Surfactant ...
  const DFT_Vec &getV_Real(int J) const { return d_[VI(J)].Real();}
  const DFT_Vec_Complex& getVweight_Four(int J) const { return d_[VI(J)].wk();}  
  
protected:
  /**
   *   @brief  This is a one-time-only evaluation of the numerical approximation to the FMT weight functions. These are all 
   *           functions w_{alpha}(i,j) = w_{alpha}(abs(i-j)). Their evaluation involves real-space integrations for which the 
   *           integration points are given in the file pointsFile. Most of the work occurs via a call to initializeWeightedDensities@
   *
   *   @param  densities: the array of weighted densities
   *   @param  hsd is the hard-sphere diameter
   *   @param  pointsFile is the file holding the integration points for integrating a spherical shell
   *   @param  Nx is the number of lattice points in the x-direction
   *   @param  Ny is the number of lattice points in the y-direction
   *   @param  Nz is the number of lattice points in the z-direction
   */        
  void generateWeights(string &pointsFile);

  /**
   *   @brief  Get the index of the "eta" partial weighted density in the array of weighted densities, d_
   *   @returns  the index.
   */        
  int EI() const {return 0;}

  /**
   *   @brief  Get the index of the scalar partial weighted density in the array of weighted densities, d_
   *   @returns  the index.
   */          
  int SI() const {return 1;}
  /**
   *   @brief  Get the index of the vector partial weighted density in the array of weighted densities, d_
   *   @returns  the index.
   */          
  int VI(int j) const {return 2+j;}
  /**
   *   @brief  Get the index of the tensor partial weighted density in the array of weighted densities, d_
   *   @returns  the index.
   */          
  int TI(int j, int k) const
  {
    if(j > k) swap(j,k);
    if(j == 0) return 5+k;
    else if (j == 1) return 7+k;
    return 10;
  }

protected:
  double hsd_ = 0.0; ///< hard sphere diameter 
  vector<FMT_Weighted_Density>  d_; ///< all weighted densities in real & fourier space
};

/*
  class VDW_Species : public FMT_Species
  {
  public:
    VDW_Species(Density& density, double hsd, string &pointsFile, Potential1& potential, double kT);

    double get_VDW_Constant() const {return a_vdw_;}

    double getInteractionEnergyAndForces()
    {
      long Ntot = density_.Ntot();
      double dV = density_.dV();

      DFT_FFT v(density_.Nx(), density_.Ny(), density_.Nz());      

      v.Four().Schur(density_.getDK(),w_att_.Four());
      v.Four().MultBy(dV*dV/Ntot);
      v.do_fourier_2_real(); 

      dF_.IncrementBy(v.Real());
    
      return 0.5*density_.getInteractionEnergy(v.Real());
    }
  
  protected:
    Potential1& potential_;
    double a_vdw_;
    DFT_FFT w_att_;
  };
*/
#endif // __LUTSKO__SPECIES__
