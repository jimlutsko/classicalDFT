#ifndef __LUTSKO__SPECIES__
#define __LUTSKO__SPECIES__

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/export.hpp>

#include "FMT_Weighted_Density.h"
#include "Potential1.h"
#include "Density.h"
#include "Fundamental_Measures.h"

class Species
{
 public:
  Species(Density &density, double mu = 0, int seq = -1) : density_(density), dF_(density.Ntot()), mu_(mu), fixedMass_(-1) { if(seq >=0) {seq_num_ = seq; SequenceNumber_ = seq+1;} else seq_num_ = SequenceNumber_++;}
  ~Species(){}

  int getSequenceNumber() const { return seq_num_;}

  void doFFT(){ density_.doFFT();}
  
  void setFixedMass(double m)            { fixedMass_          = m; if(m > 0.0) mu_ = 0.0;}
  void setFixedBackground(bool fixed)    { fixedBackground_    = fixed;}
  void setHomogeneousBoundary(bool homo) { homgeneousBoundary_ = homo;}
  
  bool is_background_fixed() const { return fixedBackground_;}
  bool is_mass_fixed()       const { return (fixedMass_ > 0);}
  bool is_fixed_boundary() const { return   fixedBackground_;}
  
  double getChemPotential() const {return mu_;}
  void   setChemPotential(double m) {mu_ = m;}

  const Lattice& getLattice() const { return density_;}
  const Density& getDensity() const { return density_;}
  const double*  get_density_data() { return density_.get_density_pointer();}
  virtual void   get_density_alias(DFT_Vec &x) const;

  virtual void convert_to_alias_deriv(DFT_Vec &x, DFT_Vec &dF_dRho) const;
  
  void doDisplay(string &title, string &file, void *param = NULL) const { density_.doDisplay(title,file, seq_num_, param);}

  void         set_density(DFT_Vec &x) {density_.set(x); density_.doFFT();}
  void         set_density(long j, double x) {density_.set(j,x);}
  virtual void set_density_from_alias(const DFT_Vec &x);

  void fft_density() { density_.doFFT();}
  
  void zeroForce() {dF_.zeros();}
  void addToForce(long i, double v) {dF_.IncrementBy(i,v);}
  void addToForce(const DFT_Vec &f) {dF_.IncrementBy(f);}
  void setForce(const DFT_Vec &f) {dF_.set(f);}
  void multForce(double f) {dF_.MultBy(f);}  
  
  double get_convergence_monitor() const { return dF_.inf_norm()/density_.dV();}
  
  DFT_Vec       &getDF() {return dF_;}
  const DFT_Vec &get_const_DF() const {return dF_;}
  void setDF(DFT_Vec &df) {return dF_.set(df);}

  // This species is held as allSpecies_[index_]
  void setIndex(int i)  { index_ = i;}
  int  getIndex() const { return index_;}
  
  /**
  *   @brief  This adds in the external field contribution (which is species dependent). Force calculation is optional
  *  
  *   @param  bCalcForces: forces are calculated if this is true
  *   @return  the contribution to the free energy
  */  
  double externalField(bool bCalcForces)
  {
    double dV = density_.dV();
    double Fx = density_.get_field_dot_density()*dV - density_.getNumberAtoms()*mu_;
    if(bCalcForces)
      {
	dF_.IncrementBy_Scaled_Vector(density_.get_external_field(),dV);
	dF_.ShiftBy(mu_*dV);
      }
    return Fx;
  }

  virtual double getHSD() const { return 0.0;}
  
  // Placeholder for FMT-specific processing: non-FMT classes do nothing
  virtual void set_fundamental_measure_derivatives(FundamentalMeasures &fm, long pos, bool needsTensor) {}
 
  //Constant particle number is enforced at the species-level. If needed, some information has to be collected before updating the forces. Note that particle number is rigorously kept constant.
  void beginForceCalculation()
  {
    if(fixedMass_ > 0.0)
      {
	density_.scale_to(fixedMass_/density_.getNumberAtoms());
	mu_ = 0.0;
      }
  }

  // Constant particle number is enforced at the species-level. If activated, the necessary corrections to the forces are applied here. Note that particle number is rigorously kept constant.
  double endForceCalculation()
  {
    if(fixedBackground_ && fixedMass_ > 0.0)
      throw std::runtime_error("Cannot have both fixed background and fixed mass .... aborting");
    
    
    if(fixedBackground_)
      {
	for(long pos = 0; pos < density_.get_Nboundary(); pos++)
	  dF_.set(density_.boundary_pos_2_pos(pos),0.0);	
      }

    if(homgeneousBoundary_)
      {
	
	double average_border_force = 0;

	for(long pos = 0; pos < density_.get_Nboundary(); pos++)
	  average_border_force += dF_.get(density_.boundary_pos_2_pos(pos));

	average_border_force /= density_.get_Nboundary();

	for(long pos = 0; pos < density_.get_Nboundary(); pos++)
	  dF_.set(density_.boundary_pos_2_pos(pos),average_border_force);
      }    

    
    if(fixedMass_ > 0.0)
      {
	mu_ = 0.0;

	double Mtarget = fixedMass_;

	for(long p=0;p<density_.Ntot();p++)
	  mu_ += dF_.get(p)*density_.get(p);
	mu_ /= Mtarget; //fixedMass_;
	for(long p=0;p<density_.Ntot();p++)
	  dF_.set(p, dF_.get(p)-mu_*density_.dV());
      }
    return 0;
  }

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive &ar, const unsigned int file_version){}
  
  template<class Archive> friend void boost::serialization::save_construct_data(Archive & ar, const Species * t, const unsigned int file_version);
  template<class Archive> friend void boost::serialization::load_construct_data(Archive & ar, Species * t, const unsigned int file_version);

  
protected:
 
  static int SequenceNumber_; ///< Give every species a unique identifier. This generates a new id for each new instance.
  
 protected:
  Density &density_;
  DFT_Vec dF_;
  double mu_;
  int seq_num_;
  double fixedMass_        = -1; // if this is > 0, then the mass of this species is being held fixed at this value.
  bool fixedBackground_    = false; // if true, forces are set to zero on the background
  bool homgeneousBoundary_ = false; 
  int index_ = -1;
};


  /**
   *  @brief Species Class: hard sphere diameters, etc.
   *
   *    This class holds an 11-dimentional array called fmt_weighted_densities . This in turn holds the weighted densities as
   *             fmt_weighted_densities = {eta(N),s(N),V1(N),...,VD(N), T11(N), T12(N), ... T1D(N), T22(N), ..., TDD(N)}
   *             where for D=3 there are 1+1+3+(3+2+1) = 11 entries.
   *             Note that each of these, e.g. eta(N), is an N-dimensional vector holding the value of the weighted density at each point.
   */

class FMT_Species : public Species
{
public:
  FMT_Species(Density& density, double hsd, double mu = 0, int seq = -1);
  FMT_Species(const FMT_Species &) = delete;
  ~FMT_Species(){}

  virtual double getHSD() const { return hsd_;}
  double getEta(long pos) const { return fmt_weighted_densities[EI()].real(pos);}
  double getS(long pos) const { return fmt_weighted_densities[SI()].real(pos);}
  double getV(int j, long pos)      const { return fmt_weighted_densities[VI(j)].real(pos);}
  double getT(int j,int k,long pos) const { return fmt_weighted_densities[TI(j,k)].real(pos);}

  virtual void set_density_from_alias(const DFT_Vec &x);
  virtual void get_density_alias(DFT_Vec &x) const;
  virtual void convert_to_alias_deriv(DFT_Vec &x, DFT_Vec &dF_dRho) const;
  
  const DFT_Vec_Complex& getWEK() const { return fmt_weighted_densities[EI()].wk();}

  /**
   *   @brief This does the convolution of the density and the weight for each weighted density after which it converts back to real space 
   *          ( so this computes the weighted densities n(r) = int w(r-r')rho(r')dr'). The results are all stored in parts of FMT_Weighted_Density
   *  
   *   @return none
   */        
  virtual void calculateFundamentalMeasures(bool needsTensor)
  {
    // reference to Fourier-space array of density
    density_.doFFT();
    const DFT_Vec_Complex &rho_k = density_.get_density_fourier();

    int imax = (needsTensor ? fmt_weighted_densities.size() : 5);

    for(int i=0;i<imax;i++)
      fmt_weighted_densities[i].convoluteWith(rho_k);      
  }

  void convolute_eta_weight_with(const DFT_FFT &v, DFT_FFT &result, bool bConjugate = false) const
  { convolute_weight_with(EI(),v,result,bConjugate);}

  void convolute_s_weight_with(const DFT_FFT &v, DFT_FFT &result, bool bConjugate = false) const
  { convolute_weight_with(SI(),v,result,bConjugate);}

  void convolute_v_weight_with(int i, const DFT_FFT &v, DFT_FFT &result, bool bConjugate = false) const
  { convolute_weight_with(VI(i),v,result, bConjugate);}

  void convolute_T_weight_with(int i, int j, const DFT_FFT &v, DFT_FFT &result, bool bConjugate = false) const
  { convolute_weight_with(TI(i,j),v,result, bConjugate);}    

  void convolute_weight_with(int pos, const DFT_FFT &v, DFT_FFT &result, bool bConjugate = false) const
  {
    if(v.get_is_dirty() == true) throw std::runtime_error("Need real and fourier parts of input to FMT_Species::convolute_weight_with(...) to be consistent");

    result.Four().Schur(v.cFour(), fmt_weighted_densities[pos].wk(),bConjugate);
    result.do_fourier_2_real();
  }  
  
  /**
   *   @brief Loop over the weighted densities and ask each one to add its contribution to dPhi
   *          In other words:   SUM_{a} SUM_j d PHI/d n_{a}(j) w_{a}(j-i)
   *  
   *   @return none
   */  
  void Accumulate_dPhi(DFT_Vec_Complex& dPhi, bool needsTensor)
  {
    int imax = (needsTensor ? fmt_weighted_densities.size() : 5);

    for(int i=0;i<imax;i++)      
      fmt_weighted_densities[i].add_weight_schur_dPhi_to_arg(dPhi);
  }
  
  // These return the real weight at position K using the extended notation: eta, s0,s1,s2,v1,v2
  double getExtendedWeight(long K, int a)
  {
    if(a == 0) return fmt_weighted_densities[EI()].getWeight(K);

    if(a == 1) return fmt_weighted_densities[SI()].getWeight(K)/(hsd_*hsd_);
    if(a == 2) return fmt_weighted_densities[SI()].getWeight(K)/hsd_;
    if(a == 3) return fmt_weighted_densities[SI()].getWeight(K);

    if(a == 4) return fmt_weighted_densities[VI(0)].getWeight(K)/hsd_;
    if(a == 5) return fmt_weighted_densities[VI(1)].getWeight(K)/hsd_;
    if(a == 6) return fmt_weighted_densities[VI(2)].getWeight(K)/hsd_;

    if(a == 7) return fmt_weighted_densities[VI(0)].getWeight(K);
    if(a == 8) return fmt_weighted_densities[VI(1)].getWeight(K);
    if(a == 9) return fmt_weighted_densities[VI(2)].getWeight(K);

    if(a == 10) return fmt_weighted_densities[TI(0,0)].getWeight(K);
    if(a == 11) return fmt_weighted_densities[TI(0,1)].getWeight(K);
    if(a == 12) return fmt_weighted_densities[TI(0,2)].getWeight(K);

    if(a == 13) return fmt_weighted_densities[TI(1,0)].getWeight(K);
    if(a == 14) return fmt_weighted_densities[TI(1,1)].getWeight(K);
    if(a == 15) return fmt_weighted_densities[TI(1,2)].getWeight(K);

    if(a == 16) return fmt_weighted_densities[TI(2,0)].getWeight(K);
    if(a == 17) return fmt_weighted_densities[TI(2,1)].getWeight(K);
    if(a == 18) return fmt_weighted_densities[TI(2,2)].getWeight(K);        

    throw std::runtime_error("Unknown index in FMT_Weighted_Density::getExtendedWeight");
  }

  
  // These return the weighted density at position K using the extended notation: eta, s0,s1,s2,v1,v2
  double getExtendedWeightedDensity(long K, int a)
  {
    if(a == 0) return fmt_weighted_densities[EI()].getDensity(K);

    if(a == 1) return fmt_weighted_densities[SI()].getDensity(K)/(hsd_*hsd_);
    if(a == 2) return fmt_weighted_densities[SI()].getDensity(K)/hsd_;
    if(a == 3) return fmt_weighted_densities[SI()].getDensity(K);

    if(a == 4) return fmt_weighted_densities[VI(0)].getDensity(K)/hsd_;
    if(a == 5) return fmt_weighted_densities[VI(1)].getDensity(K)/hsd_;
    if(a == 6) return fmt_weighted_densities[VI(2)].getDensity(K)/hsd_;

    if(a == 7) return fmt_weighted_densities[VI(0)].getDensity(K);
    if(a == 8) return fmt_weighted_densities[VI(1)].getDensity(K);
    if(a == 9) return fmt_weighted_densities[VI(2)].getDensity(K);

    if(a == 10) return fmt_weighted_densities[TI(0,0)].getDensity(K);
    if(a == 11) return fmt_weighted_densities[TI(0,1)].getDensity(K);
    if(a == 12) return fmt_weighted_densities[TI(0,2)].getDensity(K);

    if(a == 13) return fmt_weighted_densities[TI(1,0)].getDensity(K);
    if(a == 14) return fmt_weighted_densities[TI(1,1)].getDensity(K);
    if(a == 15) return fmt_weighted_densities[TI(1,2)].getDensity(K);

    if(a == 16) return fmt_weighted_densities[TI(2,0)].getDensity(K);
    if(a == 17) return fmt_weighted_densities[TI(2,1)].getDensity(K);
    if(a == 18) return fmt_weighted_densities[TI(2,2)].getDensity(K);    

    throw std::runtime_error("Unknown index in FMT_Weighted_Density::getExtendedWeightedDensity");
  }
  


  
  // FOr testing only: brute-force evaluation of weighted density at position K using the extended notation: eta, s0,s1,s2,v1,v2
  double getBruteForceWeightedDensity(int K[3], int a)
  {
    double d = 0.0;

    int Nx = density_.Nx();
    int Ny = density_.Ny();
    int Nz = density_.Nz();    

    for(int ix = 0;ix<Nx;ix++)
      for(int iy = 0;iy<Ny;iy++)
	for(int iz = 0;iz<Nz;iz++)
	  {
	    long KI = density_.get_PBC_Pos(K[0]-ix,K[1]-iy,K[2]-iz);
	    d += getExtendedWeight(KI,a)*density_.get(ix,iy,iz);
	  }
    return d;
  }
  // These return the weighted density at position K using the extended notation: eta, s0,s1,s2,v1,v2
  void getFundamentalMeasures(long K, FundamentalMeasures &fm) {getFundamentalMeasures_Helper(K,fm,fmt_weighted_densities, hsd_);}
  void getFundamentalMeasures_Helper(long K, FundamentalMeasures &fm, const vector<FMT_Weighted_Density>  &weighted_densities, double hsd) 
  {
    fm.eta = weighted_densities[EI()].getDensity(K);

    fm.s0 = weighted_densities[SI()].getDensity(K)/(hsd*hsd);
    fm.s1 = weighted_densities[SI()].getDensity(K)/hsd;
    fm.s2 = weighted_densities[SI()].getDensity(K);

    fm.v1[0] = weighted_densities[VI(0)].getDensity(K)/hsd;
    fm.v1[1] = weighted_densities[VI(1)].getDensity(K)/hsd;
    fm.v1[2] = weighted_densities[VI(2)].getDensity(K)/hsd;

    fm.v2[0] = weighted_densities[VI(0)].getDensity(K);
    fm.v2[1] = weighted_densities[VI(1)].getDensity(K);
    fm.v2[2] = weighted_densities[VI(2)].getDensity(K);

    fm.T[0][0] = weighted_densities[TI(0,0)].getDensity(K);
    fm.T[0][1] = weighted_densities[TI(0,1)].getDensity(K);
    fm.T[0][2] = weighted_densities[TI(0,2)].getDensity(K);
    
    fm.T[1][0] = weighted_densities[TI(1,0)].getDensity(K);
    fm.T[1][1] = weighted_densities[TI(1,1)].getDensity(K);
    fm.T[1][2] = weighted_densities[TI(1,2)].getDensity(K);

    fm.T[2][0] = weighted_densities[TI(2,0)].getDensity(K);
    fm.T[2][1] = weighted_densities[TI(2,1)].getDensity(K);
    fm.T[2][2] = weighted_densities[TI(2,2)].getDensity(K);

    fm.calculate_derived_quantities();      
    
  }
  
  /**
   *   @brief  Take derivatives of free energy wrt fundamental measures and set dPhi_dEta, etc. The point is that 
   *           these are, e.g., dPhi_dS2 = dPhi_dS0/(hsd*hsd) + dPhi_dS1/(hsd) + dPhi_dS2 and so differ for each species. 
   *           An alternative would be to hold each measure individually (S0, S1, S2) but this would cost more in cpu and memory (I think).
   *
   *   @param DPHI: array holding dPhi/d n_{a}(pot) where n_{a}(pot) is fundamental measure a at position pos.
   *   @param pos: spatial point (in 1D indexing)
   *   @param needsTensor: true if we need to calculate tensor quantities too.
   *  
   */  
  virtual void set_fundamental_measure_derivatives(FundamentalMeasures &DPHI, long pos, bool needsTensor)
  {
    double dPhi_dEta = DPHI.eta;
    double dPhi_dS = (DPHI.s0/(hsd_*hsd_)) + (DPHI.s1/hsd_) + DPHI.s2;
    double dPhi_dV[3] = {DPHI.v2[0] + DPHI.v1[0]/hsd_,
			  DPHI.v2[1] + DPHI.v1[1]/hsd_,
			  DPHI.v2[2] + DPHI.v1[2]/hsd_};

    fmt_weighted_densities[EI()].Set_dPhi(pos,dPhi_dEta);
    fmt_weighted_densities[SI()].Set_dPhi(pos,dPhi_dS);    

    // Note the parity factor in the vector term which is needed when we calculate forces
    for(int j=0;j<3;j++)
      {
	fmt_weighted_densities[VI(j)].Set_dPhi(pos, -dPhi_dV[j]);	
	if(needsTensor)
	  for(int k=j;k<3;k++)
	    fmt_weighted_densities[TI(j,k)].Set_dPhi(pos,(j == k ? 1 : 2)*DPHI.T[j][k]); // taking account that we only use half the entries
      }
  }

    
  FMT_Weighted_Density& getEta() { return fmt_weighted_densities[0];}

  double getWeight(int index, long pos) { return fmt_weighted_densities[index].getWeight(pos);}
  
  // Used in DFT_Surfactant ...
  const DFT_Vec &getV_Real(int J) const { return fmt_weighted_densities[VI(J)].Real();}
  const DFT_Vec_Complex& getVweight_Four(int J) const { return fmt_weighted_densities[VI(J)].wk();}  

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive &ar, const unsigned int file_version)
  {
    boost::serialization::void_cast_register<FMT_Species, Species>(static_cast<FMT_Species *>(NULL),static_cast<Species *>(NULL));
  }    
  template<class Archive> friend void boost::serialization::save_construct_data(Archive & ar, const FMT_Species * t, const unsigned int file_version);
  template<class Archive> friend void boost::serialization::load_construct_data(Archive & ar, FMT_Species * t, const unsigned int file_version);

  
protected:
  // Indices of eta, scaler, vector and tensor weighted densities
  int EI() const {return 0;}
  int SI() const {return 1;}
  int VI(int j) const {return 2+j;}
  int TI(int j, int k) const
  {
    if(j > k) swap(j,k);
    if(j == 0) return 5+k;
    else if (j == 1) return 7+k;
    return 10;
  }

protected:
  /**
   *   @brief  This is a one-time-only evaluation of the numerical approximation to the FMT weight functions. These are all 
   *           functions w_{alpha}(i,j) = w_{alpha}(abs(i-j)). 
   *
   */        
  virtual void generateWeights(double hsd, vector<FMT_Weighted_Density> &fmt_weights);

protected:
  double hsd_ = 0.0; ///< hard sphere diameter 
  vector<FMT_Weighted_Density>  fmt_weighted_densities; ///< all weighted densities in real & fourier space
};

  /**
   *  @brief Class to implement AO model
   *
   */

class FMT_AO_Species : public FMT_Species
{
public:
  /**
   *   @brief  Default  constructor for FMT_Species 
   *  
   *   @param  hsd is the hard-sphere diameter
   *   @param  Rp is the polymer radius
   *   @param  lambda_p is the polymer-density related parameter.
   *   @return nothing 
   */    
  FMT_AO_Species(Density& density, double hsd, double Rp, double reservoir_density, double mu = 0, int seq = -1);
  FMT_AO_Species(const FMT_Species &) = delete;
  ~FMT_AO_Species(){}

  virtual void set_fundamental_measure_derivatives(FundamentalMeasures &DPHI, long pos, bool needsTensor);
  virtual double free_energy_post_process(bool needsTensor);

  unsigned size() const { return fmt_weighted_densitiesAO_.size();}
  void setPSI(long pos, double val) { PSI_.Real().set(pos,val);}
  double getReservoirDensity() const {return reservoir_density_;}
  double getHSDP() const { return 2*Rp_;}
  
  void getFundamentalMeasures_AO(long K, FundamentalMeasures &fm) {getFundamentalMeasures_Helper(K,fm,fmt_weighted_densitiesAO_, 2*Rp_);}

  
  void setFundamentalMeasures_AO(long K, double eta, double s, const double v[3], const double T[3][3])
  {
    fmt_weighted_densitiesAO_[EI()].setDensity(K,eta);
    fmt_weighted_densitiesAO_[SI()].setDensity(K,s);

    fmt_weighted_densitiesAO_[VI(0)].setDensity(K,v[0]);
    fmt_weighted_densitiesAO_[VI(1)].setDensity(K,v[1]);
    fmt_weighted_densitiesAO_[VI(2)].setDensity(K,v[2]);

    fmt_weighted_densitiesAO_[TI(0,0)].setDensity(K,T[0][0]);
    fmt_weighted_densitiesAO_[TI(0,1)].setDensity(K,T[0][1]);
    fmt_weighted_densitiesAO_[TI(0,2)].setDensity(K,T[0][2]);

    fmt_weighted_densitiesAO_[TI(1,0)].setDensity(K,T[1][0]);
    fmt_weighted_densitiesAO_[TI(1,1)].setDensity(K,T[1][1]);
    fmt_weighted_densitiesAO_[TI(1,2)].setDensity(K,T[1][2]);

    fmt_weighted_densitiesAO_[TI(2,0)].setDensity(K,T[2][0]);
    fmt_weighted_densitiesAO_[TI(2,1)].setDensity(K,T[2][1]);
    fmt_weighted_densitiesAO_[TI(2,2)].setDensity(K,T[2][2]);        
  }
  
  void computeAOForceContribution()
  {
    
    for(int a=0;a<size();a++)
      {
	fmt_weighted_densitiesAO_[a].density_do_real_2_fourier();

	// This does the Schur product of the fmtr weighted density and the AO weighted density (which is actually Upsilon-bar), putting the result into PSI and transforming PSI back to real space
	fmt_weighted_densities[a].convoluteWith(fmt_weighted_densitiesAO_[a].Four(), PSI_); 

	PSI_.Real().MultBy(reservoir_density_*density_.dV()); // because all forces are multiplied by dV
	
	addToForce(PSI_.Real());
      }
    
  }

  // TODO:
  /*
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive &ar, const unsigned int file_version)
  {
    boost::serialization::void_cast_register<FMT_Species, Species>(static_cast<FMT_AO_Species *>(NULL),static_cast<Species *>(NULL));
  }    
  template<class Archive> friend void boost::serialization::save_construct_data(Archive & ar, const FMT_AO_Species * t, const unsigned int file_version);
  template<class Archive> friend void boost::serialization::load_construct_data(Archive & ar, FMT_AO_Species * t, const unsigned int file_version);
  */
  
protected:
  double Rp_ = -1; 
  double reservoir_density_ = 0.0;
  vector<FMT_Weighted_Density>  fmt_weighted_densitiesAO_; ///< all weighted densities in real & fourier space
  DFT_FFT PSI_;
};




  /**
   *  @brief Extend FMT_Species to include EOS correction
   *
   */

class FMT_Species_EOS : public FMT_Species
{
public:
  /**
   *   @brief  Default  constructor for FMT_Species 
   *  
   *   @param  hsd is the hard-sphere diameter
   *   @param  lattice describes the mesh
   *   @return nothing 
   */    
  FMT_Species_EOS(double D_EOS, Density& density, double hsd, double mu = 0, int seq = -1);

  FMT_Species_EOS(const FMT_Species &) = delete;
  
  ~FMT_Species_EOS(){}

  virtual void calculateFundamentalMeasures(bool needsTensor)
  {
    FMT_Species::calculateFundamentalMeasures(needsTensor);

    const DFT_Vec_Complex &rho_k = density_.get_density_fourier();    
    eos_weighted_density_[0].convoluteWith(rho_k);          
  }

  double get_eos_measure(long pos) const { return eos_weighted_density_[0].real(pos);}
  
  // TODO
  //  friend class boost::serialization::access;
  //  template<class Archive> void serialize(Archive &ar, const unsigned int file_version)
  //  {
  //    boost::serialization::void_cast_register<FMT_Species, Species>(static_cast<FMT_Species *>(NULL),static_cast<Species *>(NULL));
  //  }    
  //  template<class Archive> friend void boost::serialization::save_construct_data(Archive & ar, const FMT_Species * t, const unsigned int file_version);
  //  template<class Archive> friend void boost::serialization::load_construct_data(Archive & ar, FMT_Species * t, const unsigned int file_version);
  
protected:
  //  virtual void generate_additional_Weight();
  vector<FMT_Weighted_Density> eos_weighted_density_; ///< all weighted densities in real & fourier space
  double D_EOS_;
};





template<class Archive>
inline void boost::serialization::save_construct_data(Archive & ar, const Species * t, const unsigned int file_version)
{
  ar << & t->density_;
  ar << t->mu_;
  ar << t->seq_num_;
  ar << t->dF_;
  ar << t->fixedMass_;
  ar << t->SequenceNumber_;
  ar << t->index_;
}

template<class Archive>
inline void boost::serialization::load_construct_data(Archive & ar, Species * t, const unsigned int file_version)
{
    // retrieve data from archive required to construct new instance
  Density *d;
  ar >> d;

  // invoke inplace constructor to initialize instance of my_class
  double mu = 0;
  double seq = 0;
  ::new(t)Species(*d,mu,seq);

  ar >> t->mu_;
  ar >> t->seq_num_;  
  ar >> t->dF_;
  ar >> t->fixedMass_;
  ar >> t->SequenceNumber_;
  ar >> t->index_;  
}


template<class Archive>
inline void boost::serialization::save_construct_data(Archive & ar, const FMT_Species * t, const unsigned int file_version)
{
  //  ar << static_cast<const Species*>(t);
  ar << & t->density_;
  ar << t->mu_;
  ar << t->seq_num_;
  ar << t->dF_;
  ar << t->fixedMass_;
  ar << t->SequenceNumber_;
  ar << t->index_;  
  ar << t->hsd_;
  ar << t->fmt_weighted_densities;
}

template<class Archive>
inline void boost::serialization::load_construct_data(Archive & ar, FMT_Species * t, const unsigned int file_version)
{
    // retrieve data from archive required to construct new instance
  Density *d;
  ar >> d;

    // invoke inplace constructor to initialize instance of my_class
  double mu = 0;
  int seq_num = 0;
  ::new(t)FMT_Species(*d,mu,seq_num);

  ar >> t->mu_;
  ar >> t->seq_num_;  
  ar >> t->dF_;
  ar >> t->fixedMass_;
  ar >> t->SequenceNumber_;
  ar >> t->index_;  
  ar >> t->hsd_;
  ar >> t->fmt_weighted_densities;  
}

#endif // __LUTSKO__SPECIES__
