#ifndef __LUTSKO__SPECIES__
#define __LUTSKO__SPECIES__

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/export.hpp>

#include "FMT_Weighted_Density.h"
#include "Potential1.h"
#include "Density.h"

class Species
{
 public:
  Species(Density &density, double mu = 0, int seq = -1) : density_(density), dF_(density.Ntot()), mu_(mu), fixedMass_(-1) { if(seq >=0) {seq_num_ = seq; SequenceNumber_ = seq+1;} else seq_num_ = SequenceNumber_++;}
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
  void set_density(DFT_Vec &x) {density_.set(x);}
  void set_density(long j, double x) {density_.set_Density_Elem(j,x);}

  void zeroForce() {dF_.zeros();}
  void addToForce(long i, double v) {dF_.IncrementBy(i,v);}
  void addToForce(const DFT_Vec &f) {dF_.IncrementBy(f);}
  void setForce(const DFT_Vec &f) {dF_.set(f);}
  void multForce(double f) {dF_.MultBy(f);}  
  
  double get_convergence_monitor() const { return dF_.inf_norm()/density_.dV();}
  
  DFT_Vec &getDF() {return dF_;}
  void setDF(DFT_Vec &df) {return dF_.set(df);}

  // This species is held as allSpecies_[index_]
  void setIndex(int i) { index_ = i;}
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
    double Fx = density_.getExternalFieldEnergy()*dV - density_.getNumberAtoms()*mu_;
    if(bCalcForces)
      {
	dF_.IncrementBy_Scaled_Vector(density_.getExternalField(),dV);
	dF_.ShiftBy(mu_*dV);
      }
    return Fx;
  }
  /**
  *   @brief  Get the hard sphere diameter. This is a place-holder for the child objects that actually have a hsd.
  *  
  */  
  virtual double getHSD() const { return 0.0;}


  /**
   *   @brief  Constant particle number is enforced at the species-level. If needed, some information has to be collected before updating the forces. Note that particle number is rigorously kept constant.
   *  
   */    
  void beginForceCalculation()
  {
    if(fixedMass_ > 0.0)
      {
	density_.scale_to(fixedMass_/density_.getNumberAtoms());
	mu_ = 0.0;
      }
  }

  /**
   *   @brief  Constant particle number is enforced at the species-level. If activated, the necessary corrections to the forces are applied here. Note that particle number is rigorously kept constant.
   *  
   */    
  void endForceCalculation()
  {
    if(fixedMass_ > 0.0)
      {
	mu_ = 0.0;
      
	for(long p=0;p<density_.Ntot();p++)
	  mu_ += dF_.get(p)*density_.getDensity(p);
	mu_ /= fixedMass_;
	  
	for(long p=0;p<density_.Ntot();p++)
	  dF_.set(p, dF_.get(p)-mu_*density_.dV());
      }
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
  double fixedMass_; // if this is > 0, then the mass of this species is being held fixed at this value.
  int index_ = -1;
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
   *   @param  pointsFile contains the points for spherical integration: if it does not exist, the analytic weights are used
   *   @return nothing 
   */    
  FMT_Species(Density& density, double hsd, double mu = 0, int seq = -1);

  FMT_Species(const FMT_Species &) = delete;
  
  ~FMT_Species(){}

  //  void reset(string& pointsFile);
  
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
  double getEta(long pos) const { return d_[EI()].real(pos);}

  /**
   *   @brief  get value of S at pos
   *  
   *   @param  pos is the mesh position
   *   @return value of S at pos
   */      
  double getS(long pos) const { return d_[SI()].real(pos);}

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
  void Set_dPhi_Eta(long pos,double val)  {d_[EI()].Set_dPhi(pos,val);}

  /**
   *   @brief  set value of dPhi_dS at pos
   *  
   *   @param  pos is the mesh position
   *   @return none
   */        
  void Set_dPhi_S(long pos,double val) {d_[SI()].Set_dPhi(pos,val);}

  /**
   *   @brief  set value of dPhi_dV_j at pos
   *  
   *   @param  j is index of V
   *   @param  pos is the mesh position
   *   @return none
   */        
  void Set_dPhi_V(int j, long pos,double val) {d_[VI(j)].Set_dPhi(pos,val);}

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

  void convolute_eta_weight_with(const DFT_FFT &v, DFT_FFT &result, bool bConjugate = false) const
  { convolute_weight_with(EI(),v,result,bConjugate);}

  void convolute_s_weight_with(const DFT_FFT &v, DFT_FFT &result, bool bConjugate = false) const
  { convolute_weight_with(SI(),v,result,bConjugate);}

  void convolute_v_weight_with(int i, const DFT_FFT &v, DFT_FFT &result, bool bConjugate = false) const
  { convolute_weight_with(VI(i),v,result, bConjugate);}    

  void convolute_weight_with(int pos, const DFT_FFT &v, DFT_FFT &result, bool bConjugate = false) const
  {
    result.Four().Schur(v.cFour(), d_[pos].wk(),bConjugate);
    result.do_fourier_2_real();
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

  // These return the real weight at position K using the extended notation: eta, s0,s1,s2,v1,v2
  double getExtendedWeight(long K, int a)
  {
    if(a == 0) return d_[EI()].getWeight(K);

    if(a == 1) return d_[SI()].getWeight(K)/(hsd_*hsd_);
    if(a == 2) return d_[SI()].getWeight(K)/hsd_;
    if(a == 3) return d_[SI()].getWeight(K);

    if(a == 4) return d_[VI(0)].getWeight(K)/hsd_;
    if(a == 5) return d_[VI(1)].getWeight(K)/hsd_;
    if(a == 6) return d_[VI(2)].getWeight(K)/hsd_;

    if(a == 7) return d_[VI(0)].getWeight(K);
    if(a == 8) return d_[VI(1)].getWeight(K);
    if(a == 9) return d_[VI(2)].getWeight(K);

    throw std::runtime_error("Unknown index in FMT_Weighted_Density::getExtendedWeight");
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
	    d += getExtendedWeight(KI,a)*density_.getDensity(ix,iy,iz);
	  }
    return d;
  }
  // These return the weighted density at position K using the extended notation: eta, s0,s1,s2,v1,v2
  double getExtendedWeightedDensity(long K, int a)
  {
    if(a == 0) return d_[EI()].getDensity(K);

    if(a == 1) return d_[SI()].getDensity(K)/(hsd_*hsd_);
    if(a == 2) return d_[SI()].getDensity(K)/hsd_;
    if(a == 3) return d_[SI()].getDensity(K);

    if(a == 4) return d_[VI(0)].getDensity(K)/hsd_;
    if(a == 5) return d_[VI(1)].getDensity(K)/hsd_;
    if(a == 6) return d_[VI(2)].getDensity(K)/hsd_;

    if(a == 7) return d_[VI(0)].getDensity(K);
    if(a == 8) return d_[VI(1)].getDensity(K);
    if(a == 9) return d_[VI(2)].getDensity(K);

    throw std::runtime_error("Unknown index in FMT_Weighted_Density::getExtendedWeightedDensity");
  }

    
  FMT_Weighted_Density& getEta() { return d_[0];}
  
  // Used in DFT_Surfactant ...
  const DFT_Vec &getV_Real(int J) const { return d_[VI(J)].Real();}
  const DFT_Vec_Complex& getVweight_Four(int J) const { return d_[VI(J)].wk();}  

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive &ar, const unsigned int file_version)
  {
    boost::serialization::void_cast_register<FMT_Species, Species>(static_cast<FMT_Species *>(NULL),static_cast<Species *>(NULL));
  }    
  template<class Archive> friend void boost::serialization::save_construct_data(Archive & ar, const FMT_Species * t, const unsigned int file_version);
  template<class Archive> friend void boost::serialization::load_construct_data(Archive & ar, FMT_Species * t, const unsigned int file_version);

  
protected:
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

/**
   *  @brief FMT Species Class with weights generated numerically
   *
   *       Same as FMT_Species with weights generated from tabulated integration points on a sphere
   */

class FMT_Species_Numeric : public FMT_Species
{
public:
  /**
   *   @brief  Default  constructor for FMT_Species 
   *  
   *   @param  hsd is the hard-sphere diameter
   *   @param  lattice describes the mesh
   *   @param  pointsFile contains the points for spherical integration: if it does not exist, the analytic weights are used
   *   @return nothing 
   */    
  FMT_Species_Numeric(Density& density, double hsd, string &pointsFile, double mu = 0, int seq = -1);

  FMT_Species_Numeric(const FMT_Species &) = delete;
  
  ~FMT_Species_Numeric(){}

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive &ar, const unsigned int file_version)
  {
        boost::serialization::void_cast_register<FMT_Species_Numeric, FMT_Species>(static_cast<FMT_Species_Numeric *>(NULL),static_cast<FMT_Species *>(NULL));
  }     
  template<class Archive> friend void boost::serialization::save_construct_data(Archive & ar, const FMT_Species_Numeric * t, const unsigned int file_version);
  template<class Archive> friend void boost::serialization::load_construct_data(Archive & ar, FMT_Species_Numeric * t, const unsigned int file_version);
  
 protected:
  /**
   *   @brief  This is a one-time-only evaluation of the numerical approximation to the FMT weight functions. These are all 
   *           functions w_{alpha}(i,j) = w_{alpha}(abs(i-j)). Their evaluation involves real-space integrations for which the 
   *           integration points are given in the file pointsFile. Most of the work occurs via a call to initializeWeightedDensities@
   *
   *   @param  pointsFile is the file holding the integration points for integrating a spherical shell
   */        
  void generateWeights(string &pointsFile);



};



/**
   *  @brief An FMT Species Class with analytic weight generation 
   *
   *    Same as FMT_Species except that the weights are generated from analytic formulas. Tensors are not yet implemented.
   */

class FMT_Species_Analytic : public FMT_Species
{
public:
  /**
   *   @brief  Default  constructor for FMT_Species 
   *  
   *   @param  hsd is the hard-sphere diameter
   *   @param  lattice describes the mesh
   *   @return nothing 
   */    
  FMT_Species_Analytic(Density& density, double hsd, double mu = 0, int seq = -1);

  FMT_Species_Analytic(const FMT_Species &) = delete;
  
  ~FMT_Species_Analytic(){}

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive &ar, const unsigned int file_version)
  {
    boost::serialization::void_cast_register<FMT_Species_Analytic, FMT_Species>(static_cast<FMT_Species_Analytic *>(NULL),static_cast<FMT_Species *>(NULL));
  }      
  template<class Archive> friend void boost::serialization::save_construct_data(Archive & ar, const FMT_Species_Analytic * t, const unsigned int file_version);
  template<class Archive> friend void boost::serialization::load_construct_data(Archive & ar, FMT_Species_Analytic * t, const unsigned int file_version);


  
protected:
  /**
   *   @brief  This is a one-time-only evaluation of the numerical approximation to the FMT weight functions. These are all 
   *           functions w_{alpha}(i,j) = w_{alpha}(abs(i-j)). 
   *
   */        
  void generateWeights();
};

/**
   *  @brief An FMT Species Class with analytic weight generation 
   *
   *    Same as FMT_Species except that the weights are generated from analytic formulas. Tensors are not yet implemented.
   */

class FMT_Species_Analytic_2 : public FMT_Species
{
public:
  /**
   *   @brief  Default  constructor for FMT_Species 
   *  
   *   @param  hsd is the hard-sphere diameter
   *   @param  lattice describes the mesh
   *   @return nothing 
   */    
  FMT_Species_Analytic_2(Density& density, double hsd, double mu = 0, int seq = -1);

  FMT_Species_Analytic_2(const FMT_Species &) = delete;
  
  ~FMT_Species_Analytic_2(){}

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive &ar, const unsigned int file_version)
  {
    boost::serialization::void_cast_register<FMT_Species_Analytic_2, FMT_Species>(static_cast<FMT_Species_Analytic *>(NULL),static_cast<FMT_Species *>(NULL));
  }      
  template<class Archive> friend void boost::serialization::save_construct_data(Archive & ar, const FMT_Species_Analytic_2 * t, const unsigned int file_version);
  template<class Archive> friend void boost::serialization::load_construct_data(Archive & ar, FMT_Species_Analytic_2 * t, const unsigned int file_version);


  
protected:
  /**
   *   @brief  This is a one-time-only evaluation of the numerical approximation to the FMT weight functions. These are all 
   *           functions w_{alpha}(i,j) = w_{alpha}(abs(i-j)). 
   *
   */        
  void generateWeights();
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
  ar << t->d_;
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
  ar >> t->d_;  
}

template<class Archive>
inline void boost::serialization::save_construct_data(Archive & ar, const FMT_Species_Analytic * t, const unsigned int file_version)
{
  //  ar << static_cast<const FMT_Species*>(t);
  ar << & t->density_;
  ar << t->mu_;
  ar << t->seq_num_;
  ar << t->dF_;
  ar << t->fixedMass_;
  ar << t->SequenceNumber_;
  ar << t->index_;
  
  ar << t->hsd_;
  ar << t->d_;  
}

template<class Archive>
inline void boost::serialization::load_construct_data(Archive & ar, FMT_Species_Analytic * t, const unsigned int file_version)
{
    // retrieve data from archive required to construct new instance
  Density *d;
  ar >> d;

    // invoke inplace constructor to initialize instance of my_class
  double mu = 0;
  int seq_num = 0;
  double hsd = 1; // place holder
  ::new(t)FMT_Species_Analytic(*d,hsd, mu, seq_num);
  
  ar >> t->mu_;
  ar >> t->seq_num_;  
  ar >> t->dF_;
  ar >> t->fixedMass_;
  ar >> t->SequenceNumber_;
  ar >> t->index_;
  
  ar >> t->hsd_;
  ar >> t->d_;  
}

template<class Archive>
inline void boost::serialization::save_construct_data(Archive & ar, const FMT_Species_Analytic_2 * t, const unsigned int file_version)
{
  //  ar << static_cast<const FMT_Species*>(t);
  ar << & t->density_;
  ar << t->mu_;
  ar << t->seq_num_;
  ar << t->dF_;
  ar << t->fixedMass_;
  ar << t->SequenceNumber_;
  ar << t->index_;
  
  ar << t->hsd_;
  ar << t->d_;  
}

template<class Archive>
inline void boost::serialization::load_construct_data(Archive & ar, FMT_Species_Analytic_2 * t, const unsigned int file_version)
{
    // retrieve data from archive required to construct new instance
  Density *d;
  ar >> d;

    // invoke inplace constructor to initialize instance of my_class
  double mu = 0;
  int seq_num = 0;
  double hsd = 1; // place holder
  ::new(t)FMT_Species_Analytic_2(*d,hsd, mu, seq_num);
  
  ar >> t->mu_;
  ar >> t->seq_num_;  
  ar >> t->dF_;
  ar >> t->fixedMass_;
  ar >> t->SequenceNumber_;
  ar >> t->index_;
  
  ar >> t->hsd_;
  ar >> t->d_;  
}

template<class Archive>
inline void boost::serialization::save_construct_data(Archive & ar, const FMT_Species_Numeric * t, const unsigned int file_version)
{
  //  ar << static_cast<const FMT_Species*>(t);
  ar << & t->density_;
  ar << t->mu_;
  ar << t->seq_num_;
  ar << t->dF_;
  ar << t->fixedMass_;
  ar << t->SequenceNumber_;
  ar << t->index_;
  
  ar << t->hsd_;
  ar << t->d_;  
}

template<class Archive>
inline void boost::serialization::load_construct_data(Archive & ar, FMT_Species_Numeric * t, const unsigned int file_version)
{
    // retrieve data from archive required to construct new instance
  Density *d;
  ar >> d;

  // invoke inplace constructor to initialize instance of my_class
  double mu = 0;
  double hsd = 1;
  int seq_num = 0;
  string noFile;
  ::new(t)FMT_Species_Numeric(*d, hsd, noFile, 0, seq_num);

  ar >> t->mu_;
  ar >> t->seq_num_;    
  ar >> t->dF_;
  ar >> t->fixedMass_;
  ar >> t->SequenceNumber_;
  ar >> t->index_;
  
  ar >> t->hsd_;
  ar >> t->d_;  
}

#endif // __LUTSKO__SPECIES__
