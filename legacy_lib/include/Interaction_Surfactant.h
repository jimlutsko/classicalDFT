#ifndef __classicalDFT_INTERATION_SURFACTANT__
#define  __classicalDFT_INTERATION_SURFACTANT__

#include "Potential1.h"
#include "Interaction.h"

class Interaction_Surfactant : public Interaction_Interpolation_QF
{
 public:
 Interaction_Surfactant(FMT_Species *water, Species *surfactant, Potential1 *v, double kT)
   :   Interaction_Interpolation_QF(water,surfactant,v,kT) {}

  // Over-ride bulk contributions since these are zero in the absence of gradients
  virtual double Mu(const vector<double> &x, int species) {return 0.0;}
  virtual double Fhelmholtz(const vector<double> &x) const {return 0.0;}


  /**
   *   @brief  This executes the basic functionality of computing energies and forces
   *  
   *   @returns The mean-field contribution to the (total) free energy divided by kT.
   */    
  virtual double getInteractionEnergyAndForces()
  {
    // Species 1 is the **water** and species 2 is the surfactant.

    FMT_Species &water = *((FMT_Species*) s1_);
    Species &surfactant = *((FMT_Species*) s2_);
  
    const Density& surf_density = surfactant.getDensity();

    int Nx = surf_density.Nx();
    int Ny = surf_density.Ny();
    int Nz = surf_density.Nz();

    long Ntot = Nx*Ny*Nz;
  
    double dV = surf_density.dV();

    ////////////////////////////////////////////////////////////////////////////
    // Construct v2(r) for the water & do fft

    const DFT_Vec& vRealx = water.getV_Real(0);
    const DFT_Vec& vRealy = water.getV_Real(1);
    const DFT_Vec& vRealz = water.getV_Real(2);

    DFT_FFT v2_water(Nx, Ny, Nz);      
    for(long i=0;i<Ntot;i++)
      v2_water.Real().set(i, vRealx.get(i)*vRealx.get(i)+vRealy.get(i)*vRealy.get(i)+vRealz.get(i)*vRealz.get(i));
    v2_water.do_real_2_fourier();

    //////////////////////////////////////////////////////////////////////////
    // Surfactant force is just convoluton of the potential and v2_water
    // Lets do this the dumb way first:
    DFT_FFT dF(Nx, Ny, Nz);      

    dF.Four().Schur(w_att_.Four(),v2_water.Four(),false); //true);
    dF.Four().MultBy(1.0/dF.Real().size());  // FFTW's convention
    dF.do_fourier_2_real();
    dF.Real().MultBy(dV*dV); // this is d F/d rho_i and so the factors of dV survive

    surfactant.addToForce(dF.Real());
  
    /////////////////////////////////////////////////////////////////////////////////
    // Contribution to the free energy is this dotted with the surfactant density

    double F = surf_density.getInteractionEnergy(dF.Real());
    
    ////////////////////////////////////////////////////////////////////////
    // Now we need the derivative wrt the water density

    for(int i=0;i<3;i++)
      { 
	dF.zeros();

	dF.Four().Schur(surf_density.getDK(), w_att_.Four(),false);
	dF.Four().MultBy(1.0/dF.Real().size());  // FFTW's convention
	dF.do_fourier_2_real();  

	dF.Real().Schur(dF.Real(),water.getV_Real(i));

	dF.do_real_2_fourier();
	dF.Four().Schur(dF.Four(),water.getVweight_Four(i),false); // this includes the FFTW scale factor already
	dF.do_fourier_2_real();

	dF.Real().MultBy(2*dV*dV);
      
	water.addToForce(dF.cReal());
      }

    return F;
  }

 protected:
  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Interaction_Interpolation_QF>(*this);
    boost::serialization::void_cast_register<Interaction_Surfactant, Interaction_Interpolation_QF>(static_cast<Interaction_Surfactant *>(NULL),static_cast<Interaction_Interpolation_QF *>(NULL));    
  }    
};





#endif  // __classicalDFT_INTERATION_SURFACTANT__
