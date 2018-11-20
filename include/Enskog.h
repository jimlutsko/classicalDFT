#ifndef __ENSKOG__H
#define __ENSKOG__H


/** This class encapsulates the properties of the hard-sphere fluid.
 */

/**
  *  @brief  Thermodynamics and other properties of hard-sphere fluids
  *
  *  @detailed This class enapsulates various analytic properties and models for hard-sphere fluids (PY and Carnahan-Starling thermodyynamics, viscosities in Enskog approximation, ... Note that the hard-sphere diameter is taken to be unity.
  */  

class Enskog
{
 public:
   /**
    *   @brief  Constructor
    *
    *   @param  denisty is the density reduced by the hard-sphere diameter $\\rho d^3$
    */     
  Enskog(double density) {setDensity(density);}

  /**
    *   @brief  Change the density
    *
    *   @param  denisty is the density reduced by the hard-sphere diameter $\\rho d^3$
    */    
  void setDensity(double density){ n = density; e = M_PI*density/6.0; c = (1-e/2)/(1-3*e+3*e*e-e*e*e);}

  /**
   *   @brief  Carnahan-Starling chi (pdf at contact)
   *
   *   @return chi = (1-e/2)/(1-3*e+3*e*e-e*e*e).
   */
  
  double getChiCS(){ return c;}

  /**
   *   @brief  PY excess free energy (virial route) per particle
   *
   *   @return excess free energy per particle for d = kT = 1
   */      
  double exFreeEnergyPYV(){return 2*log(1-e)+6*e/(1-e);}
  /**
   *   @brief  PY excess free energy (compressibility) per particle
   *
   *   @return excess free energy per particle for d = kT = 1
   */        
  double exFreeEnergyPYC(){return -log(1-e)+1.5*e*(2-e)/(1-2*e+e*e);}

  /**
   *   @brief  CS excess free energy per particle
   *
   *   @return excess free energy per particle for d = kT = 1
   */        
  double exFreeEnergyCS(){return e*(4.0-3.0*e)/(1-2*e+e*e);}

  /**
   *   @brief first density derivative d/dn for PY excess free energy (virial route) per particle
   *
   *   @return first density derivative d/dn for excess free energy per particle for d = kT = 1
   */        
  double dexFreeEnergyPYVdRho(){return (M_PI/6.0)*(4.0+2.0*e)/(1-2*e+e*e);}
    /**
   *   @brief  first density derivative d/dn for PY excess free energy (compressibility) per particle
   *
   *   @return first density derivative d/dn for excess free energy per particle for d = kT = 1
   */        
  double dexFreeEnergyPYCdRho(){return (M_PI/6.0)*(4.0-2.0*e+e*e)/(1-3*e+3*e*e-e*e*e);}
  /**
   *   @brief  first density derivative d/dn for CS excess free energy per particle
   *
   *   @return first density derivative d/dn for excess free energy per particle for d = kT = 1
   */          
  double dexFreeEnergyCSdRho(){return (M_PI/6.0)*(4.0-2.0*e)/(1-3*e+3*e*e-e*e*e);}

  /**
   *   @brief  second density derivative d2/dn2 for PY excess free energy (virial route) per particle
   *
   *   @return second density derivative d2/dn2 for excess free energy per particle for d = kT = 1
   */          
  double d2exFreeEnergyPYVdRho2(){return (M_PI/6.0)*(M_PI/6.0)*(10.0+2.0*e)/(1-3*e+3*e*e-e*e*e);}
  /**
   *   @brief  second density derivative d2/dn2 for PY excess free energy (compressibility) per particle
   *
   *   @return second density derivative d2/dn2 for excess free energy per particle for d = kT = 1
   */          
  double d2exFreeEnergyPYCdRho2(){return (M_PI/6.0)*(M_PI/6.0)*(10.0-2.0*e+e*e)/(1-4*e+6*e*e-4*e*e*e+e*e*e*e);}
  /**
   *   @brief  second density derivative d2/dn2 of CS excess free energy per particle
   *
   *   @return second density derivative d2/dn2 of excess free energy per particle for d = kT = 1
   */          
  double d2exFreeEnergyCSdRho2(){return (M_PI/6.0)*(M_PI/6.0)*(10.0-4.0*e)/(1-4*e+6*e*e-4*e*e*e+e*e*e*e);}

  /**
   *   @brief  third density derivative d3/dn3 for PY excess free energy (compressibility) per particle
   *
   *   @return third density derivative d3/dn3 for excess free energy per particle for d = kT = 1
   */            
  double d3exFreeEnergyPYCdRho3(){return (M_PI/6.0)*(M_PI/6.0)*(M_PI/6.0)*(38-4*e+2*e*e)*pow(1-e,-5);}

  /**
   *   @brief  third density derivative d3/dn3 of CS excess free energy per particle
   *
   *   @return third density derivative d3/dn3 of excess free energy per particle for d = kT = 1
   */            
  double d3exFreeEnergyCSdRho3(){return (M_PI/6.0)*(M_PI/6.0)*(M_PI/6.0)*12*(3-e)*pow(1-e,-5);}

  /**
   *   @brief  PY total free energy (virial) per particle
   *
   *   @return total free energy per particle for d = kT = 1
   */          
  double freeEnergyPYV(){return log(n)-1+exFreeEnergyPYV();}

  /**
   *   @brief  PY total free energy (compressibility) per particle
   *
   *   @return total free energy per particle for d = kT = 1
   */            
  double freeEnergyPYC(){return log(n)-1+exFreeEnergyPYC();}

  /**
   *   @brief  CS total free energy per particle
   *
   *   @return total free energy per particle for d = kT = 1
   */            
  double freeEnergyCS(){return log(n)-1+exFreeEnergyCS();}

  /**
   *   @brief  PY pressure (virial) relative to ideal gas
   *
   *   @return P/nKT
   */            
  double pressurePYV(){return (1+2*e+3*e*e)/(1-2*e+e*e);}

  /**
   *   @brief  PY pressure (compressibility) relative to ideal gas
   *
   *   @return P/nKT
   */              
  double pressurePYC(){return (1+e+e*e)/(1-3*e+3*e*e-e*e*e);}

  /**
   *   @brief  CS pressure relative to ideal gas
   *
   *   @return P/nKT
   */              
  double pressureCS(){return (1+e+e*e-e*e*e)/(1-3*e+3*e*e-e*e*e);}

  /**
   *   @brief  Determines density from pressure (CS eos)
   *
   *   @return density times d^3
   */              
  double Density_From_Pressure_CS(double p) const;

  /**
   *   @brief  virial PY chemical potential
   *
   *   @return mu/kT
   */            
  double chemPotentialPYV(){ return freeEnergyPYV()+pressurePYV();}

  /**
   *   @brief  compressibility PY chemical potential
   *
   *   @return mu/kT
   */              
  double chemPotentialPYC(){ return freeEnergyPYC()+pressurePYC();}

  /**
   *   @brief  CSchemical potential
   *
   *   @return mu/kT
   */              
  double chemPotentialCS(){ return freeEnergyCS()+pressureCS();}

  /**
   *   @brief  Heat capacity at constant volume
   *
   *   @return dimensionless heat capacity
   */            
  double cv() { return 3.0/2.0;}

  /**
   *   @brief  Heat capacity at constant pressure
   *
   *   @return dimensionless heat capacity
   */            
  double cp() { return 5.0/2.0;}

  /**
   *   @brief  ratio of heat capacities cp/cv
   *
   *   @return ratio cp/cv
   */              
  double gamma() { return cp()/cv();}

  /**
   *   @brief  CS pressure divided by nkT P/nkT
   *
   *   @return CS pressure divided by kT P/nkT
   */            
  double pressure() { return 1+4*e*c;}


  /**
    *   @brief  Enskog expression for shear viscosity (reduced by temperature)
    *
    *   @param  
    *   @return shear viscosity with $d$ and $kT = 1$.
    */    
  double shearViscosity()
    {
      return (5.0/(16.0*sqrt(M_PI)*c))*((1+4*M_PI*n*c/15.0)*(1+4*M_PI*n*c/15.0))+0.6*bulkViscosity();
    }

  /**
    *   @brief  Enskog expression for bulk viscosity (reduced by temperature)
    *
    *   @param  
    *   @return bulk viscosity with $d$ and $kT = 1$.
    */      
  double bulkViscosity()
    {
      return (5.0/(16.0*sqrt(M_PI)))*((64.0/45.0)*M_PI*n*n*c);
    }

    /**
    *   @brief  Enskog expression for thermal conductivity(reduced by temperature)
    *
    *   @param  
    *   @return thermal conductivity with $d$ and $kT = 1$.
    */      
  double thermalConductivity()
    {
      return (75.0/(64.0*sqrt(M_PI)*c))*(1+0.4*M_PI*n*c)*(1+0.4*M_PI*n*c)+(128.0/225.0)*M_PI*n*n*c;
    }

    /**
    *   @brief  Enskog expression for sound dampint constant (reduced by temperature)
    *
    *   @param  
    *   @return sound damping constant for $d$ and $kT = 1$.
    */      
  double soundDampingConstant()
    {
      return (((gamma()-1)/cp())*thermalConductivity()+(4.0/3.0)*shearViscosity()+bulkViscosity())/n;
    }



  
 protected:
  double n,e,c;
};



#endif
