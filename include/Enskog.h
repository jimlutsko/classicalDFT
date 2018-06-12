#ifndef __ENSKOG__H
#define __ENSKOG__H


/** This class encapsulates the properties of the hard-sphere fluid.
 */
class Enskog
{
 public:
  Enskog(double density) {setDensity(density);}

  void setDensity(double density){ n = density; e = M_PI*density/6.0; c = (1-e/2)/(1-3*e+3*e*e-e*e*e);}

  double shearViscosity()
    {
      return (5.0/(16.0*sqrt(M_PI)*c))*((1+4*M_PI*n*c/15.0)*(1+4*M_PI*n*c/15.0))+0.6*bulkViscosity();
    }

  double bulkViscosity()
    {
      return (5.0/(16.0*sqrt(M_PI)))*((64.0/45.0)*M_PI*n*n*c);
    }
      
  double thermalConductivity()
    {
      return (75.0/(64.0*sqrt(M_PI)*c))*(1+0.4*M_PI*n*c)*(1+0.4*M_PI*n*c)+(128.0/225.0)*M_PI*n*n*c;
    }

  double soundDampingConstant()
    {
      return (((gamma()-1)/cp())*thermalConductivity()+(4.0/3.0)*shearViscosity()+bulkViscosity())/n;
    }

  double getChiCS(){ return c;}

  double exFreeEnergyPYV(){return 2*log(1-e)+6*e/(1-e);}
  double exFreeEnergyPYC(){return -log(1-e)+1.5*e*(2-e)/(1-2*e+e*e);}
  double exFreeEnergyCS(){return e*(4.0-3.0*e)/(1-2*e+e*e);}

  double dexFreeEnergyPYVdRho(){return (M_PI/6.0)*(4.0+2.0*e)/(1-2*e+e*e);}
  double dexFreeEnergyPYCdRho(){return (M_PI/6.0)*(4.0-2.0*e+e*e)/(1-3*e+3*e*e-e*e*e);}
  double dexFreeEnergyCSdRho(){return (M_PI/6.0)*(4.0-2.0*e)/(1-3*e+3*e*e-e*e*e);}

  double d2exFreeEnergyPYVdRho2(){return (M_PI/6.0)*(M_PI/6.0)*(10.0+2.0*e)/(1-3*e+3*e*e-e*e*e);}
  double d2exFreeEnergyPYCdRho2(){return (M_PI/6.0)*(M_PI/6.0)*(10.0-2.0*e+e*e)/(1-4*e+6*e*e-4*e*e*e+e*e*e*e);}
  double d2exFreeEnergyCSdRho2(){return (M_PI/6.0)*(M_PI/6.0)*(10.0-4.0*e)/(1-4*e+6*e*e-4*e*e*e+e*e*e*e);}

  double d3exFreeEnergyPYCdRho3(){return (M_PI/6.0)*(M_PI/6.0)*(M_PI/6.0)*(38-4*e+2*e*e)*pow(1-e,-5);}
  double d3exFreeEnergyCSdRho3(){return (M_PI/6.0)*(M_PI/6.0)*(M_PI/6.0)*12*(3-e)*pow(1-e,-5);}

  double freeEnergyPYV(){return log(n)-1+exFreeEnergyPYV();}
  double freeEnergyPYC(){return log(n)-1+exFreeEnergyPYC();}
  double freeEnergyCS(){return log(n)-1+exFreeEnergyCS();}

  double pressurePYV(){return (1+2*e+3*e*e)/(1-2*e+e*e);}
  double pressurePYC(){return (1+e+e*e)/(1-3*e+3*e*e-e*e*e);}
  double pressureCS(){return (1+e+e*e-e*e*e)/(1-3*e+3*e*e-e*e*e);}

  double Density_From_Pressure_CS(double p) const;

  double chemPotentialPYV(){ return freeEnergyPYV()+pressurePYV();}
  double chemPotentialPYC(){ return freeEnergyPYC()+pressurePYC();}
  double chemPotentialCS(){ return freeEnergyCS()+pressureCS();}

  double cv() { return 3.0/2.0;}

  double cp() { return 5.0/2.0;}

  double gamma() { return cp()/cv();}
  // P/kT
  double pressure() { return 1+4*e*c;}

  // protected:
  double n,e,c;
};



#endif
