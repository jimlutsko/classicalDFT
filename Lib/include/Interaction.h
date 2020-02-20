/* This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
 * To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
 *
 * Author: James F. Lutsko
 * www.lutsko.com
 */

#ifndef __LUTSKO__INTERACTION__
#define __LUTSKO__INTERACTION__

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cstring>
#include <complex>
#include <complex.h>

#include "Species.h"
#include "Log.h"
#include "myColor.h"

/**
  *  @brief This class encapsulates the interaction between two species (or one species with itself)
  */  


class Interaction_Base
{
 public:

 Interaction_Base(Species &s1, Species &s2, Potential1 &v, double kT, Log &log) : s1_(s1), s2_(s2), v_(v), kT_(kT), log_(log), initialized_(false) {}

  void initialize();

  virtual bool checkWeightsFile(ifstream &in) = 0;
  virtual void generateWeights(Potential1 &v, stringstream &ss, Log& log) = 0;  
  // Note that the matrix w already contains a factor of dV
  double getInteractionEnergyAndForces();

  double checkCalc(int jx, int jy, int jz);

  double getW(long s)  { return w_att_.Real().get(s);}

  double Mu(const vector<double> &x, int species) const
  {
    double mu = 0.0;

    if(s1_.getSequenceNumber() == species)
      mu += 0.5*a_vdw_*x[s2_.getSequenceNumber()];

    if(s2_.getSequenceNumber() == species)
      mu += 0.5*a_vdw_*x[s1_.getSequenceNumber()];    

    return mu;
  }

  double Fhelmholtz(const vector<double> &x) const {return 0.5*a_vdw_*x[s1_.getSequenceNumber()]*x[s2_.getSequenceNumber()];}  

  double getVDWParameter() const { if(!initialized_) throw std::runtime_error("Interaction object must be initialized before calling getVDWParameter()"); return a_vdw_;}

 protected:
  virtual bool readWeights(stringstream &ss1);
  
 protected:
  Species &s1_;
  Species &s2_;

  double a_vdw_;
  DFT_FFT w_att_;

  bool initialized_ = false;
  Potential1 &v_;
  double kT_;
  Log& log_;
};



class Interaction : public Interaction_Base
{
 public:
 Interaction(Species &s1, Species &s2, Potential1 &v, double kT, Log &log, string& pointsFile) :
  Interaction_Base(s1,s2,v,kT,log), pointsFile_(pointsFile) {};

  virtual bool checkWeightsFile(ifstream &in);
  virtual void generateWeights(Potential1 &v, stringstream &ss, Log& log);
  
 protected:
  string pointsFile_;
};

class Interaction_Full : public Interaction_Base
{
 public:

 Interaction_Full(Species &s1, Species &s2, Potential1 &v, double kT, Log &log, int Ngauss) :
  Interaction_Base(s1,s2,v,kT,log), Ngauss_(Ngauss) {}

  // TODO:
  virtual bool checkWeightsFile(ifstream &in) {return true;}

  virtual void generateWeights(Potential1 &v, stringstream &ss, Log& log);
 protected:
  int Ngauss_;
};


class Interaction_Linear_Interpolation : public Interaction_Base
{
 public:

 Interaction_Linear_Interpolation(Species &s1, Species &s2, Potential1 &v, double kT, Log &log) :
  Interaction_Base(s1,s2,v,kT,log) {}  

  virtual void generateWeights(Potential1 &v, stringstream &ss, Log& log);
  virtual bool checkWeightsFile(ifstream &in) {return false;}
  
 protected:
  virtual bool readWeights(stringstream &ss1) { return false;} 
};





#endif // __LUTSKO__INTERACTION__
