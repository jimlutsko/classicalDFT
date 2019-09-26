#ifndef __LUTSKO_STRING_METHOD__
#define __LUTSKO_STRING_METHOD__

#include "Minimizer.h"
#include "Grace.h"

/**
  *  @brief Nucleation Pathways via the String Method and DDFT dynamics. Practically obselete in favor of StringMethod_MPI.
  */  


class StringMethod
{
 public:
 StringMethod(DDFT &ddft, vector<Density*> string, double Time_Step_Max = 1e-2, Grace *g = NULL,bool freeEnd = false) : ddft_(ddft), string_(string), grace_(g), Time_Step_Max_(Time_Step_Max), freeEnd_(freeEnd), Jclimb_(-1)
  {
    DT_.resize(string_.size(),0.0);
  }
  
  ~StringMethod(){} 

  void setClimbingImage(int J) { Jclimb_ = J; tangent_.zeros(string_[0]->Ntot());}

  void setMu(double m) { mu_ = m;}
  
  void run(string& logfile);

  void archive(string &filename) const;
  void read(string &filename);
  
 private:
  void rescale(double Ntarget);
  void rescale_linear();
  void Draw(Density& density, int image_number, double F);
  void Display(vector<double> &F, int dataSet, double dFmax = 0.0, double dFav = 0.0);
  
 private:
  DDFT &ddft_;
  vector<Density*> string_;
  vector<double> DT_;
  DFT_Vec tangent_;
  long step_counter_;
  Grace *grace_;
  bool freeEnd_;

  int Jclimb_;
  
  double mu_;
  double Time_Step_Max_;
};



#endif //#ifndef __LUTSKO_STRING_METHOD__
