#ifndef __LUTSKO_STRING_METHOD__
#define __LUTSKO_STRING_METHOD__

#include "Minimizer.h"
#include "Grace.h"

class StringMethod_MPI
{
 public:
 StringMethod_MPI(bool freeEnd) : freeEnd_(freeEnd), step_counter_(0){};
  ~StringMethod_MPI(){};

  void setMu(double m) { mu_ = m;}
  virtual void run(string& logfile) = 0;

  
 protected:
  double mu_;
  bool freeEnd_;
  long step_counter_;
};


class StringMethod_MPI_Master : public StringMethod_MPI
{
 public:
 StringMethod_MPI_Master(int Nimages, long Ntot, Grace *g = NULL,bool freeEnd = false) 
   : StringMethod_MPI(freeEnd), grace_(g)
    {
      Images_.resize(Nimages);
    }

  virtual void run(string& logfile);

  void Display(vector<double> &F, int dataSet, double dFmax = 0.0, double dFav = 0.0);


  void addTask(int images) { taskList.push_back(images);}
  
 private:
    Grace *grace_;
    vector<int> taskList;

    vector< vector<double> > Images_;
    long Ntot_;

};


class StringMethod_MPI_Slave : public StringMethod_MPI
{
 public:
 StringMethod_MPI_Slave(DDFT &ddft, vector<Density*> string, bool freeEnd = false) 
   : StringMethod_MPI(freeEnd), ddft_(ddft), string_(string), Jclimb_(-1)
  {
    int N = string_.size();  

    gr_ = new mglGraph;
    DT_.resize(N,0.0);
    oldF_.resize(N,0.0);
    N_.resize(N,0.0);

    // Initialize free energies  
    for(int J=0;J<N;J++)
      {
	oldF_[J] = ddft_.F_string(*(string_[J])) - mu_*string_[J]->getNumberAtoms();
	DT_[J] = 0.01; 
      }

    // this holds a copy of the string
    string_copy_.resize(N);
    for(int J=0;J<N;J++)
      string_copy_[J].zeros(string_[0]->Ntot());
    
    //We just need one density to set up the array which is then shared by all
    string_[N-1]->initialize_2D_data(data_2D_);
  }
  
  ~StringMethod_MPI_Slave(){if(gr_) delete gr_;}

  virtual void run(string& logfile);

  void setClimbingImage(int J) { Jclimb_ = J; tangent_.zeros(string_[0]->Ntot());}

  void archive(string &filename) const;
  void read(string &filename);
  
 private:
  void rescale(double Ntarget);
  void rescale_linear();
  void Draw(Density& density, int image_number, double F);
  
 private:
  DDFT &ddft_;
  vector<Density*> string_;
  vector<DFT_Vec> string_copy_;   // this holds a copy of the string

  vector<double> DT_;
  vector<double> oldF_;
  vector<double> N_;
  DFT_Vec tangent_;
  bool freeEnd_;

  int Jclimb_;
  
  mglGraph *gr_;
  mglData data_2D_;
};



#endif //#ifndef __LUTSKO_STRING_METHOD__
