#ifndef __LUTSKO__EOS__
#define __LUTSKO__EOS__

class EOS
{
 public:
 EOS(double kT):  kT_(kT) {}
  ~EOS(){}

  virtual double fex(double density) {return 0.0;}
  virtual double dfexdx(double density) {return 0.0;}
  virtual double d2fexdx2(double density) {return 0.0;}

 protected:
  double kT_ = 1;
  
};



#endif // __LUTSKO__EOS__
