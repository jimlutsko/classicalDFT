#ifndef __LUTSKO__EOS__
#define __LUTSKO__EOS__

class EOS
{
 public:
  EOS(){}
  ~EOS(){}

  virtual double f(double density, double kT) {return 0.0;}
  virtual double dfdx(double density, double kT) {return 0.0;}
  virtual double d2fd2xf(double density, double kT) {return 0.0;}

};



#endif // __LUTSKO__EOS__
