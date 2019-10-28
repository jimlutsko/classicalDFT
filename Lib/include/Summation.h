/* This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. 
 * To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
 *
 * Author: James F. Lutsko
 * www.lutsko.com
 */

#ifndef __LUTSKO__SUMMATION__
#define __LUTSKO__SUMMATION__
#include <quadmath.h>

/**
  *  @brief Implementation of the improved Kahan–Babuška algorithm of Neumaier (Zeitschrift für Angewandte Mathematik und Mechanik (in German). 54 (1): 39–51).
  */  
class Summation
{
 public:
 Summation() : sum_(0.0), c_(0.0) {}
  ~Summation(){}

  void add(double val)
  {
    double t = sum_ + val;
    if(fabs(sum_) >= fabs(val))
      c_ += (sum_-t) + val;
    else c_ += (val-t) + sum_;
    sum_ = t;
  }
  double sum() const { return double(sum_ + c_);}

  Summation& operator+=(const double rhs)
    {
      this->add(rhs);
      return *this;
    }

  operator double() const { return sum(); }


  Summation& operator += (const Summation & other)
    {
      sum_ += other.sum_;
      c_ += other.c_;
      return *this;
    }

  Summation& operator=(double val) 
    {
      sum_ = val;
      c_= 0.0; // what other option?
      return *this;
    }
  
 protected:
  double  sum_;
  double c_;
};











#endif  // __LUTSKO__SUMMATION__
