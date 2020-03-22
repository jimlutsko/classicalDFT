#ifndef __LUTSKO_POLYNOMIAL__
#define __LUTSKO_POLYNOMIAL__


#include <cmath>
#include <vector>


/**
  *  @brief Polynomial Class. This is a helper class for manipulating polynomials. Its (initial) goal was to encode a reliable method of determining the real roots of a polynomial with real coefficients.It may not be efficient, but it is reliable. 
  *
  *  @detailed The method of finding roots is taken from   . 
  */  
class Polynomial
{
public:
  /**
   *   @brief  Constructor for Polynomial object
   *  
   *   @param c: vector containing the coefficients: p(x) = c[0] + c[1] x + c[2] x^2 ...
   */  
  Polynomial(vector<double> &c) : c_(c){}

  /**
   *   @brief  Size of polynomial: 1 + its order 
   *  
   *   @returns: size of array of coefficients defining the polynomial
   */  
  size_t size() const { return c_.size();}


  /**
   *   @brief  Evaluate polynomial at x
   *  
   *   @param: x is the input vaue
   *   @returns: p(x)
   */    
  double eval(double x)
  {
    double ret = 0;  
    for(int i=c_.size()-1;i>=0;i--)
      ret = c_[i]+x*ret;
    return ret;
  }

  /**
   *   @brief  Determine a root by bisection in the interval (a,b]
   *  
   *   @param: a is the lower bound of the interval
   *   @param: b is the upper bound of the interval
   *   @param: tol is the tolerance
   *   @returns: the root
   */    
  
  double bisect(double a, double b, double tol)
  {
    double pa = eval(a);
    double pb = eval(b);
    
    if(fabs(pb) < tol) return b;

    if(pa*pb < 0) throw std::runtime_error("Polynomial class asked to find root in invalid interval");
    
    double c,pc;  
    do{
      c = (a+b)/2;
      pc = eval(c);
      if(pc*pb > 0) {b = c; pb = pc;}
      else {a = c; pa = pc;}
    } while(fabs(pc) > tol);
    return c;
  }

  /**
   *   @brief  Returns the derivative of the polynomial
   *  
   *   @returns: q(x) = dp(x)/dx
   */    
  
  Polynomial deriv()
  {
    vector<double> d;    
    for(int i=0;i<c_.size()-1;i++)
      d.push_back(c_[i+1]*(i+1));
    return Polynomial(d);
  }
  
  /**
   *   @brief  Dermine the real roots of the polynomial
   *  
   *   @param: tol is the tolerence.
   *   @returns: array of real roots
   */    
  vector<double> real_roots(double tol)
  {
    vector< Polynomial > derivs;
    derivs.push_back(*this);
    while(derivs.back().size() > 2)
      derivs.push_back(derivs.back().deriv());
    
    vector<double> roots;
    for(auto rit = derivs.rbegin(); rit!= derivs.rend(); rit++)
      roots = (*rit).solve(roots, tol);
    
    return roots;  
  }

 protected:
  
  /**
   *   @brief  Solve for roots given inflection points
   *  
   *   @param: inflections is an array of inflection points (where p'(x) = 0)
   *   @param: tol is the tolerence
   *   @returns: array of roots
   */    
  vector<double> solve(vector<double> &inflections, double tol)
  {
    if(c_.size() < 2) throw std::runtime_error("c_.size < 2 ???");

    vector<double> roots;
  
    if(c_.size() == 2) // linear case: trivial
      {
	if(fabs(c_[1]) > tol) roots.push_back(-c_[0]/c_[1]);
	return roots;
      }

    // Now work with intervals [-inf, inflections[0], ..., inflections.back(), inf]
    double pinf = (c_.back() > 0 ? 1.0 : -1);
    double minf = (2*int((c_.size()+1)/2) == c_.size()+1 ? pinf : -pinf);
    
    if(inflections.size() == 0)
      {
	if(pinf*minf > 0) return roots;
	double b = 1; 
	while(eval(b)*pinf < 0) b *= 2;
	double a = -1;
	while(eval(a)*minf < 0) a *= 2;
      
	roots.push_back(bisect(a,b,tol));
      } else {
      // (-inf, i0]
      double pj = eval(inflections[0]);
      if(fabs(pj) < tol) roots.push_back(inflections[0]);
      else if(pj*minf < 0)
	{
	  double a = min(0.0,inflections[0])-1;
	  while(eval(a)*minf < 0) a *= 2;      
	  roots.push_back(bisect(a,inflections[0],tol));
	}
      for(int j=1;j<inflections.size();j++)
	{
	  double pjm = pj;
	  pj = eval(inflections[j]);
	  if(fabs(pj) < tol) roots.push_back(inflections[j]);
	  else if(pj*pjm < 0) roots.push_back(bisect(inflections[j-1], inflections[j],tol));
	}
      // (i_pmax, inf]
      if(pj*pinf < 0)
	{
	  double b = max(1.0,inflections.back())+1;
	  while(eval(b)*pinf < 0) b *= 2;      
	  roots.push_back(bisect(inflections.back(),b,tol));
	}      
    }
    return roots;
  }

private:
  vector<double> c_; ///< Array of coefficients defining the polynomial: p(x) = c[0] + c[1] x + c[2] x^2 +...
};

#endif // __LUTSKO_POLYNOMIAL__
