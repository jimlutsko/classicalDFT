#ifndef __INTEGRATOR__
#define __INTEGRATOR__

#include <gsl/gsl_integration.h>


/** This class is a wrapper for the GSL integration 
 * routines qags and q
 *
 * Syntax for coding of member function pointers
 * is taken from 
 *  http://www.parashift.com/c++-faq-lite/pointers-to-members.html
 */

template <class T>
class Integrator
{
 public:
  Integrator(T* const obj , double (T::*fun)(double x) const, double epsRel = 1e-5, double epsAbs = 1e-5)
    : epsRel_(epsRel), epsAbs_(epsAbs), Nspace_(10000), obj_(obj), fun_(fun)
    {
      F_.function = &IntegratorHelper;
      F_.params = this;
      w1_ = gsl_integration_workspace_alloc(Nspace_);
//      w2_ = gsl_integration_cquad_workspace_alloc(Nspace_);
    }

  ~Integrator()
    {
      gsl_integration_workspace_free(w1_);
//      gsl_integration_cquad_workspace_free(w2_);
    }


  double integrateFinite(double a, double b) const
    {
      double result;
      gsl_integration_qags(&F_, a, b, epsRel_, epsAbs_, Nspace_, w1_, &result, &error_); 
      //      size_t nint = 0;
//      gsl_integration_cquad(&F_, a, b, epsRel_, epsAbs_, w2_, &result, &error_, &nint); 
      return result;
    }

  double integrateFiniteFast(double a, double b)
  {
      double result;
      size_t neval;
      gsl_integration_qng(&F_, a, b, epsAbs_, epsRel_, &result, &error_, &neval);
      return result;
  }

  double integrateSemiInfinite(double a) const
    {
      double result;
      gsl_integration_qagiu(&F_, a, epsRel_, epsAbs_, Nspace_, w1_, &result, &error_); 
      return result;
    }


  double toIntegrate(double x) const {return (obj_->*fun_)(x);}


  static double IntegratorHelper(double x, void* param) 
    {
      Integrator * i = (Integrator*) param;
      return i->toIntegrate(x);
    }


  double getError() const {return error_;}

 private:
  /// Bound on relative error
  double epsRel_;
  /// Bound on absolute error
  double epsAbs_;
  /// Error returned from gsl routine - access with getError function
  mutable double error_;

  /// gsl_function wrapper : This is only mutable because the gsl-calls won't take a (const gsl_function *)
  mutable gsl_function F_;
  /// Workspace for the gsl integration routines
  gsl_integration_workspace * w1_;
//  gsl_integration_cquad_workspace * w2_;
  /// Size of gsl workspace object
  int Nspace_;

  /// Object whose member function is to be integrated
  T* obj_;
  /// Member function to integrate
  double (T::*fun_)(double x) const;

};


class Integrator1 : public Integrator<Integrator1>
{
 public:
  Integrator1(double (*f)(double, void*), void *params = NULL, double epsRel = 1e-5, double epsAbs = 1e-5)
    : Integrator<Integrator1>(this, &Integrator1::helper, epsRel, epsAbs), func_(f), params_(params){}

  double helper(double r) const { return (*func_)(r, params_);}

 private:
  double (*func_)(double, void*);
  void * params_;
};


template <class T>
class IntegratorFlag
{
 public:
IntegratorFlag(T* const obj , double (T::*fun)(double x, int flag) const, double epsRel = 1e-5, double epsAbs = 1e-5)
    : epsRel_(epsRel), epsAbs_(epsAbs), Nspace_(10000), obj_(obj), fun_(fun)
    {
      F_.function = &IntegratorHelper;
      F_.params = this;
      w1_ = gsl_integration_workspace_alloc(Nspace_);
//      w2_ = gsl_integration_cquad_workspace_alloc(Nspace_);
    }

  ~IntegratorFlag()
    {
      gsl_integration_workspace_free(w1_);
//      gsl_integration_cquad_workspace_free(w2_);
    }


  double integrateFinite(double a, double b, int flag) const
    {
	flag_ = flag;
      double result;
      gsl_integration_qags(&F_, a, b, epsRel_, epsAbs_, Nspace_, w1_, &result, &error_); 
      //      size_t nint = 0;
//      gsl_integration_cquad(&F_, a, b, epsRel_, epsAbs_, w2_, &result, &error_, &nint); 
      return result;
    }

  double integrateFiniteFast(double a, double b, int flag) const
  {
      double result;
      size_t neval;
      flag_ = flag;
      gsl_integration_qng(&F_, a, b, epsAbs_, epsRel_, &result, &error_, &neval);
      return result;
  }


  double integrateSemiInfinite(double a, int flag) const
  {
      flag_ = flag;
      double result;
      gsl_integration_qagiu(&F_, a, epsRel_, epsAbs_, Nspace_, w1_, &result, &error_); 
      return result;
    }


  double toIntegrate(double x) const {return (obj_->*fun_)(x, flag_);}

  static double IntegratorHelper(double x, void* param) 
    {
      IntegratorFlag * i = (IntegratorFlag*) param;
      return i->toIntegrate(x);
    }


  double getError() const {return error_;}

 private:
  /// Bound on relative error
  double epsRel_;
  /// Bound on absolute error
  double epsAbs_;
  /// Error returned from gsl routine - access with getError function
  mutable double error_;

  /// gsl_function wrapper : This is only mutable because the gsl-calls won't take a (const gsl_function *)
  mutable gsl_function F_;
  /// Workspace for the gsl integration routines
  gsl_integration_workspace * w1_;
//  gsl_integration_cquad_workspace * w2_;
  /// Size of gsl workspace object
  int Nspace_;

  /// Object whose member function is to be integrated
  T* obj_;
  /// Member function to integrate
  double (T::*fun_)(double x, int flag) const;

  mutable int flag_;
};






#endif // sentinal
