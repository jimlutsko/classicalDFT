#ifndef CLASSICALDFT_INTEGRATION_H
#define CLASSICALDFT_INTEGRATION_H

#include <gsl/gsl_integration.h>
#include <functional>
#include <memory>
#include "dft_lib/utils/console.h"

namespace dft_core
{
namespace numerics
{

//region Constants:

// Convenient definitions for the default initialisation of standard parameters:
const double DEFAULT_RELATIVE_ERROR_TOLERANCE = 1e-5;
const double DEFAULT_ABSOLUTE_ERROR_TOLERANCE = 1e-5;
const double DEFAULT_INITIAL_ERROR_VALUE = 123456.789;
const double DEFAULE_INITIAL_RESULT_VALUE = 123456.789;
const int DEFAULT_GSL_WORKING_SPACE_SIZE = 1000;

//endregion

//region Alias:

template<class T, class return_type, class input_type = return_type>
using class_method = std::function<return_type(const T&, input_type)>;

//endregion

//region Integrator class:

/**
 * @brief  UTILITY: wrapper for GSL integration routines.

 * @detailed Provides a simplified, c++ interface for using some of the GSL integration routines.
 *  This version is intended to be the parent of an object (of type `T`) with a member
 *  function we wish to integrate, say `x_type func(x_type x)`, with x_type = double.
 *  The following integration methods are wrapped:
 *
 *  - gsl_integration_qags: Quadrature Adaptive General-integrand with Singularities (QAGS).
 *    This integrator uses the Gauss-Kronrod 21-point rule.
 *
 *  - gsl_integration_qng: Quadrature Non-adaptive General-integrand (QNG).
 *    The integrator uses Gauss-Kronrod 10-point, 21-point, 43-point and 87-point integration rules.
 *
 *  - gsl_integration_qagiu: Quadrature Adaptive General-integration on Infinite interval (QAGI).
 *    This integrator works on the positive semi-infinite axis [0,+inf].
 */
template<class T, class x_type = double>
class Integrator
{
 private:
  //region Attributes:

  /**
   * @brief Determines the number of double precision intervals, their integration results
   * and error estimates, needed to be allocated and saved.
   */
  int working_space_size_ = DEFAULT_GSL_WORKING_SPACE_SIZE;

  /**
   * @brief The relative error tolerance to be used as stopping criterion
   */
  double relative_error_ = DEFAULT_RELATIVE_ERROR_TOLERANCE;

  /**
   * @brief The absolute error tolerance to be used as stopping criterion
   */
  double absolute_error_ = DEFAULT_ABSOLUTE_ERROR_TOLERANCE;

  /**
   * @brief The numerical error achieved with the integration method
   * @remark The `mutable` keyword is required to allow gsl_function the updating of this
   * attribute within a const environment
   */
  mutable double error_ = DEFAULT_INITIAL_ERROR_VALUE;

  /**
   * @brief The numerical result achieved with the integration method
   * @remark The `mutable` keyword is required to allow gsl_function the updating of this
   * attribute within a const environment
   */
  mutable double result_ = DEFAULE_INITIAL_RESULT_VALUE;

  /**
   * @brief (const) Reference to the problem object which contains the method to be integrated
   */
  const T& problem_;

  /**
   * @brief Stores the callable object (`function`) of a given problem object (of type `T`)
   * @remark For more information about `std::function` see https://en.cppreference.com/w/cpp/utility/functional/function
   */
  class_method<T, x_type> method_;

  //region GSL attributes:

  /**
   * @brief The gsl_function structure to pass to the integration method (https://www.gnu.org/software/gsl/doc/html/roots.html?highlight=gsl_function#c.gsl_function)
   */
  gsl_function gsl_function_{};

  /**
   * @brief Smart pointer which will hold the `gsl_integration_workspace` pointer required
   * by GSL integrators (https://www.gnu.org/software/gsl/doc/html/integration.html?highlight=gsl_integration_workspace#c.gsl_integration_workspace)
   */
  std::unique_ptr<gsl_integration_workspace> working_space_;

  //endregion
  //endregion

  //region Initialisers:

  /**
   * @brief Initialises the `private: gsl_function` structure with:
   *    - function -> integrand_function
   *    - params -> this
   */
  void InitialiseGslFunction()
  {
    gsl_function_.function = &this->integrand_function;
    gsl_function_.params = this;
  }

  /**
   * @brief Initialises the `private: working_space_` pointer by calling `gsl_integration_workspace_alloc`
   * (https://www.gnu.org/software/gsl/doc/html/integration.html?highlight=gsl_integration_workspace_alloc#c.gsl_integration_workspace_alloc)
   */
  void InitialiseGslWorkSpace()
  {
    this->working_space_.reset(gsl_integration_workspace_alloc(working_space_size_));
  }

  //endregion

 public:

  //region Cttors:

  /**
   * @brief Cttor
   *
   * @param problem
   * @param method
   * @param relative_error
   * @param absolute_error
   * @param work_space_size
   */
  explicit Integrator(
      const T& problem,
      class_method<T, x_type> method,
      double relative_error = DEFAULT_RELATIVE_ERROR_TOLERANCE,
      double absolute_error = DEFAULT_ABSOLUTE_ERROR_TOLERANCE,
      int work_space_size = DEFAULT_GSL_WORKING_SPACE_SIZE
  ): problem_(problem), method_(method), working_space_size_(work_space_size)
  {
    this->InitialiseGslWorkSpace();
    this->InitialiseGslFunction();
  }

  //endregion

  //region Inspectors:

  ///Gets the numerical error achieved with the integration method
  double numerical_error() const { return error_; };
  ///Gets the numerical result achieved with the integration method
  double numerical_result() const { return result_; };

  //endregion

  //region Mutators:

  void SetAbsoluteError(double absolute_error) { absolute_error_ = absolute_error; }
  void SetRelativeError(double relative_error) { relative_error_ = relative_error; }
  void SetWorkingSpaceSize(int work_space_size) { working_space_size_ = work_space_size; }

  //endregion

  //region Methods:

  /**
   * @brief Binds the method pointer to the `problem_` object and evaluates it at `x`
   * @param x The x-value where the method will be evaluated at
   * @return The result of evaluating the method at `x`
   */
  x_type function(x_type x) const
  {
    return method_(problem_, x);
  }

  /**
   * @brief The integrand function to pass to gsl methods, which needs to have a signature: (x, *params).
   *
   * @param x The x-value where the function will be evaluated at
   * @param param A void-ptr which is needed for legacy reasons
   *    This needs to be the reference to the Integrator object we are working with
   * @return: x_type = double
   *    The value of the integrand function evaluated at `x`
   */
  static x_type integrand_function(x_type x, void* param)
  {
    auto integrator_obj = (Integrator *) param;
    return integrator_obj->function(x);
  }

  /**
   * @brief Performs the numerical integration of the `integrand_function`on the interval
   *    [`limit_inferior`, `limit_superior`] using an adaptive method
   *
   * @param limit_inferior The lower limit of the integral
   * @param limit_superior The upper limit of the integral
   * @return The numerical approximation of the integral
   */
  double DefiniteIntegral(x_type limit_inferior, x_type limit_superior) const
  {
    // Quadrature Adaptive General-integrand with Singularities (QAGS):
    // For mor information: https://www.gnu.org/software/gsl/doc/html/integration.html?highlight=gsl_integration_qags#c.gsl_integration_qags
    gsl_integration_qags(
        &gsl_function_,
        limit_inferior,
        limit_superior,
        this->absolute_error_,
        this->relative_error_,
        this->working_space_size_,
        this->working_space_.get(),
        &this->result_,
        &this->error_
    );

    return this->result_;
  }

  /**
   * @brief Performs the numerical integration of the `integrand_function`on the interval
   *    [`limit_inferior`, `limit_superior`] using a non-adaptive method
   *
   * @param limit_inferior The lower limit of the integral
   * @param limit_superior The upper limit of the integral
   * @return The numerical approximation of the integral
   */
  double DefiniteIntegralFast(x_type limit_inferior, x_type limit_superior) const
  {
    size_t number_of_evaluations = 0;
    // Quadrature Non-adaptive General-integrand (QNG):
    // For mor information: https://www.gnu.org/software/gsl/doc/html/integration.html?highlight=gsl_integration_qng#c.gsl_integration_qng
    gsl_integration_qng(
        &gsl_function_,
        limit_inferior,
        limit_superior,
        this->absolute_error_,
        this->relative_error_,
        &this->result_,
        &this->error_,
        &number_of_evaluations
    );

    return this->result_;
  }

  /**
   * @brief Performs the numerical integration of the `integrand_function`on the interval
   *    [`limit_inferior`, +inf]
   *
   * @param limit_inferior The lower limit of the integral
   * @return The numerical approximation of the integral
   */
  double UpperSemiInfiniteIntegral(x_type limit_inferior)
  {
    // Quadrature Adaptive General-integration on Infinite Upper interval (QAGIU):
    // For more information: https://www.gnu.org/software/gsl/doc/html/integration.html?highlight=gsl_integration_qagiu#c.gsl_integration_qagiu
    gsl_integration_qagiu(
        &gsl_function_,
        limit_inferior,
        this->absolute_error_,
        this->relative_error_,
        this->working_space_size_,
        this->working_space_.get(),
        &this->result_,
        &this->error_
    );

    return this->result_;
  }

  /**
   * @brief Performs the numerical integration of the `integrand_function`on the interval
   *    [-inf, `limit_superior`]
   *
   * @param limit_superior The upper limit of the integral
   * @return The numerical approximation of the integral
   */
  double LowerSemiInfiniteIntegral(x_type limit_superior)
  {
    // Quadrature Adaptive General-integration on Infinite Lower interval (QAGIL):
    // For more information: https://www.gnu.org/software/gsl/doc/html/integration.html?highlight=gsl_integration_qagiu#c.gsl_integration_qagil
    gsl_integration_qagil(
        &gsl_function_,
        limit_superior,
        this->absolute_error_,
        this->relative_error_,
        this->working_space_size_,
        this->working_space_.get(),
        &this->result_,
        &this->error_
    );

    return this->result_;
  }

  /**
   * @brief Performs the numerical integration of the `integrand_function`on the interval [-inf, +inf]
   *
   * @return The numerical approximation of the integral
   */
  double FullInfiniteIntegral()
  {
    // Quadrature Adaptive General-integration on Infinite Lower interval (QAGIL):
    // For more information: https://www.gnu.org/software/gsl/doc/html/integration.html?highlight=gsl_integration_qagiu#c.gsl_integration_qagil
    gsl_integration_qagi(
        &gsl_function_,
        this->absolute_error_,
        this->relative_error_,
        this->working_space_size_,
        this->working_space_.get(),
        &this->result_,
        &this->error_
    );

    return this->result_;
  }

  //endregion
};


//endregion

}}

#endif  // CLASSICALDFT_INTEGRATION_H
