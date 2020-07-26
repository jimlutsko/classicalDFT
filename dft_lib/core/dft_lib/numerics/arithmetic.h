#ifndef CLASSICALDFT_ARITHMETIC_H
#define CLASSICALDFT_ARITHMETIC_H

#include <numeric>
#include <vector>
#include <tuple>
#include <functional>

#include "dft_lib/utils/console.h"

namespace dft_core
{
namespace numerics
{
namespace arithmetic
{
namespace summation
{

// region types:

enum class Type { KahanBabuska=0, KahanBabuskaNeumaier, KahanBabuskaKlein };

// endregion

// region alias:
using x_type = double;
using error_type = std::vector<x_type>;
using return_type = std::tuple<x_type, error_type>;

template<class return_type = return_type, class input_type = x_type, class sum_type = x_type, class err_type = error_type>
using summation_method = std::function<return_type(input_type, x_type, err_type)>;
// endregion

// region free functions:

/**
 * @brief Returns the standard sum of a `std::vector<T>`'s components
 *
 * @detailed The STL's vector class is quite useful as a container, but sometimes is
 *  insufficient as a numerical class, given that it does not provide basic numerical capabilities.
 *  The reason for this being that the vector class is nothing but a container class,
 *  despite its name. This free-function intends to provide the most basic functionality,
 *  that of adding up the components of a container when such components are numerical values.
 *
 * @tparam x_type the type of the components of the container
 * @param x_input the input std::vector
 * @return the standard sum of the components carried out via the std::accumulate algorithm,
 *  which uses the natural `+` operator.
 */
template <
    typename x_type=double,
    // This meta-function makes this function only available for numerical types
    typename = typename std::enable_if<std::is_arithmetic<x_type>::value, x_type>::type
>
x_type StandardVectorSum(const std::vector<x_type>& x_input)
{
  x_type result = std::accumulate(x_input.begin(), x_input.end(), static_cast<x_type>(0.0));
  return result;
}

/**
 * @brief Standard implementation of the Kahan-Babuška summation algorithm
 *
 * @detailed To deep dive into the numerical details of the K-B summation algorithm we
 *  recommend to read the wiki link: https://en.wikipedia.org/wiki/Kahan_summation_algorithm.
 *  There you can find multiple and useful references which provide details about the formal
 *  derivation, precision, convergence theorems etc.
 *
 * @param x_series the std::vector<double> we want to sum
 * @param sum_ini an initial value for the summation algorithm (default = 0.0)
 *  This is useful when the summation is carried out in several steps, so that we can
 *  continue the summation from where it was left.
 * @param error_ini an initial value for the compensation (default = 0.0).
 *  This is useful when the summation is carried out in several steps, so that we can
 *  continue the summation from where it was left.
 *
 * @return tuple of double values which contains:
 *  1) the sum value, and
 *  2) the status of the compensation in case we need to continue the sum later
 */
return_type KahanBabuskaSum(const std::vector<x_type>& x_series, const x_type& sum_ini = 0.0, const error_type & error_ini = {0.0});
return_type KahanBabuskaAdd(x_type x, x_type sum, error_type error);


/**
 * @brief Standard implementation of the Kahan-Babuška-Neumaier summation algorithm (https://www.mat.univie.ac.at/~neum/scan/01.pdf)
 *
 * @detailed To deep dive into the numerical details of the K-B-N summation algorithm we
 *  recommend to read the wiki link: https://en.wikipedia.org/wiki/Kahan_summation_algorithm#Further_enhancements.
 *  There you can find multiple and useful references which provide details about the formal
 *  derivation, accuracy, convergence theorems etc.
 *
 * @param x_series the std::vector<double> we want to sum
 * @param sum_ini an initial value for the summation algorithm (default = 0.0)
 *  This is useful when the summation is carried out in several steps, so that we can
 *  continue the summation from where it was left.
 * @param error_ini an initial value for the compensation (default = 0.0).
 *  This is useful when the summation is carried out in several steps, so that we can
 *  continue the summation from where it was left.
 *
 * @return tuple of double values which contains:
 *  1) the sum value, and
 *  2) the status of the compensation in case we need to continue the sum later
 */
return_type KahanBabuskaNeumaierSum(const std::vector<x_type>& x_series, const x_type& sum_ini = 0.0, const error_type & error_ini = {0.0});
return_type KahanBabuskaNeumaierAdd(x_type x, x_type sum, error_type error);

/**
 * @brief Standard implementation of the Kahan-Babuška-Klein summation algorithm (https://link.springer.com/article/10.1007/s00607-005-0139-x)
 *
 * @detailed To deep dive into the numerical details of the K-B-K summation algorithm we
 *  recommend to read the wiki link: https://en.wikipedia.org/wiki/Kahan_summation_algorithm#Further_enhancements.
 *  There you can find multiple and useful references which provide details about the formal
 *  derivation, accuracy, convergence theorems etc.
 *
 * @param x_series the std::vector<double> we want to sum
 * @param sum_ini an initial value for the summation algorithm (default = 0.0)
 *  This is useful when the summation is carried out in several steps, so that we can
 *  continue the summation from where it was left.
 * @param error_ini an initial value for the compensation (default = 0.0).
 *  This is useful when the summation is carried out in several steps, so that we can
 *  continue the summation from where it was left.
 *
 * @return tuple of double values which contains:
 *  1) the sum value, and
 *  2) the status of the compensation in case we need to continue the sum later
 */
return_type KahanBabuskaKleinSum(const std::vector<x_type>& x_series, const x_type& sum_ini = 0.0, const error_type & error_ini = {0.0, 0.0});
return_type KahanBabuskaKleinAdd(x_type x, x_type sum, error_type error);

// endregion

class CompensatedSum {
 private:
  // region Constants:
  // endregion

  // region Attributes:
  /**
   * @brief The cumulative-sum value
   */
  x_type sum_ = 0.0;

  /**
   * @brief The state of the correction components needed
   */
  error_type error_;

  /**
   * @brief The type of compensated sum to be used, from the enum class `summation::Type`.
   *    By default set to `Type::KahanBabuskaNeumaier`
   */
  summation::Type type_ = summation::Type::KahanBabuskaNeumaier;
  // endregion

  /**
   * @brief The underlying free method to be used (either Kahan-Babuška, Kahan-Babuška-Neumaier,
   * or Kahan-Babuška-Klein). By default = Kahan-Babuška-Neumaier
   */
  summation_method<return_type, x_type, x_type, error_type> method_;

  /**
   * @brief A helper function which basically calls the underlying `method_` to carry out the sum
   * @param value The value to be added to the sum
   */
  void add(x_type value);

 public:
  // region Cttors:
  /**
   * @brief A wrapper for the different Kahan-Babuška summation algorithms
   * @param type Specifies the specific variation of the Kahan-Babuška algorithm to be used.
   *    By default the type of compensated summation is set to `KahanBabuskaNeumaier`
   */
  explicit CompensatedSum(summation::Type type = summation::Type::KahanBabuskaNeumaier);
  // endregion

  // region Inspectors:
  /// Gets the error_ vector
  const error_type& error() const;
  /// Gets the sum type_
  const summation::Type& type() const;
  // endregion

  // region Methods:
  /// Gets the value of the sum given the current error-vector state
  x_type Sum() const;
  // endregion

  // region Overloads:
  explicit operator double() const;

  CompensatedSum& operator+=(x_type x);
  CompensatedSum& operator-=(x_type x);

  CompensatedSum& operator+=(const std::vector<x_type>& x);
  CompensatedSum& operator-=(const std::vector<x_type>& x);

  CompensatedSum& operator+=(const CompensatedSum& other);
  CompensatedSum& operator-=(const CompensatedSum& other);

  CompensatedSum& operator+(x_type x);
  CompensatedSum& operator-(x_type x);

  CompensatedSum& operator+(const CompensatedSum& other);
  CompensatedSum& operator-(const CompensatedSum& other);

  CompensatedSum& operator=(x_type x);

  friend std::ostream& operator<<(std::ostream &output, const CompensatedSum&obj);
  // endregion
};

}}}}

#endif  // CLASSICALDFT_ARITHMETIC_H
