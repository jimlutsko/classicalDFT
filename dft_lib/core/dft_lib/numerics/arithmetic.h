#ifndef CLASSICALDFT_ARITHMETIC_H
#define CLASSICALDFT_ARITHMETIC_H

#include <numeric>
#include <vector>

namespace dft_core
{
namespace numerics
{
namespace arithmetic
{
namespace summation {

// region types:

enum class Type { KahanBabuska=0, KahanBabuskaNeumaier, KahanBabuskaKlein };

// endregion

// region free functions

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
std::tuple<double, double> KahanBabuskaSum(const std::vector<double>& x_series, const double& sum_ini = 0.0, const double& error_ini = 0.0);

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
std::tuple<double, double> KahanBabuskaNeumaierSum(const std::vector<double>& x_series, const double& sum_ini = 0.0, const double& error_ini = 0.0);

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
std::tuple<double, std::vector<double>> KahanBabuskaKleinSum(const std::vector<double>& x_series, const double& sum_ini = 0.0, const std::vector<double>& error_ini = {0, 0});

class CompensatedSum
{
  // region Constants:
  // endregion

  // region Attributes:
  // endregion

  // region Cttors:
  // endregion

  // region Methods:
  // endregion
};

}}}}

#endif  // CLASSICALDFT_ARITHMETIC_H
