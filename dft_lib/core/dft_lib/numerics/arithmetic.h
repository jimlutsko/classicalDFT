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

enum class Type { Kahan = 0, KahanBabuska, KahanBabuskaNeumaier, KahanBabuskaKlein };

// endregion

// region free functions

template <class x_type=double>
x_type StandardVectorSum(std::vector<x_type> x_input)
{
  x_type result = std::accumulate(x_input.begin(), x_input.end(), static_cast<x_type>(0.0));
  return result;
}

std::tuple<double, double> KahanBabuskaSum(std::vector<double>& x_series, const double& sum_ini = 0.0, const double& error_ini = 0.0);
std::tuple<double, double> KahanBabuskaNeumaierSum(std::vector<double>& x_series, const double& sum_ini = 0.0, const double& error_ini = 0.0);
std::tuple<double, double> KahanBabuskaKleinSum(std::vector<double>& x_series, const double& sum_ini = 0.0, const double& error_ini = 0.0);

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
