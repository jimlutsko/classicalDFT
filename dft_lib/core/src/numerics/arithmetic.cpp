#include "dft_lib/numerics/arithmetic.h"

#include <vector>

#include <dft_lib/utils/console.h>

namespace dft_core
{
namespace numerics
{
namespace arithmetic
{
namespace summation {

std::tuple<double, double> KahanBabuskaSum(std::vector<double>& x_series, const double& sum_ini, const double& error_ini)
{
  // The method is given in here:
  // wikipedia: https://en.wikipedia.org/wiki/Kahan_summation_algorithm#The_algorithm

  double sum = sum_ini;
  double error = error_ini;

  for (const auto& x : x_series)
  {
    auto x_corrected = x - error;
    auto t = sum + x_corrected;
    error = (t - sum) - x_corrected;
    sum = t;
  }

  // When catching the result we need to use std::tie(x, y) [https://en.cppreference.com/w/cpp/utility/tuple]
  return std::tuple<double, double>{ sum, error};
}


}}}}