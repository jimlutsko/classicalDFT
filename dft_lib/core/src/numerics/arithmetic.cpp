#include "dft_lib/numerics/arithmetic.h"

#include <vector>
#include <cmath>

namespace dft_core
{
namespace numerics
{
namespace arithmetic
{
namespace summation {

std::tuple<double, double> KahanBabuskaAdd(double x, double sum,  double error)
{
  auto x_corrected = x - error;
  auto t = sum + x_corrected;
  error = (t - sum) - x_corrected;
  sum = t;
  return std::tuple<double, double>{ sum, error };
}

std::tuple<double, double> KahanBabuskaSum(const std::vector<double>& x_series, const double& sum_ini, const double& error_ini)
{
  // The method is given in here:
  // wikipedia: https://en.wikipedia.org/wiki/Kahan_summation_algorithm#The_algorithm

  double sum = sum_ini;
  double error = error_ini;

  for (const auto& x : x_series)
  {
    std::tie(sum, error) = KahanBabuskaAdd(x, sum, error);
  }

  // When catching the result we need to use std::tie(x, y) [https://en.cppreference.com/w/cpp/utility/tuple]
  return std::tuple<double, double>{ sum, error };
}

std::tuple<double, double> KahanBabuskaNeumaierSum(const std::vector<double>& x_series, const double& sum_ini, const double& error_ini)
{
  // The method is given in here:
  // wikipedia: https://en.wikipedia.org/wiki/Kahan_summation_algorithm#Further_enhancements

  double sum = sum_ini;
  double error = error_ini;

  for (const auto& x : x_series)
  {
    auto t = sum + x;
    if (std::abs(sum) >= std::abs(x))
    {
      error += (sum - t) + x;
    }
    else
    {
      error += (x - t) + sum;
    }
    sum = t;
  }

  sum += error;
  return std::tuple<double, double>{ sum, error};
}

std::tuple<double, std::vector<double>> KahanBabuskaKleinSum(const std::vector<double>& x_series, const double& sum_ini, const std::vector<double>& error_ini)
{
  // The method is given in here:
  // wikipedia: https://en.wikipedia.org/wiki/Kahan_summation_algorithm#Further_enhancements

  double sum = sum_ini;
  std::vector<double> error = error_ini;

  for (const auto& x : x_series)
  {
    auto t = sum + x;
    double c;

    if (std::abs(sum) >= std::abs(x)) {
      c = (sum - t) + x;
    }
    else
    {
      c = (x - t) + sum;
    }

    double cc;
    sum = t;
    t = error[0] + c;

    if (std::abs(error[0]) >= std::abs(c))
    {
      cc = (error[0] - t) + c;
    }
    else
    {
      cc = (c - t) + error[0];
    }

    error[0] = t;
    error[1] += cc;
  }

  sum += (error[0] + error[1]);

  return std::tuple<double, std::vector<double>>{ sum, error };
}

}}}}