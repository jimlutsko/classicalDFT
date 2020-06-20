#include <gtest/gtest.h>
#include <cmath> // std::abs

#include "dft_lib/numerics/arithmetic.h"

//region free functions:

using namespace dft_core::numerics;

TEST(numerics, kahan_babuska_sum_works_ok)
{
  auto x_input = std::vector<double>{1.0, 2.0, 4.0, 3.0};

  double expected = 10.0;
  double actual, err;
  std::tie(actual, err) = arithmetic::summation::KahanBabuskaSum(x_input);
  EXPECT_DOUBLE_EQ(expected,actual);
  EXPECT_LE(std::abs(std::abs(expected)-std::abs(actual)), std::abs(err));
}

TEST(numerics, kahan_babuska_accuracy_ok)
{
  auto x_input = std::vector<double>{1.0, 2.0, 4.0, 3.0};

  int expected = 10;
  double actual, err;
  std::tie(actual, err) = arithmetic::summation::KahanBabuskaSum(x_input);
  std::cout << err << std::endl;
  std::cout << DBL_EPSILON << std::endl;
  std::cout << 2*DBL_EPSILON*expected << std::endl;
  EXPECT_LE(std::abs(err), 2*DBL_EPSILON*expected);
}

//endregion