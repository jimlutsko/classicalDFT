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


//endregion