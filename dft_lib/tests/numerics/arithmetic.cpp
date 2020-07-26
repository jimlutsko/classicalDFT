#include <gtest/gtest.h>
#include <cmath> // std::abs

#include "dft_lib/numerics/arithmetic.h"

//region free functions:

using namespace dft_core::numerics;
using error_type = dft_core::numerics::arithmetic::summation::error_type;
using x_type = dft_core::numerics::arithmetic::summation::x_type;

TEST(numerics, kahan_babuska_sum_works_ok)
{
  auto x_input = std::vector<double>{1.0, 2.0, 4.0, 3.0};

  int expected = 10;
  x_type actual;
  error_type err;
  std::tie(actual, err) = arithmetic::summation::KahanBabuskaSum(x_input);
  EXPECT_DOUBLE_EQ(expected,actual);
  EXPECT_LE(std::abs(std::abs(expected)-std::abs(actual)), std::abs(err.front()));
}

TEST(numerics, kahan_babuska_accuracy_ok)
{
  auto x_input = std::vector<double>{1.0, 2.0, 4.0, 3.0};

  int expected = 10;
  x_type actual;
  error_type err;
  std::tie(actual, err) = arithmetic::summation::KahanBabuskaSum(x_input);
  EXPECT_LE(std::abs(err.front()), 2*DBL_EPSILON*expected);
}

TEST(numerics, kahan_babuska_neumaier_sum_works_ok)
{
  auto x_input = std::vector<double>{1.0, 2.0, 4.0, 3.0};

  int expected = 10;
  x_type actual;
  error_type err;
  std::tie(actual, err) = arithmetic::summation::KahanBabuskaNeumaierSum(x_input);
  EXPECT_DOUBLE_EQ(expected,actual);
  EXPECT_LE(std::abs(std::abs(expected)-std::abs(actual)), std::abs(err.front()));
}

TEST(numerics, kahan_babuska_neumaier_accuracy_ok)
{
  auto x_input = std::vector<double>{1.0, 2.0, 4.0, 3.0};

  int expected = 10;
  x_type actual;
  error_type err;

  std::tie(actual, err) = arithmetic::summation::KahanBabuskaNeumaierSum(x_input);
  EXPECT_LE(std::abs(err.front()), 2*DBL_EPSILON*expected);
}

TEST(numerics, kahan_babuska_klein_sum_works_ok)
{
  auto x_input = std::vector<double>{1.0, 2.0, 4.0, 3.0};

  int expected = 10.0;
  x_type actual;
  error_type err;
  std::tie(actual, err) = arithmetic::summation::KahanBabuskaKleinSum(x_input);
  EXPECT_DOUBLE_EQ(expected, actual);
}

TEST(numerics, kahan_babuska_klein_accuracy_ok)
{
  auto x_input = std::vector<double>{1.0, 2.0, 4.0, 3.0};

  int expected = 10;
  x_type actual;
  error_type err;
  std::tie(actual, err) = arithmetic::summation::KahanBabuskaKleinSum(x_input);
  EXPECT_LE(std::abs(err[0])+std::abs(err[1]), 2*DBL_EPSILON*expected);
}

TEST(numerics, compensated_summation_default_works_ok)
{
  auto x_input = std::vector<double>{1.0, 2.0, 4.0, 3.0};

  x_type expected = 10;
  auto actual = arithmetic::summation::CompensatedSum();

  actual += x_input;
  EXPECT_DOUBLE_EQ(expected, actual.Sum());
  EXPECT_LE(std::abs(actual.error().front()), 2*DBL_EPSILON*expected);
}

TEST(numerics, compensated_summation_kb_works_ok)
{
  auto x_input = std::vector<double>{1.0, 2.0, 4.0, 3.0};

  x_type expected = 10;
  auto actual = arithmetic::summation::CompensatedSum(arithmetic::summation::Type::KahanBabuska);

  actual += x_input;
  EXPECT_DOUBLE_EQ(expected, actual.Sum());
  EXPECT_LE(std::abs(actual.error().front()), 2*DBL_EPSILON*expected);
}

TEST(numerics, compensated_summation_kbn_works_ok)
{
  auto x_input = std::vector<double>{1.0, 2.0, 4.0, 3.0};

  x_type expected = 10;
  auto actual = arithmetic::summation::CompensatedSum(arithmetic::summation::Type::KahanBabuskaNeumaier);

  actual += x_input;
  EXPECT_DOUBLE_EQ(expected, actual.Sum());
  EXPECT_LE(std::abs(actual.error().front()), 2*DBL_EPSILON*expected);
}

TEST(numerics, compensated_summation_kbk_works_ok)
{
  auto x_input = std::vector<double>{1.0, 2.0, 4.0, 3.0};

  x_type expected = 10;
  auto actual = arithmetic::summation::CompensatedSum(arithmetic::summation::Type::KahanBabuskaKlein);

  actual += x_input;
  EXPECT_DOUBLE_EQ(expected, actual.Sum());
  EXPECT_LE(std::abs(actual.error()[0])+std::abs(actual.error()[1]), 2*DBL_EPSILON*expected);
}
//endregion