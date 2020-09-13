#include <gtest/gtest.h>

#include "dft_lib/utils/functions.h"

#include <boost/range/combine.hpp>

namespace example
{
class TestClass
{
 public:
  TestClass() = default;
  double fn_cubic(double x) const
  {
    return x*x*x;
  }
};

double fn_sqr(double x)
{
  return x*x;
}
}

using namespace dft_core::utils::functions;

//region Methods
TEST(functions, apply_vector_wise_fn_works_ok)
{
  auto x_vec = std::vector<double>{1.0, 2.0, 3.0, 1.765};
  auto actual = apply_vector_wise<double>(&example::fn_sqr, x_vec);
  auto expected = std::vector<double>{1.0, 4.0, 9.0, 3.115225};
  for (auto tup : boost::combine(expected, actual))
  {
    double x, y;
    boost::tie(x,y) = tup;
    EXPECT_DOUBLE_EQ(x, y);
  }
}

TEST(functions, apply_vector_wise_method_works_ok)
{
  auto x_vec = std::vector<double>{1.0, 2.0, 3.0, 1.765};
  auto obj = example::TestClass();
  auto actual = apply_vector_wise<example::TestClass, double>(obj, &example::TestClass::fn_cubic, x_vec);
  auto expected = std::vector<double>{1.0, 8.0, 27.0, 5.498372125};
  for (auto tup : boost::combine(expected, actual))
  {
    double x, y;
    boost::tie(x,y) = tup;
    EXPECT_DOUBLE_EQ(x, y);
  }
}
//endregion