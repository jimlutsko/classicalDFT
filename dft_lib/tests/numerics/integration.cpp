#include <gtest/gtest.h>

#include "dft_lib/numerics/integration.h"

//region Default parameters:

using namespace dft_core::numerics;

TEST(numerics, default_relative_error_test)
{
  EXPECT_DOUBLE_EQ(
      1e-5,
      integration::DEFAULT_RELATIVE_ERROR_TOLERANCE
  );
}

TEST(numerics, default_absolute_error_test)
{
  EXPECT_DOUBLE_EQ(
      1e-5,
      integration::DEFAULT_ABSOLUTE_ERROR_TOLERANCE
  );
}


TEST(numerics, default_initial_error_value_test)
{
  EXPECT_DOUBLE_EQ(
      123456.789,
      integration::DEFAULT_INITIAL_ERROR_VALUE
  );
}

TEST(numerics, default_initial_result_value_test)
{
  EXPECT_DOUBLE_EQ(
      123456.789,
      integration::DEFAULE_INITIAL_RESULT_VALUE
  );
}


TEST(numerics, default_gsl_working_space_size_test)
{
  EXPECT_EQ(
      1000,
      integration::DEFAULT_GSL_WORKING_SPACE_SIZE
  );
}

//endregion


class TestProblemClass
{
 private:
  double param_ = 0.5;

 public:
  explicit TestProblemClass(double param = 0.5): param_(param) {}
  double NegativeExp(double x) const { return param_ * exp(-x); }
  double PositiveExp(double x) const { return param_ * exp(x); }
  double NormalDist(double x) const { return param_ * exp(-x*x*0.5)/sqrt(2*M_PI); }
};

//region Integrator class:
TEST(integrator, cttor_works_ok_test)
{
  auto problem = TestProblemClass(1.0);
  auto integrator = integration::Integrator<TestProblemClass>(problem, &TestProblemClass::NegativeExp);

  EXPECT_DOUBLE_EQ(
      integrator.absolute_error_tolerance(),
      integration::DEFAULT_ABSOLUTE_ERROR_TOLERANCE
  );

  EXPECT_DOUBLE_EQ(
      integrator.relative_error_tolerance(),
      integration::DEFAULT_RELATIVE_ERROR_TOLERANCE
  );

  EXPECT_DOUBLE_EQ(
      integrator.numerical_error(),
      integration::DEFAULT_INITIAL_ERROR_VALUE
  );

  EXPECT_DOUBLE_EQ(
      integrator.numerical_result(),
      integration::DEFAULE_INITIAL_RESULT_VALUE
  );

  EXPECT_EQ(
      integrator.gsl_working_space_size(),
      integration::DEFAULT_GSL_WORKING_SPACE_SIZE
  );
}

TEST(integrator, setters_works_ok_test) {
  auto problem = TestProblemClass(1.0);
  auto integrator = integration::Integrator<TestProblemClass>(problem, &TestProblemClass::NegativeExp);

  auto expected_val = 1e-7;
  integrator.SetAbsoluteError(expected_val);
  EXPECT_DOUBLE_EQ(
      integrator.absolute_error_tolerance(),
      expected_val
  );

  integrator.SetRelativeError(expected_val);
  EXPECT_DOUBLE_EQ(
      integrator.relative_error_tolerance(),
      expected_val
  );

  auto expected_size = 15000;
  integrator.SetWorkingSpaceSize(expected_size);
  EXPECT_EQ(
      integrator.gsl_working_space_size(),
      expected_size
  );
}

TEST(integrator, definite_integral_works_ok_test)
{
  auto problem = TestProblemClass(1.0);
  auto integrator = integration::Integrator<TestProblemClass>(problem, &TestProblemClass::NegativeExp);

  auto expected_result = 0.5;
  auto actual_result = integrator.DefiniteIntegral(0, -log(0.5));

  EXPECT_DOUBLE_EQ(
      expected_result,
      actual_result
  );
  EXPECT_DOUBLE_EQ(
      expected_result,
      integrator.numerical_result()
  );
  EXPECT_LE(
      integrator.numerical_error(),
      integrator.relative_error_tolerance()
  );
}

TEST(integrator, definite_integral_fast_works_ok_test)
{
  auto problem = TestProblemClass(1.0);
  auto integrator = integration::Integrator<TestProblemClass>(problem, &TestProblemClass::NegativeExp);

  auto expected_result = 0.5;
  auto actual_result = integrator.DefiniteIntegralFast(0, -log(0.5));

  EXPECT_DOUBLE_EQ(
      expected_result,
      actual_result
  );
  EXPECT_DOUBLE_EQ(
      expected_result,
      integrator.numerical_result()
  );
  EXPECT_LE(
      integrator.numerical_error(),
      integrator.relative_error_tolerance()
  );
}

TEST(integrator, upper_semi_infinite_integral_works_ok_test)
{
  auto problem = TestProblemClass(1.0);
  auto integrator = integration::Integrator<TestProblemClass>(problem, &TestProblemClass::NegativeExp);

  double expected_result = 1.0;
  auto actual_result = integrator.UpperSemiInfiniteIntegral(0);

  EXPECT_NEAR(
      expected_result,
      actual_result,
      integrator.numerical_error()
  );
  EXPECT_NEAR(
      expected_result,
      integrator.numerical_result(),
      integrator.numerical_error()
  );
  EXPECT_LE(
      integrator.numerical_error(),
      integrator.relative_error_tolerance()
  );
}

TEST(integrator, lower_semi_infinite_integral_works_ok_test)
{
  auto problem = TestProblemClass(1.0);
  auto integrator = integration::Integrator<TestProblemClass>(problem, &TestProblemClass::PositiveExp);

  double expected_result = 1.0;
  auto actual_result = integrator.LowerSemiInfiniteIntegral(0);

  EXPECT_NEAR(
      expected_result,
      actual_result,
      integrator.numerical_error()
  );
  EXPECT_NEAR(
      expected_result,
      integrator.numerical_result(),
      integrator.numerical_error()
  );
  EXPECT_LE(
      integrator.numerical_error(),
      integrator.relative_error_tolerance()
  );
}

TEST(integrator, infinite_integral_works_ok_test)
{
  auto problem = TestProblemClass(1.0);
  auto integrator = integration::Integrator<TestProblemClass>(problem, &TestProblemClass::NormalDist);

  double expected_result = 1.0;
  auto actual_result = integrator.FullInfiniteIntegral();

  EXPECT_NEAR(
      expected_result,
      actual_result,
      integrator.numerical_error()
  );
  EXPECT_NEAR(
      expected_result,
      integrator.numerical_result(),
      integrator.numerical_error()
  );
  EXPECT_LE(
      integrator.numerical_error(),
      integrator.relative_error_tolerance()
  );
}

TEST(integrator, class_method_passing_works_ok_test) {
  auto problem = TestProblemClass(1.0);
  auto integrator = integration::Integrator<TestProblemClass>(problem, &TestProblemClass::NormalDist);

  double expected_value = problem.NormalDist(2*M_PI);
  EXPECT_DOUBLE_EQ(
      integrator.function(2*M_PI),
      expected_value
  );
  EXPECT_DOUBLE_EQ(
      integration::Integrator<TestProblemClass>::integrand_function(2*M_PI, &integrator),
      expected_value
  );
}

//endregion
