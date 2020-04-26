#include <gtest/gtest.h>

#include "dft_lib/graph/grace.h"
#include "dft_lib/exceptions/grace_exception.h"

TEST(grace_plot, min_axis_value_test)

{
  EXPECT_DOUBLE_EQ(dft_core::grace_plot::min_axis_value, 0.0);
}

TEST(grace_plot, max_axis_value_test)
{
  EXPECT_DOUBLE_EQ(dft_core::grace_plot::max_axis_value, 10.0);
}

TEST(grace_plot, default_dataset_number_test)
{
  EXPECT_EQ(dft_core::grace_plot::default_dataset_number, 0);
}

TEST(grace_plot, default_x_y_size_test)
{
  EXPECT_EQ(dft_core::grace_plot::default_x_size, 800);
  EXPECT_EQ(dft_core::grace_plot::default_y_size, 600);
}

TEST(grace_plot, default_number_of_graphs_test)
{
  EXPECT_EQ(dft_core::grace_plot::default_number_of_graphs, 1);
}

TEST(grace, cttor_show_false_test)
{
  auto g = dft_core::grace_plot::Grace(0,0,0,false);
  EXPECT_EQ(g.x_min(), 0.0);
  EXPECT_EQ(g.x_max(), 10.0);
  EXPECT_EQ(g.y_min(), 0.0);
  EXPECT_EQ(g.y_max(), 10.0);
  ASSERT_FALSE(g.is_initialised());
}

TEST(grace_plot, send_command_throws_grace_not_opened_exception_test)
{
  EXPECT_THROW(
      dft_core::grace_plot::SendCommand("hello world"),
      dft_core::grace_plot::GraceNotOpenedException
      );
}

TEST(grace_plot, error_parsing_function_prints_ok_test)
{
  testing::internal::CaptureStdout();
  dft_core::grace_plot::ErrorParsingFunction("test");
  std::string output = testing::internal::GetCapturedStdout();
  auto expected_str = std::string("Grace library message: \"test\"\n");
  ASSERT_STREQ(output.c_str(), expected_str.c_str());
}


TEST(grace_plot, start_grace_communication_negative_parameters_throws_grace_test)
{
  EXPECT_THROW(
      dft_core::grace_plot::StartGraceCommunication(-10, 10, 1000),
      dft_core::grace_plot::GraceException
  );

  EXPECT_THROW(
      dft_core::grace_plot::StartGraceCommunication(10, -10, 1000),
      dft_core::grace_plot::GraceException
  );

  EXPECT_THROW(
      dft_core::grace_plot::StartGraceCommunication(10, 10, -1000),
      dft_core::grace_plot::GraceException
  );
}

