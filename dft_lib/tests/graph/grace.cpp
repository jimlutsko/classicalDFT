#include <gtest/gtest.h>

#include "dft_lib/graph/grace.h"
#include "dft_lib/exceptions/grace_exception.h"

TEST(grace_plot, min_axis_value_test)
{
  EXPECT_DOUBLE_EQ(dft_core::grace_plot::default_min_axis_value, 0.0);
}

TEST(grace_plot, max_axis_value_test)
{
  EXPECT_DOUBLE_EQ(dft_core::grace_plot::default_max_axis_value, 10.0);
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

TEST(grace_plot, default_offset_test)
{
  EXPECT_FLOAT_EQ(dft_core::grace_plot::default_offset, 0.1);
}

TEST(grace_plot, default_hspace_test)
{
  EXPECT_FLOAT_EQ(dft_core::grace_plot::default_horizontal_space, 0.15);
}

TEST(grace_plot, default_vspace_test)
{
  EXPECT_FLOAT_EQ(dft_core::grace_plot::default_vertical_space, 0.2);
}

TEST(grace, cttor_show_false_test)
{
  auto g = dft_core::grace_plot::Grace(0,0,0,false);
  ASSERT_DOUBLE_EQ(g.x_min(), 0.0);
  ASSERT_DOUBLE_EQ(g.x_max(), 10.0);
  ASSERT_DOUBLE_EQ(g.y_min(), 0.0);
  ASSERT_DOUBLE_EQ(g.y_max(), 10.0);
  ASSERT_FLOAT_EQ(g.offset(), 0.1);
  ASSERT_FLOAT_EQ(g.horizontal_space(), 0.15);
  ASSERT_FLOAT_EQ(g.vertical_space(), 0.2);
  ASSERT_FALSE(g.is_initialised());
}

TEST(grace_plot, send_command_throws_grace_not_opened_exception_test)
{
  EXPECT_THROW(
      dft_core::grace_plot::SendCommand("hello world"),
      dft_core::exception::GraceNotOpenedException
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

TEST(gace_plot_options, free_constant)
{
  ASSERT_STREQ(dft_core::grace_plot::option::FREE.c_str(), "-free");
}

TEST(gace_plot_options, nosafe_constant)
{
  ASSERT_STREQ(dft_core::grace_plot::option::NO_SAFE.c_str(), "-nosafe");
}

TEST(gace_plot_options, geometry_constant)
{
  ASSERT_STREQ(dft_core::grace_plot::option::GEOMETRY.c_str(), "-geometry");
}

TEST(grace_plot, arrange_command_ok_test)
{
  std::string actual = dft_core::grace_plot::command::Arrange(1, 2, 0.3, 0.4, 0.5);
  std::string expected = "ARRANGE(1, 2, 0.300000, 0.400000, 0.500000)";
  ASSERT_STREQ(actual.c_str(), expected.c_str());
}

TEST(grace_plot, x_min_command_ok_test)
{
  std::string actual = dft_core::grace_plot::command::SetXMinCommand(1.5);
  std::string expected = "WORLD XMIN 1.500000";
  ASSERT_STREQ(actual.c_str(), expected.c_str());
}

TEST(grace_plot, x_max_command_ok_test)
{
  std::string actual = dft_core::grace_plot::command::SetXMaxCommand(1.5);
  std::string expected = "WORLD XMAX 1.500000";
  ASSERT_STREQ(actual.c_str(), expected.c_str());
}

TEST(grace_plot, y_min_command_ok_test)
{
  std::string actual = dft_core::grace_plot::command::SetYMinCommand(1.5);
  std::string expected = "WORLD YMIN 1.500000";
  ASSERT_STREQ(actual.c_str(), expected.c_str());
}

TEST(grace_plot, y_max_command_ok_test)
{
  std::string actual = dft_core::grace_plot::command::SetYMaxCommand(1.5);
  std::string expected = "WORLD YMAX 1.500000";
  ASSERT_STREQ(actual.c_str(), expected.c_str());
}

TEST(grace_plot, start_grace_communication_negative_parameters_throws_grace_test)
{
  EXPECT_THROW(
      dft_core::grace_plot::StartGraceCommunication(-10, 10, 1000),
      dft_core::exception::GraceException
  );

  EXPECT_THROW(
      dft_core::grace_plot::StartGraceCommunication(10, -10, 1000),
      dft_core::exception::GraceException
  );

  EXPECT_THROW(
      dft_core::grace_plot::StartGraceCommunication(10, 10, -1000),
      dft_core::exception::GraceException
  );
}

TEST(grace_plot, get_number_of_rows_test)
{
  auto actual = dft_core::grace_plot::GetNumberOfRows(1);
  auto expected = 1;
  ASSERT_EQ(actual, expected);

  actual = dft_core::grace_plot::GetNumberOfRows(2);
  expected = 1;
  ASSERT_EQ(actual, expected);

  actual = dft_core::grace_plot::GetNumberOfRows(10);
  expected = 2;
  ASSERT_EQ(actual, expected);
}

TEST(grace_plot, get_number_of_cols_test)
{
  auto actual = dft_core::grace_plot::GetNumberOfColumns(1, 1);
  auto expected = 1;
  ASSERT_EQ(actual, expected);

  actual = dft_core::grace_plot::GetNumberOfColumns(2, 1);
  expected = 2;
  ASSERT_EQ(actual, expected);

  actual = dft_core::grace_plot::GetNumberOfColumns(10, 2);
  expected = 5;
  ASSERT_EQ(actual, expected);
}

TEST(grace_plot, get_numer_of_rows_neg_or_zero_throws_except_test)
{
  EXPECT_THROW(
      dft_core::grace_plot::GetNumberOfRows(-10),
      dft_core::exception::GraceException
      );
  EXPECT_THROW(
      dft_core::grace_plot::GetNumberOfRows(0),
      dft_core::exception::GraceException
  );
}

TEST(grace_plot, get_numer_of_cols_throws_except_test)
{
  EXPECT_THROW(
      dft_core::grace_plot::GetNumberOfColumns(-10, 2),
      dft_core::exception::GraceException
  );
  EXPECT_THROW(
      dft_core::grace_plot::GetNumberOfColumns(0, 2),
      dft_core::exception::GraceException
  );
  EXPECT_THROW(
      dft_core::grace_plot::GetNumberOfColumns(2, 0),
      dft_core::exception::GraceException
  );
  EXPECT_THROW(
      dft_core::grace_plot::GetNumberOfColumns(2, -10),
      dft_core::exception::GraceException
  );
}
