#include <gtest/gtest.h>

#include "dft_lib/graph/grace.h"
#include "dft_lib/exceptions/grace_exception.h"

template<class T = double>
using Vector = std::vector<T>;

//region Default values
TEST(grace_plot, default_min_axis_value_test)
{
  EXPECT_DOUBLE_EQ(dft_core::grace_plot::default_min_axis_value, 0.0);
}

TEST(grace_plot, default_max_axis_value_test)
{
  EXPECT_DOUBLE_EQ(dft_core::grace_plot::default_max_axis_value, 10.0);
}

TEST(grace_plot, default_dataset_number_test)
{
  EXPECT_EQ(dft_core::grace_plot::default_dataset_id, 0);
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
//endregion

//region Commands
TEST(grace_plot_command, arrange_command_ok_test)
{
  std::string actual = dft_core::grace_plot::command::ArrangeCommand(1, 2, 0.3, 0.4, 0.5);
  std::string expected = "ARRANGE(1, 2, 0.300000, 0.400000, 0.500000)";
  ASSERT_STREQ(actual.c_str(), expected.c_str());
}

TEST(grace_plot_command, x_min_command_ok_test)
{
  std::string actual = dft_core::grace_plot::command::SetXMinCommand(1.5);
  std::string expected = "WORLD XMIN 1.500000";
  ASSERT_STREQ(actual.c_str(), expected.c_str());
}

TEST(grace_plot_command, x_max_command_ok_test)
{
  std::string actual = dft_core::grace_plot::command::SetXMaxCommand(1.5);
  std::string expected = "WORLD XMAX 1.500000";
  ASSERT_STREQ(actual.c_str(), expected.c_str());
}

TEST(grace_plot_command, y_min_command_ok_test)
{
  std::string actual = dft_core::grace_plot::command::SetYMinCommand(1.5);
  std::string expected = "WORLD YMIN 1.500000";
  ASSERT_STREQ(actual.c_str(), expected.c_str());
}

TEST(grace_plot_command, y_max_command_ok_test)
{
  std::string actual = dft_core::grace_plot::command::SetYMaxCommand(1.5);
  std::string expected = "WORLD YMAX 1.500000";
  ASSERT_STREQ(actual.c_str(), expected.c_str());
}

TEST(grace_plot_command, add_point_command_ok_test)
{
  std::string actual = dft_core::grace_plot::command::AddPointCommand(1, 2, 0, 0);
  std::string expected = "G0.S0 POINT 1.000000,2.000000";
  ASSERT_STREQ(actual.c_str(), expected.c_str());
}

TEST(grace_plot_command, redraw_command_ok_test)
{
  std::string actual = dft_core::grace_plot::command::RedrawCommand();
  std::string expected = "REDRAW";
  ASSERT_STREQ(actual.c_str(), expected.c_str());
}

TEST(grace_plot_command, autoticks_command_ok_test)
{
  std::string actual = dft_core::grace_plot::command::AutoTicksCommand();
  std::string expected = "AUTOTICKS";
  ASSERT_STREQ(actual.c_str(), expected.c_str());
}

TEST(grace_plot_command, autoscale_command_ok_test)
{
  std::string actual = dft_core::grace_plot::command::AutoScaleCommand();
  std::string expected = "AUTOSCALE";
  ASSERT_STREQ(actual.c_str(), expected.c_str());
}

TEST(grace_plot_command, focus_command_ok_test)
{
  std::string actual = dft_core::grace_plot::command::FocusCommand(0);
  std::string expected = "FOCUS G0";
  ASSERT_STREQ(actual.c_str(), expected.c_str());
}

TEST(grace_plot_command, kill_set_command_ok_test)
{
  std::string actual = dft_core::grace_plot::command::KillSetCommand(1, 2);
  std::string expected = "KILL G2.S1";
  ASSERT_STREQ(actual.c_str(), expected.c_str());
}

TEST(grace_plot_command, set_legend_ok_test)
{
  std::string actual = dft_core::grace_plot::command::SetLegendCommand("test", 1, 0);
  std::string expected = "G0.S1 LEGEND \"test\"";
  ASSERT_STREQ(actual.c_str(), expected.c_str());
}

TEST(grace_plot_command, set_line_color_id_ok_test)
{
  std::string actual = dft_core::grace_plot::command::SetLineColorCommand(dft_core::grace_plot::Color::RED, 1, 0);
  std::string expected = "G0.S1 LINE COLOR 2";
  ASSERT_STREQ(actual.c_str(), expected.c_str());
}

TEST(grace_plot_command, set_label_xaxis_ok_test)
{
  std::string expected = "XAXIS LABEL \"test\"";
  std::string actual = dft_core::grace_plot::command::SetAxisLabelCommand("test", dft_core::grace_plot::Axis::X);
  ASSERT_STREQ(actual.c_str(), expected.c_str());
}

TEST(grace_plot_command, set_label_yaxis_ok_test)
{
  std::string expected = "YAXIS LABEL \"test\"";
  std::string actual = dft_core::grace_plot::command::SetAxisLabelCommand("test", dft_core::grace_plot::Axis::Y);
  ASSERT_STREQ(actual.c_str(), expected.c_str());
}

TEST(grace_plot_command, set_title_ok_test)
{
  std::string expected = "TITLE \"test\"";
  std::string actual = dft_core::grace_plot::command::SetTitleCommand("test");
  ASSERT_STREQ(actual.c_str(), expected.c_str());
}

TEST(grace_plot_command, set_subtitle_ok_test)
{
  std::string expected = "SUBTITLE \"subtest\"";
  std::string actual = dft_core::grace_plot::command::SetSubtitleCommand("subtest");
  ASSERT_STREQ(actual.c_str(), expected.c_str());
}

TEST(grace_plot_command, set_ticks_ok_test)
{
  std::string expected = "XAXIS TICK MAJOR 0.010000";
  std::string actual = dft_core::grace_plot::command::SetTicksCommand(0.01, dft_core::grace_plot::Axis::X);
  ASSERT_STREQ(expected.c_str(), actual.c_str());
}

TEST(grace_plot_command, set_symbol_ok_test)
{
  std::string expected = "G0.S1 SYMBOL 2";
  std::string actual = dft_core::grace_plot::command::SetSymbolCommand(dft_core::grace_plot::Symbol::SQUARE, 1, 0);
  ASSERT_STREQ(expected.c_str(), actual.c_str());
}

TEST(grace_plot_command, set_symbol_color_ok_test)
{
  std::string expected = "G0.S1 SYMBOL COLOR 1";
  std::string actual = dft_core::grace_plot::command::SetSymbolColorCommand(dft_core::grace_plot::Color::BLACK, 1, 0);
  ASSERT_STREQ(expected.c_str(), actual.c_str());
}

TEST(grace_plot_command, set_symbol_fill_color_ok_test)
{
  std::string expected = "G0.S1 SYMBOL FILL COLOR 1";
  std::string actual = dft_core::grace_plot::command::SetSymbolColorFillCommand(dft_core::grace_plot::Color::BLACK, 1, 0);
  ASSERT_STREQ(expected.c_str(), actual.c_str());
}

TEST(grace_plot_command, set_symbol_fill_pattern_ok_test)
{
  std::string expected = "G0.S1 SYMBOL FILL PATTERN 1";
  std::string actual = dft_core::grace_plot::command::SetSymbolColorFillPatternCommand(1, 1, 0);
  ASSERT_STREQ(expected.c_str(), actual.c_str());
}

TEST(grace_plot_command, set_symbol_line_type_ok_test)
{
  std::string expected = "G0.S1 LINE LINESTYLE 1";
  std::string actual = dft_core::grace_plot::command::SetLineStyleCommand(
      dft_core::grace_plot::LineStyle::LINE, 1, 0);
  ASSERT_STREQ(expected.c_str(), actual.c_str());
}

TEST(grace_plot_command, set_symbol_size_ok_test)
{
  std::string expected = "G0.S1 SYMBOL SIZE 0.500000";
  std::string actual = dft_core::grace_plot::command::SetSymbolSizeCommand(0.5, 1, 0);
  ASSERT_STREQ(expected.c_str(), actual.c_str());
}

TEST(grace_plot_command, set_format_ok_test)
{
  std::string expected = "HARDCOPY DEVICE \"PDF\"";
  std::string actual = dft_core::grace_plot::command::SetFormatCommand(dft_core::grace_plot::ExportFormat::PDF);
  ASSERT_STREQ(expected.c_str(), actual.c_str());

  expected = "HARDCOPY DEVICE \"PNG\"";
  actual = dft_core::grace_plot::command::SetFormatCommand(dft_core::grace_plot::ExportFormat::PNG);
  ASSERT_STREQ(expected.c_str(), actual.c_str());

  expected = "HARDCOPY DEVICE \"JPG\"";
  actual = dft_core::grace_plot::command::SetFormatCommand(dft_core::grace_plot::ExportFormat::JPG);
  ASSERT_STREQ(expected.c_str(), actual.c_str());

  expected = "HARDCOPY DEVICE \"EPS1\"";
  actual = dft_core::grace_plot::command::SetFormatCommand(dft_core::grace_plot::ExportFormat::EPS1);
  ASSERT_STREQ(expected.c_str(), actual.c_str());

  expected = "HARDCOPY DEVICE \"PS\"";
  actual = dft_core::grace_plot::command::SetFormatCommand(dft_core::grace_plot::ExportFormat::PS);
  ASSERT_STREQ(expected.c_str(), actual.c_str());
}

TEST(grace_plot_command, print_to_file_test_ok)
{
  std::string expected = "PRINT TO \"test.t\"";
  std::string actual = dft_core::grace_plot::command::PrintToFileCommand("test.t");
  ASSERT_STREQ(expected.c_str(), actual.c_str());
}

TEST(grace_plot_command, print_test_ok)
{
  std::string expected = "PRINT";
  std::string actual = dft_core::grace_plot::command::PrintCommand();
  ASSERT_STREQ(expected.c_str(), actual.c_str());
}
//endregion

//region Methods
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
//endregion

//region Options
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
//endregion

//region Grace class
TEST(grace_class, cttor_show_false_test)
{
  auto g = dft_core::grace_plot::Grace(0,0,0,false);
  ASSERT_DOUBLE_EQ(g.x_min(), 0.0);
  ASSERT_DOUBLE_EQ(g.x_max(), 10.0);
  ASSERT_DOUBLE_EQ(g.y_min(), 0.0);
  ASSERT_DOUBLE_EQ(g.y_max(), 10.0);
  ASSERT_FLOAT_EQ(g.offset(), 0.1);
  ASSERT_EQ(g.last_dataset_id(), 0);
  ASSERT_EQ(g.number_of_graphs(), 0);
  ASSERT_FLOAT_EQ(g.horizontal_space(), 0.15);
  ASSERT_FLOAT_EQ(g.vertical_space(), 0.2);
  ASSERT_FALSE(g.is_initialised());
}

TEST(grace_class, getters_setter_xmin_ok_test)
{
  auto g = dft_core::grace_plot::Grace(10,10,0,false);
  double x_min_initial = g.x_min();

  auto expected_min = -10;
  g.SetXMin(expected_min);

  double x_min_new = g.x_min();
  ASSERT_NE(x_min_initial, x_min_new);
  ASSERT_DOUBLE_EQ(x_min_new, expected_min);
}

TEST(grace_class, getters_setter_xmax_ok_test)
{
  auto g = dft_core::grace_plot::Grace(10,10,0,false);
  double x_max_initial = g.x_max();

  auto expected_max = 123;
  g.SetXMax(expected_max);

  double x_max_new = g.x_max();
  ASSERT_NE(x_max_initial, x_max_new);
  ASSERT_DOUBLE_EQ(x_max_new, expected_max);
}

TEST(grace_class, getters_setter_ymin_ok_test)
{
  auto g = dft_core::grace_plot::Grace(10,10,0,false);
  double y_min_initial = g.y_min();

  auto expected_min = -10;
  g.SetYMin(expected_min);

  double y_min_new = g.y_min();
  ASSERT_NE(y_min_initial, y_min_new);
  ASSERT_DOUBLE_EQ(y_min_new, expected_min);
}

TEST(grace_class, getters_setter_ymax_ok_test)
{
  auto g = dft_core::grace_plot::Grace(10,10,0,false);
  double y_max_initial = g.y_max();

  auto expected_max = 123;
  g.SetYMax(expected_max);

  double y_max_new = g.y_max();
  ASSERT_NE(y_max_initial, y_max_new);
  ASSERT_DOUBLE_EQ(y_max_new, expected_max);
}

TEST(grace_class, getters_setter_x_limits_ok_test)
{
  auto g = dft_core::grace_plot::Grace(10,10,0,false);
  auto expected_min = -10;
  auto expected_max = 20;
  g.SetXLimits(expected_min, expected_max);

  double x_min = g.x_min();
  double x_max = g.x_max();
  ASSERT_DOUBLE_EQ(x_min, expected_min);
  ASSERT_DOUBLE_EQ(x_max, expected_max);
}

TEST(grace_class, getters_setter_x_limits_throws_excp_test)
{
  auto g = dft_core::grace_plot::Grace(10,10,0,false);
  auto expected_min = 20;
  auto expected_max = 10;
  EXPECT_THROW(
      g.SetXLimits(expected_min, expected_max),
      dft_core::exception::GraceException
  );
}

TEST(grace_class, getters_setter_y_limits_ok_test)
{
  auto g = dft_core::grace_plot::Grace(10,10,0,false);
  auto expected_min = -10;
  auto expected_max = 20;
  g.SetYLimits(expected_min, expected_max);

  double y_min = g.y_min();
  double y_max = g.y_max();
  ASSERT_DOUBLE_EQ(y_min, expected_min);
  ASSERT_DOUBLE_EQ(y_max, expected_max);
}

TEST(grace_class, getters_setter_y_limits_throws_excp_test)
{
  auto g = dft_core::grace_plot::Grace(10,10,0,false);
  auto expected_min = 20;
  auto expected_max = 10;
  EXPECT_THROW(
      g.SetYLimits(expected_min, expected_max),
      dft_core::exception::GraceException
  );
}

TEST(grace_class, getters_setter_limits_ok_test)
{
  auto g = dft_core::grace_plot::Grace(10,10,0,false);

  double expected_x_min = -5; double expected_x_max = 12;
  double expected_y_min = -10; double expected_y_max = 20;

  std::vector<double> x_limits { expected_x_min, expected_x_max };
  std::vector<double> y_limits { expected_y_min, expected_y_max };

  g.SetLimits(x_limits, y_limits);

  auto actual_x_min = g.x_min(); auto actual_x_max = g.x_max();
  auto actual_y_min = g.y_min(); auto actual_y_max = g.y_max();

  ASSERT_DOUBLE_EQ(actual_x_min, expected_x_min);
  ASSERT_DOUBLE_EQ(actual_y_min, expected_y_min);
  ASSERT_DOUBLE_EQ(actual_x_max, expected_x_max);
  ASSERT_DOUBLE_EQ(actual_y_max, expected_y_max);
}

TEST(grace_class, add_point_throws_excp_test)
{
  auto g = dft_core::grace_plot::Grace(10,10,0,false);

  EXPECT_THROW(
      g.AddPoint(1, 1, 0, 3),
      dft_core::exception::GraceException
  );
}

TEST(grace_class, redraw_throws_excp_graph_id_test)
{
  auto g = dft_core::grace_plot::Grace(10,10,0,true);

  EXPECT_THROW(
      g.Redraw(false, false, 12),
      dft_core::exception::GraceException
  );

  EXPECT_THROW(
      g.Redraw(false, false, -10),
      dft_core::exception::GraceException
  );
}

TEST(grace_class, add_dataset_throws_excp_xy_size_test)
{
  auto g = dft_core::grace_plot::Grace(10,10,0,true);

  auto x = Vector<>(10);
  auto y = Vector<>(12);

  EXPECT_THROW(
      g.AddDataset(x, y),
      dft_core::exception::GraceException
  );
}

TEST(grace_class, add_dataset_throws_excp_graph_id_test)
{
  auto g = dft_core::grace_plot::Grace(10,10,0,true);

  auto x = Vector<>(10);
  auto y = Vector<>(10);

  EXPECT_THROW(
      g.AddDataset(x, y, 12),
      dft_core::exception::GraceException
  );
}

TEST(grace_class, add_dataset_works_ok)
{
  auto g = dft_core::grace_plot::Grace(10,10);

  auto x = Vector<>(10);
  auto y = Vector<>(10);

  int dataset_id = g.AddDataset(x, y);
  ASSERT_EQ(g.last_dataset_id(), dataset_id);
}

TEST(grace_class, replace_dataset_throws_excp_dataset_id_test)
{
  auto g = dft_core::grace_plot::Grace(10,10,1,true);

  auto x = Vector<>(10);
  auto y = Vector<>(10);

  g.AddPoint(0, 0);
  EXPECT_THROW(
      g.ReplaceDataset(x, y, 1),
      dft_core::exception::GraceException
  );
}

TEST(grace_class, replace_dataset_throws_excp_graph_id_test)
{
  auto g = dft_core::grace_plot::Grace(10,10,1,true);

  auto x = Vector<>(10);
  auto y = Vector<>(10);

  g.AddPoint(0, 0);
  EXPECT_THROW(
      g.ReplaceDataset(x, y, 0, 2),
      dft_core::exception::GraceException
  );
}

TEST(grace_class, replace_dataset_works_ok)
{
  auto g = dft_core::grace_plot::Grace(10,10,1,true);

  auto x = Vector<>(10);
  auto y = Vector<>(10);

  g.AddPoint(0, 0);
  g.ReplaceDataset(x, y, 0, 0);

  EXPECT_EQ(0, g.last_dataset_id());
}

TEST(grace_class, delete_dataset_throws_excp_dataset_id)
{
  auto g = dft_core::grace_plot::Grace(10,10,1,false);

  auto x = Vector<>(10);
  auto y = Vector<>(10);

  g.AddPoint(0, 0);
  EXPECT_THROW(
      g.DeleteDataset(1, 0),
      dft_core::exception::GraceException
  );
}

TEST(grace_class, delete_dataset_throws_excp_graph_id)
{
  auto g = dft_core::grace_plot::Grace(10,10,1,false);

  auto x = Vector<>(10);
  auto y = Vector<>(10);

  g.AddPoint(0, 0);
  EXPECT_THROW(
      g.DeleteDataset(0, 1),
      dft_core::exception::GraceException
  );
}

TEST(grace_class, set_legend_throws_excp_dataset_id)
{
  auto g = dft_core::grace_plot::Grace(10,10,1,false);
  g.AddPoint(0, 0);

  EXPECT_THROW(
      g.SetLegend("test", 1, 0),
      dft_core::exception::GraceException
  );
}

TEST(grace_class, set_legend_throws_excp_graph_id)
{
  auto g = dft_core::grace_plot::Grace(10,10,1,false);
  g.AddPoint(0, 0);

  EXPECT_THROW(
      g.SetLegend("test", 0, 1),
      dft_core::exception::GraceException
  );
}

TEST(grace_class, set_color_throws_excp_dataset_id)
{
  auto g = dft_core::grace_plot::Grace(10,10,1,false);
  g.AddPoint(0, 0);

  EXPECT_THROW(
      g.SetColor(dft_core::grace_plot::Color::BLUE, 1, 0),
      dft_core::exception::GraceException
  );
}

TEST(grace_class, set_color_throws_excp_graph_id)
{
  auto g = dft_core::grace_plot::Grace(10,10,1,false);
  g.AddPoint(0, 0);

  EXPECT_THROW(
      g.SetColor(dft_core::grace_plot::Color::RED, 0, 1),
      dft_core::exception::GraceException
  );
}

TEST(grace_class, set_lable_throws_excp_graph_id)
{
  auto g = dft_core::grace_plot::Grace(10,10,1,false);
  g.AddPoint(0, 0);

  EXPECT_THROW(
      g.SetLabel("test", dft_core::grace_plot::Axis::X, 12),
      dft_core::exception::GraceException
  );

  EXPECT_THROW(
      g.SetLabel("test", dft_core::grace_plot::Axis::Y, 12),
      dft_core::exception::GraceException
  );
}

TEST(grace_class, set_ticks_throws_excp_negative_sep)
{
  auto g = dft_core::grace_plot::Grace(10,10,1,false);
  g.AddPoint(0, 0);

  EXPECT_THROW(
      g.SetTicks(-0.01, 0.02),
      dft_core::exception::GraceException
  );

  EXPECT_THROW(
      g.SetTicks(0.01, -0.02),
      dft_core::exception::GraceException
  );

  EXPECT_THROW(
      g.SetTicks(0.01, -0.02, 2),
      dft_core::exception::GraceException
  );
}
//endregion