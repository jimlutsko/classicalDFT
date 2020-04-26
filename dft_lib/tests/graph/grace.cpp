#include "dft_lib/graph/grace.h"
#include <gtest/gtest.h>

TEST(graceTest, testMinAxisValue)
{
  EXPECT_DOUBLE_EQ(dft_core::grace_plot::min_axis_value, 0.0);
}

TEST(graceTest, testMaxAxisValue)
{
  EXPECT_DOUBLE_EQ(dft_core::grace_plot::max_axis_value, 10.0);
}

TEST(graceTest, testDefaultDatasetNumber)
{
  EXPECT_EQ(dft_core::grace_plot::default_dataset_number, 0);
}

TEST(graceTest, testDefaultSize)
{
  EXPECT_EQ(dft_core::grace_plot::default_x_size, 800);
  EXPECT_EQ(dft_core::grace_plot::default_y_size, 600);
}

TEST(graceTest, testDefaultNumberOfGraphs)
{
  EXPECT_EQ(dft_core::grace_plot::default_number_of_graphs, 1);
}

TEST(graceTest, testCttorShowFalse)
{
  auto g = dft_core::grace_plot::Grace(0,0,0,false);
  EXPECT_EQ(g.get_x_min(), 0.0);
  EXPECT_EQ(g.get_x_max(), 10.0);
  EXPECT_EQ(g.get_y_min(), 0.0);
  EXPECT_EQ(g.get_y_max(), 10.0);
  ASSERT_FALSE(g.is_initialised());
}
