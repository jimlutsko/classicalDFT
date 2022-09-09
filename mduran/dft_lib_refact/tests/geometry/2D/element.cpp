#include "dft_lib/geometry/2D/element.h"

#include <gtest/gtest.h>

using namespace dft_core::geometry;

// region Cttors:
TEST(geometry_2d_square_box, two_dim_squarebox_default_cttor_test)
{
  auto e = two_dimensional::SquareBox();

  ASSERT_EQ(4, e.number_of_vertices());
  ASSERT_EQ(2, e.dimension());
  ASSERT_DOUBLE_EQ(1, e.volume());
}

TEST(geometry_2d_square_box, two_dim_squarebox_indexer)
{
  auto e = two_dimensional::SquareBox();

  ASSERT_EQ(0, e[0].coordinates()[0]);
  ASSERT_EQ(0, e[0].coordinates()[0]);
}
// endregion