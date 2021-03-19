#include "dft_lib/geometry/2D/mesh.h"

#include <gtest/gtest.h>

using namespace dft_core::geometry;

// region Cttors:

TEST(geometry_2d_suqmesh, two_dim_suqmesh_default_cttor_test)
{
  auto m = two_dimensional::Lattice(0.25, {0, 0}, {1, 1});
  ASSERT_DOUBLE_EQ(4, m.number_vertices());
  ASSERT_DOUBLE_EQ(0.0625, m.volume());
}

// endregion