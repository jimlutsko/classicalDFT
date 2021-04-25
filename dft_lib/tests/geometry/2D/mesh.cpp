#include "dft_lib/geometry/2D/mesh.h"

#include <gtest/gtest.h>

using namespace dft_core::geometry;

// region Cttors:

TEST(geometry_2d_suqmesh, two_dim_suqmesh_default_cttor_test)
{
  auto origin = std::vector<double>{0, 0};
  auto dimensions = std::vector<double>{1, 1};
  auto m = two_dimensional::Lattice(0.5, dimensions, origin);
  ASSERT_DOUBLE_EQ(9, m.number_vertices());
  ASSERT_DOUBLE_EQ(1, m.volume());
}

// endregion