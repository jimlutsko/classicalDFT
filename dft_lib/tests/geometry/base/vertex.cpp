#include "dft_lib/geometry/base/vertex.h"

#include <gtest/gtest.h>

#include <boost/range/combine.hpp>

// region Cttors:

TEST(geometry_vertex, point_default_cttor_test)
{
  auto p = dft_core::geometry::Vertex();

  auto actual_dim = p.dimension();
  auto expected_dim = 0;

  std::cout << expected_dim << std::endl;

  ASSERT_EQ(expected_dim, actual_dim);

  auto vec_actual = p.coordinates();
  auto vec_expected = std::vector<double>{};
  EXPECT_EQ(vec_actual.size(), vec_expected.size());
}

TEST(geometry_vertex, point_const_vec_cttor_test)
{
  auto x = std::vector<double>{0, 1, 2};
  auto p = dft_core::geometry::Vertex(x);

  auto actual_dim = p.dimension();
  auto expected_dim = x.size();
  std::cout << expected_dim << std::endl;
  ASSERT_EQ(expected_dim, actual_dim);

  for (auto tup : boost::combine(x, p.coordinates()))
  {
    double z, y; boost::tie(z,y) = tup;
    EXPECT_DOUBLE_EQ(z, y);
  }
}

TEST(geometry_vertex, point_move_vec_cttor_test)
{
  auto x = std::vector<double>{0, 1, 2};
  auto expected_dim = x.size();

  auto p = dft_core::geometry::Vertex(std::move(x));
  auto actual_dim = p.dimension();

  ASSERT_EQ(expected_dim, actual_dim);
  ASSERT_EQ(0, x.size());

  for (auto tup : boost::combine(std::vector<double>{0, 1, 2}, p.coordinates()))
  {
    double y, z; boost::tie(y,z) = tup;
    EXPECT_NEAR(z, y, 1E-10);
  }
}

TEST(geometry_vertex, point_initializer_list_cttor_test)
{
  auto x = std::vector<double>{0, 1, 2};
  auto expected_dim = x.size();

  // Implicit conversion of the initializer_list (this is why using explicit in this case
  // turns out being a bad idea, as we'd not be able to exploit this advantage.
  dft_core::geometry::Vertex p = {0,1,2};
  auto actual_dim = p.dimension();

  ASSERT_EQ(expected_dim, actual_dim);

  for (auto tup : boost::combine(std::vector<double>{0, 1, 2}, p.coordinates()))
  {
    double y, z; boost::tie(y,z) = tup;
    EXPECT_NEAR(z, y, 1E-10);
  }
}

TEST(geometry_vertex, vertex_indexer_works_ok)
{
  dft_core::geometry::Vertex p = {0.0, 1.0, 2.0};
  ASSERT_DOUBLE_EQ(0.0, p[0]);
  ASSERT_DOUBLE_EQ(1.0, p[1]);
  ASSERT_DOUBLE_EQ(2.0, p[2]);
}

TEST(geometry_vertex, vertex_sum_overload_ok)
{
  dft_core::geometry::Vertex p1 = {0.0, 1.0, 2.0};
  dft_core::geometry::Vertex p2 = {2.0, 1.0, 0.0};
  auto p_sum = p1 + p2;
  ASSERT_DOUBLE_EQ(2.0 , p_sum[0]);
  ASSERT_DOUBLE_EQ(2.0 , p_sum[1]);
  ASSERT_DOUBLE_EQ(2.0 , p_sum[2]);
}

TEST(geometry_vertex, vertex_sum_overload_throws_error)
{
  dft_core::geometry::Vertex p1 = {0.0, 1.0, 2.0};
  dft_core::geometry::Vertex p2 = {2.0, 1.0};
  EXPECT_THROW(
      p1 + p2,
      std::runtime_error
  );
}

TEST(geometry_vertex, vertex_sub_overload_ok)
{
  dft_core::geometry::Vertex p1 = {0.0, 1.0, 2.0};
  dft_core::geometry::Vertex p2 = {2.0, 1.0, 0.0};
  auto p_sub = p2-p1;
  ASSERT_DOUBLE_EQ(2.0 , p_sub[0]);
  ASSERT_DOUBLE_EQ(0.0 , p_sub[1]);
  ASSERT_DOUBLE_EQ(-2.0 , p_sub[2]);
}

TEST(geometry_vertex, vertex_sub_overload_throws_error)
{
  dft_core::geometry::Vertex p1 = {0.0, 1.0, 2.0};
  dft_core::geometry::Vertex p2 = {2.0, 1.0};
  EXPECT_THROW(
      p1 - p2,
      std::runtime_error
  );
}
// endregion
