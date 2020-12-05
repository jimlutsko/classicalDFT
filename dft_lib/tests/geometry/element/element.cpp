#include "dft_lib/geometry/element.h"

#include <gtest/gtest.h>

#include <boost/range/combine.hpp>

// region Cttors:

using namespace dft_core::geometry;

TEST(geometry_element, element_default_cttor_test)
{
  auto e = Element();

  ASSERT_EQ(0, e.number_of_vertices());
  ASSERT_EQ(0, e.dimension());
  ASSERT_EQ(0, e.volume());
}

TEST(geometry_element, element_initializer_list_cttor_test)
{
  auto e = Element({{1, 2, 3},{4, 5, 6},{7, 8, 9}});

  ASSERT_EQ(3, e.number_of_vertices());

  int e_idx = 0;
  for (int k = 0; k < 9; k++)
  {
    e_idx += ((k!=0) && (k%3 == 0)) ? 1 : 0;
    ASSERT_EQ(k+1, e[e_idx][k%3]);
  }
}