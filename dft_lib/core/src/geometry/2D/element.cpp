#include "dft_lib/geometry/2D/element.h"

namespace dft_core {
namespace geometry {
namespace two_dimensional {

static vertex_vec _generate_vertices(double dx, const std::vector<double>& x0)
{
  // bottom-up anti-clock-wise:
  return vertex_vec {
      Vertex(x0),
      Vertex({x0.at(0)+ dx, x0.at(1)}),
      Vertex({x0.at(0)+ dx, x0.at(1)+ dx}),
      Vertex({x0.at(0), x0.at(1)+ dx}),
  };
}

SquareBox::SquareBox()
{
  dimension_ = DEFAULT_DIMENSION;
  auto origin = DEFAULT_ORIGIN;
  vertices_raw_ = _generate_vertices(DEFAULT_LENGTH, origin);
  _initialise_element();
}

SquareBox::SquareBox(double length, const std::vector<double>& origin)
{
  dimension_ = DEFAULT_DIMENSION;
  vertices_raw_ = _generate_vertices(length, origin);
  _initialise_element();
}

}}}