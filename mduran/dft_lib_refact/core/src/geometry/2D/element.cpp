#include "dft_lib/geometry/2D/element.h"

#include <dft_lib/geometry/base/mesh.h>

namespace dft_core {
namespace geometry {
namespace two_dimensional {

vertex_vec SquareBox::generate_vertices(double dx, const std::vector<double>& x0)
{
  // bottom-up anti-clock-wise:
  // [3] - - - - [2]
  //  |           |
  //  |           |
  //  |           |
  // [0] - - - - [1]

  auto X = static_cast<unsigned long>(Direction::X);
  auto Y = static_cast<unsigned long>(Direction::Y);

  return vertex_vec {
      Vertex(x0),
      Vertex({x0.at(X) + dx, x0.at(Y)}),
      Vertex({x0.at(X) + dx, x0.at(Y) + dx}),
      Vertex({x0.at(X), x0.at(Y)+ dx}),
  };
}

SquareBox::SquareBox()
{
  origin_ = std::vector<double>{0, 0};
  length_ = DEFAULT_SQUAREBOX_LENGTH;
  vertices_raw_ = generate_vertices(length_, origin_);
  initialise(length_, origin_);
}

SquareBox::SquareBox(double length, const std::vector<double>& origin)
{
  vertices_raw_ = generate_vertices(length, origin);
  initialise(length, origin);
}

SquareBox::SquareBox(vertex_vec && vertices)
{
  if (vertices.size() == 4)
  {
    vertices_raw_ = std::move(vertices);
    initialise_element();
  } else {
    throw std::runtime_error("2D square-box needs 4 vertices to be initialised");
  }
}

}}}