#include "dft_lib/geometry/3D/element.h"

#include <dft_lib/geometry/base/mesh.h>

namespace dft_core {
namespace geometry {
namespace three_dimensional {

//region SquareBox:

vertex_vec SquareBox::generate_vertices(double dx, const std::vector<double>& x0)
{
  // bottom-up anti-clock-wise:
  //       [4] - - - -[5]
  //      / |        / |
  //     /  |       /  |
  //    /   |      /   |
  //   /   [3] - -/- -[2]
  // [7] - /- - [6]   /
  //  |   /      |   /
  //  |  /       |  /
  //  | /        | /
  // [0] - - - -[1]

  auto X = static_cast<unsigned long>(Direction::X);
  auto Y = static_cast<unsigned long>(Direction::Y);
  auto Z = static_cast<unsigned long>(Direction::Z);

  auto x0_x = x0.at(X);
  auto x0_y = x0.at(Y);
  auto x0_z = x0.at(Z);

  return vertex_vec {
      Vertex(x0),
      Vertex({x0_x + dx, x0_y,      x0_z}),
      Vertex({x0_x + dx, x0_y + dx, x0_z}),
      Vertex({x0_x,      x0_y + dx, x0_z}),
      Vertex({x0_x,      x0_y + dx, x0_z + dx}),
      Vertex({x0_x + dx, x0_y + dx, x0_z + dx}),
      Vertex({x0_x + dx, x0_y,      x0_z + dx}),
      Vertex({x0_x,      x0_y,      x0_z + dx})
  };
}

SquareBox::SquareBox()
{
  origin_ = std::vector<double>{0, 0, 0};
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
  if (vertices.size() == 8)
  {
    vertices_raw_ = std::move(vertices);
    initialise_element();
  } else {
    throw std::runtime_error("2D square-box needs 4 vertices to be initialised");
  }
}

//endregion

}}}