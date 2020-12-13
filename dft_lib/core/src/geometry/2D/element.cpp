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

void SquareBox::_initialise(double dx, const std::vector<double>& origin)
{
  vertices_raw_ = _generate_vertices(dx, origin);
  length_ = dx;
  _initialise_element();
}

SquareBox::SquareBox()
{
  _initialise(DEFAULT_BOX_LENGTH, DEFAULT_2D_BOX_ORIGIN);
}

SquareBox::SquareBox(double length, const std::vector<double>& origin)
{
  _initialise(length, origin);
}

SquareBox::SquareBox(vertex_vec && vertices)
{
  if (vertices.size() == 4)
  {
    vertices_raw_ = std::move(vertices);
    _initialise_element();
  } else {
    throw std::runtime_error("2D square-box needs 4 vertices to be initialised");
  }
}

double SquareBox::volume() const {
  return length_ * length_;
}

}}}