#include "dft_lib/geometry/2D/mesh.h"
#include "dft_lib/geometry/2D/element.h"
#include "dft_lib/graph/grace.h"

#include <numeric>

namespace dft_core {
namespace geometry {
namespace two_dimensional {

void SUQMesh::initialise(double dx)
{
  //region Vertices init:
  auto vertex_index = 0;
  auto element_index = 0;

  auto x = origin_[0]; auto y = origin_[1];
  vertices_raw_ = vertex_vec(number_vertices_);
  elements_raw_ = sqbox_vec(number_elements_);

  for (auto i_idx = 0; i_idx < shape_[0]; i_idx++)
  {
    for (auto j_idx = 0; j_idx < shape_[1]; j_idx++)
    {
      vertices_raw_[vertex_index] = Vertex({x, y});
      vertices_.insert({vertex_index, std::ref(vertices_raw_[vertex_index])});

      if ((i_idx < idx_max_[0]) && (j_idx < idx_max_[1]))
      {
        elements_raw_[element_index] = SquareBox(dx, {x, y});
        elements_.insert({element_index, std::ref(elements_raw_[element_index])});
        element_index += 1;
      }

      vertex_index += 1;
      y += dx;
    }
    y = 0; x += dx;
  }
  //endregion
}


SUQMesh::SUQMesh(double dx, std::vector<double>& dimensions, std::vector<double>& origin)
    : dft_core::geometry::SUQMesh(dx, dimensions, origin)
{
  this->initialise(dx);
}

void SUQMesh::plot() const
{
  auto g = dft_core::grace_plot::Grace();
  for (const auto& v : vertices_raw_)
  { g.AddPoint(v.coordinates()[0], v.coordinates()[1]); }

  auto dx = std::vector<double>{
    0.1 * dimensions_[0], 0.1 * dimensions_[1]
  };

  g.SetLimits(
      {origin_[0]-dx[0], (dimensions_[0]+origin_[0]) + dx[0]},
      {origin_[1]-dx[1], (dimensions_[1]+origin_[1]) + dx[1]}
      );

  g.SetLineType(dft_core::grace_plot::LineStyle::NO_LINE, 0);
  g.SetSymbol(dft_core::grace_plot::Symbol::SQUARE, 0);
  g.SetSymbolColor(dft_core::grace_plot::Color::BLUE, 0);
  g.SetSymbolFill(dft_core::grace_plot::Color::DARKGREEN, 0);

  g.RedrawAndWait();
}

const std::vector<SquareBox>& SUQMesh::elements() const {return elements_raw_;}

double SUQMesh::element_volume() const { return elements().front().volume(); }

}}}