#include "dft_lib/geometry/3D/mesh.h"
#include "dft_lib/geometry/3D/element.h"
#include "dft_lib/graph/grace.h"

#include <numeric>

namespace dft_core {
namespace geometry {
namespace three_dimensional {

void SUQMesh::initialise(double dx)
{
  //region Vertices init:
  auto vertex_index = 0;
  auto element_index = 0;

  auto x = origin_[0]; auto y = origin_[1]; auto z = origin_[2];
  vertices_raw_ = vertex_vec(number_vertices_);
  elements_raw_ = sqbox_vec(number_elements_);

  for (auto i_idx = 0; i_idx < shape_[0]; i_idx++) {
    for (auto j_idx = 0; j_idx < shape_[1]; j_idx++) {
      for (auto k_idx = 0; k_idx < shape_[2]; k_idx++) {

        vertices_raw_[vertex_index] = Vertex({x, y, z});
        vertices_.insert({vertex_index, std::ref(vertices_raw_[vertex_index])});

        if ((i_idx < idx_max_[0]) && (j_idx < idx_max_[1]) && (k_idx < idx_max_[2])) {
          elements_raw_[element_index] = SquareBox(dx, {x, y, z});
          elements_.insert({element_index, std::ref(elements_raw_[element_index])});
          element_index += 1;
        }

        vertex_index += 1; z += dx;
      }
      z = 0; y += dx;
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
  throw std::runtime_error("Not implemented yet!");
}

const std::vector<SquareBox>& SUQMesh::elements() const {return elements_raw_;}

double SUQMesh::element_volume() const { return elements().front().volume(); }

}}}