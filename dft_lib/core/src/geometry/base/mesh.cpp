
#include <dft_lib/geometry/base/mesh.h>

#include "dft_lib/geometry/2D/mesh.h"

namespace dft_core {
namespace geometry {

const std::vector<long>& dft_core::geometry::Mesh::shape() const { return shape_; }
const std::vector<double>& dft_core::geometry::Mesh::dimensions() const { return dimensions_; }
long Mesh::number_vertices() const { return number_vertices_; }
const vertex_vec& Mesh::vertices() const { return vertices_raw_; }

std::ostream& operator<<(std::ostream& os, const Mesh& mesh)
{
  os << "Mesh object [" << std::addressof(mesh) << "]" << std::endl;
  os << "\u2b91 Volume: " << mesh.volume() << std::endl;
  os << "\u2b91 Number of vertices: " << mesh.number_vertices() << std::endl;
  os << "\u2b91 Shape: [" << mesh.shape()[0] << ", " << mesh.shape()[1] << "]" << std::endl;
  os << "\u2b91 Dimensions: [" << mesh.dimensions()[0] << ", " << mesh.dimensions()[1] << "]" << std::endl;
  return os;
}

}}