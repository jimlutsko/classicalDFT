#include "dft_lib/geometry/base/mesh.h"

const std::vector<int>& dft_core::geometry::Mesh::shape() const { return shape_; }
const std::vector<double>& dft_core::geometry::Mesh::length() const { return length_; }
long dft_core::geometry::Mesh::number_vertices() const { return vertices_raw_.size(); }
