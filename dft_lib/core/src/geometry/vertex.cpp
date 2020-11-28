#include "dft_lib/geometry/vertex.h"

namespace dft_core {
namespace geometry {

Vertex::Vertex(const vec_type& x): dimension_(x.size()), coordinates_(x) {}
Vertex::Vertex(vec_type&& x): dimension_(x.size()), coordinates_(std::move(x)) {}
Vertex::Vertex(const std::initializer_list<x_type>& x): dimension_(x.size()), coordinates_(x) {}

int Vertex::dimension() const { return dimension_; }
const std::vector<x_type>& Vertex::coordinates() const { return coordinates_; }

}}