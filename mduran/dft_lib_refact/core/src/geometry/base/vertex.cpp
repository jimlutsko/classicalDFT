#include "dft_lib/geometry/base/vertex.h"

#include <boost/range/combine.hpp>

namespace dft_core {
namespace geometry {

// region Cttors:
Vertex::Vertex(const vec_type& x): dimension_(x.size()), coordinates_(x) {}
Vertex::Vertex(vec_type&& x): dimension_(x.size()), coordinates_(std::move(x)) {}
Vertex::Vertex(const std::initializer_list<x_type>& x): dimension_(x.size()), coordinates_(x) {}
// endregion

// region Methods:
int Vertex::dimension() const { return dimension_; }
const std::vector<x_type>& Vertex::coordinates() const { return coordinates_; }
// endregion

// region Overloads:
const x_type& Vertex::operator[](int k) const { return coordinates_.at(k); }

std::ostream& operator<<(std::ostream& os, const Vertex& vertex)
{
  os << "(dimensions = " << vertex.dimension() << "): ";
  if (vertex.dimension() > 0)
  {
    os << "(";
    for (int k = 0; k < vertex.dimension()-1; k++) { os << vertex.coordinates()[k] << ", "; }
    os << vertex.coordinates()[vertex.dimension()-1] << ")";
  }
  else {
    os << "()";
  }
  os << " [front \u27fc " << std::addressof(vertex.coordinates_.front()) << "]";
  return os;
}

Vertex operator+(const Vertex& a, const Vertex& b)
{
  if (a.dimension_ != b.dimension_) { throw std::runtime_error("The vertices you're trying to add don't have the same dimension");}

  vec_type x = vec_type();
  for (auto tup : boost::combine(a.coordinates(), b.coordinates()))
  {
    double y, z; boost::tie(y,z) = tup;
    x.push_back(y+z);
  }

  return Vertex(std::move(x));
}

Vertex operator-(const Vertex& a, const Vertex& b)
{
  if (a.dimension_ != b.dimension_) { throw std::runtime_error("The vertices you're trying to add don't have the same dimension");}

  vec_type x = vec_type();
  for (auto tup : boost::combine(a.coordinates(), b.coordinates()))
  {
    double y, z; boost::tie(y,z) = tup;
    x.push_back(y-z);
  }

  return Vertex(std::move(x));
}
// endregion
}}