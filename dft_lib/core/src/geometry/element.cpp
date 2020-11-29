#include "dft_lib/geometry/element.h"

#include <iostream>
#include <memory>
#include <iterator>
#include <stdexcept>

namespace dft_core {
namespace geometry {

// region Cttors:

Element::Element(const std::vector<Vertex>& vertices): vertices_raw_(vertices) {
  for (auto k = 0; k < vertices_raw_.size(); ++k)
  {
    vertices_.insert({k, std::ref(vertices_raw_[k])});
  }
}

Element::Element(std::vector<Vertex>&& vertices): vertices_raw_(std::move(vertices))
{
  for (auto k = 0; k < vertices_raw_.size(); ++k)
  {
    vertices_.insert({k, std::ref(vertices_raw_[k])});
  }
}

// endregion

// region Inspectors:

const vertex_vec& Element::vertices_raw() const { return vertices_raw_; }

const std::unordered_map<int, vertex_refwrap>& Element::vertices() const { return vertices_; }

int Element::number_of_vertices() const
{
  return vertices_raw().size();
}
// endregion

// region Overloads:

const Vertex& Element::operator[](int idx) const
{
  if (vertices_.find(idx) != vertices_.end())
  {
    return vertices_.at(idx).get();
  }
  else {
    throw std::out_of_range("The index {" + std::to_string(idx) + "} is out of range!");
  }
}

std::ostream& operator<<(std::ostream& os, const Element& element)
{
  os << "Element [" << std::addressof(element) << "]" << std::endl;
  os << "\u2b91 Number of vertices: " << element.number_of_vertices() << std::endl;

  for (const auto& v : element.vertices())
  {
    os << "\u2b91 v[" << v.first << "] = " << v.second.get() << std::endl;
  }

  return os;
}

// endregion

}}
