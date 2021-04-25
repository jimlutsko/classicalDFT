#include "dft_lib/geometry/base/element.h"

#include <string>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <cmath>

namespace dft_core {
namespace geometry {

// region Cttors:

void Element::initialise_vertex_map(vertex_vec& v_vec, vertex_map& vertex_map)
{
  auto dimension = v_vec.front().dimension();
  for (auto k = 0; k < v_vec.size(); ++k)
  {
    auto dk = v_vec[k].dimension();
    if (dk != dimension) {
      std::string msg = "The component: " + std::to_string(k) + " of the vertex list has non-consistent dimensions";
      throw std::runtime_error(msg);
    }
    vertex_map.insert({k, std::ref(v_vec[k])});
  }
}

void Element::initialise_dimension()
{
  dimension_ = vertices_raw_.empty() ? 0 : vertices_raw_.front().dimension();
}

void Element::initialise_element()
{
  initialise_vertex_map(vertices_raw_, vertices_);
  initialise_dimension();
}

Element::Element(const std::vector<Vertex>& vertices): vertices_raw_(vertices)
{
  initialise_element();
}

Element::Element(std::vector<Vertex>&& vertices): vertices_raw_(std::move(vertices))
{
  initialise_element();
}

Element::Element(const std::initializer_list<Vertex>& vertices): vertices_raw_(vertices)
{
  initialise_element();
}

// endregion

// region Inspectors:

const vertex_vec& Element::vertices_raw() const { return vertices_raw_; }

const std::unordered_map<int, vertex_refwrap>& Element::vertices() const { return vertices_; }

int Element::number_of_vertices() const
{
  return vertices_.size();
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
  os << "\u2b91 Volume: " << element.volume() << std::endl;

  for (const auto& v : element.vertices())
  {
    os << "\u2b91 v[" << v.first << "] = " << v.second.get() << std::endl;
  }

  return os;
}

// endregion

// region SquareBox:

void SquareBox::initialise(double dx, const std::vector<double>& origin)
{
  length_ = dx;
  origin_ = origin;
  initialise_element();
}

double SquareBox::volume() const
{
  return std::pow(length_, dimension_);
}

// endregion
}}
