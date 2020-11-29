#ifndef CLASSICALDFT_ELEMENT_H
#define CLASSICALDFT_ELEMENT_H

#include "dft_lib/geometry/vertex.h"

#include <vector>
#include <unordered_map>
#include <functional>

namespace dft_core {
namespace geometry {

const std::vector<Vertex> DEFAULT_VERTICES_RAW = std::vector<Vertex>{};
typedef std::reference_wrapper<Vertex> vertex_refwrap;
typedef std::vector<Vertex> vertex_vec;

class Element {
 private:
  vertex_vec vertices_raw_ = DEFAULT_VERTICES_RAW;
  std::unordered_map<int, vertex_refwrap> vertices_ = {};

 public:
  // region Cttors:

  Element() = default;
  //Element(const std::initializer_list<Vertex>& x);

  explicit Element(const std::vector<Vertex>& vertices);
  explicit Element(std::vector<Vertex>&& vertices);

  // endregion

  // region Inspectors:

  const vertex_vec& vertices_raw() const;

  const std::unordered_map<int, vertex_refwrap>& vertices() const;

  int number_of_vertices() const;

  virtual int dimension() const { return vertices_raw().empty() ? 0 : vertices_raw().front().dimension(); };
  virtual double length() const { return 0.0; };
  virtual double volume() const { return 0.0; };
  // endregion

  // region Overloads:
  const Vertex& operator[](int idx) const;

  friend std::ostream& operator<<(std::ostream& os, const Element& element);
  // endregion
};

}}

#endif
