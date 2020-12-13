#ifndef CLASSICALDFT_ELEMENT_H
#define CLASSICALDFT_ELEMENT_H

#include <functional>
#include <unordered_map>
#include <vector>

#include "vertex.h"

namespace dft_core {
namespace geometry {

const std::vector<Vertex> DEFAULT_VERTICES_RAW = std::vector<Vertex>{};

typedef std::reference_wrapper<Vertex> vertex_refwrap;
typedef std::vector<Vertex> vertex_vec;
typedef std::unordered_map<int, vertex_refwrap> vertex_map;

class Element {
 protected:
  int dimension_ = DEFAULT_DIMENSION;
  vertex_vec vertices_raw_ = DEFAULT_VERTICES_RAW;
  vertex_map vertices_ = {};

  void _initialise_dimension();
  void _initialise_vertex_map(vertex_vec& v_vec, vertex_map& vertex_map);
  void _initialise_element();

 public:
  // region Cttors:

  Element() = default;
  Element(const std::initializer_list<Vertex>& vertices);

  explicit Element(const std::vector<Vertex>& vertices);
  explicit Element(std::vector<Vertex>&& vertices);

  // endregion

  // region Inspectors:

  const vertex_vec& vertices_raw() const;

  const std::unordered_map<int, vertex_refwrap>& vertices() const;

  int number_of_vertices() const;

  virtual int dimension() const { return dimension_; };
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
