#ifndef CLASSICALDFT_MESH_H
#define CLASSICALDFT_MESH_H

#include "vertex.h"
#include "element.h"

#include <vector>

namespace dft_core {
namespace geometry {

enum class Direction {
  X,
  Y,
  Z
};

typedef std::reference_wrapper<Element> element_refwrap;
typedef std::vector<Element> element_vec;
typedef std::unordered_map<int, element_refwrap> element_map;

class Mesh {

 protected:
  vertex_vec vertices_raw_ = {};
  element_vec elements_raw_ = {};
  vertex_map vertices_ = {};
  element_map elements_ = {};
  std::vector<int> shape_ = {};
  std::vector<double> length_ = {};

 public:
  // region Cttors:

  Mesh() = default;

  // endregion

  // region Inspectors:

  const std::vector<int>& shape() const;
  const std::vector<double>& length() const;
  long number_vertices() const;

  // endregion

  // region Methods:

  virtual double volume() const = 0;

  // endregion

  // region Overloads:

  virtual const Vertex& operator[](const std::initializer_list<int>& idx) const = 0;

  // endregion
};

}}
#endif
