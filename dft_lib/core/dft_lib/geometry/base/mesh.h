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

class Mesh {

 protected:
  std::vector<long> shape_ = {};
  std::vector<double> dimensions_;
  vertex_vec vertices_raw_ = {};
  vertex_map vertices_ = {};

  long number_vertices_ = 0;

 public:
  // region Cttors:

  Mesh() = default;

  // endregion

  // region Inspectors:

  const std::vector<long>& shape() const;
  const std::vector<double>& dimensions() const;
  long number_vertices() const;

  const vertex_vec& vertices() const;

  // endregion

  // region Methods:

  virtual double volume() const = 0;

  // endregion

  // region Overloads:

  virtual const Vertex& operator[](const std::vector<long>& idx) const = 0;

  friend std::ostream& operator<<(std::ostream& os, const Mesh& mesh);

  // endregion
};

}}
#endif
