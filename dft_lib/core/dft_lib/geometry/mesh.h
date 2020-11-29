#ifndef CLASSICALDFT_MESH_H
#define CLASSICALDFT_MESH_H

#include "dft_lib/geometry/vertex.h"

#include <vector>

namespace dft_core {
namespace geometry {

const std::vector<Vertex> DEFAULT_VERTICES = std::vector<Vertex>{};

class Mesh {

 private:
  std::vector<Vertex> vertices_ = DEFAULT_VERTICES;

 public:
  Mesh() = default;

};


class Lattice: Mesh {
 private:
  //

 public:
  Lattice() = default;
};

}}
#endif
