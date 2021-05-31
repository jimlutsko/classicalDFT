#ifndef CLASSICALDFT_3D_ELEMENT_H
#define CLASSICALDFT_3D_ELEMENT_H

#include "dft_lib/geometry/base/element.h"

namespace dft_core {
namespace geometry {
namespace three_dimensional {

/**
 * two_dimensional::SquareBox is the particular implementation of the SquareBox
 * in 3D.
 */
class SquareBox : public dft_core::geometry::SquareBox {
 private:
  /**
   * Private methods which generates the vertices which comprise the element
   * @param dx Length of the square box
   * @param x0 The coordinates of the first vertex of the element
   * @return Array (vector) of vertices
   */
  vertex_vec generate_vertices(double dx, const std::vector<double>& x0);

 public:
  /**
   * Default constructor of the class where basic initialisation is carried out
   */
  SquareBox();

  /**
   * Specific constructor which builds all the vertices comprising the element given the
   * length/width/height of the box, and its origin's coordinates.
   * @param length The length/width/height of the box
   * @param origin The element's origin coordinates
   */
  explicit SquareBox(double length, const std::vector<double>& origin);

  /**
   * Specific constructor which builds all the vertices based on move semantics
   * @param vertices Array of vertices to be moved to the inner attribute vertices_raw_
   */
  explicit SquareBox(std::vector<Vertex>&& vertices);
};

}}}

#endif
