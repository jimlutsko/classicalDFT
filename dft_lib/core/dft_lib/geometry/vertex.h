#ifndef CLASSICALDFT_POINT_H
#define CLASSICALDFT_POINT_H

#include <vector>
#include <iostream>

namespace dft_core {
namespace geometry {

using x_type = double;
using vec_type = std::vector<x_type>;

const int DEFAULT_DIMENSION = 0;
const double DEFAULT_COORDINATES = 0.0;

/**
 * @brief: Vertex is a convenient class to represent a spatial point.
 */
class Vertex {
 private:
  int dimension_ = DEFAULT_DIMENSION;
  vec_type coordinates_ = vec_type();

 public:
  // region Cttors:

  /**
   * @brief: Initialises the Vertex object with dimension_ = 0 and empty coordinates_
   */
  Vertex() = default;

  /**
   * @brief: Initialises the Vertex object with the coordinates passed in the initalizer_list.
   * @param x Set of values representing the coordinates (x1,...,xN)
   */
  Vertex(const std::initializer_list<x_type>& x);

  /**
   * @brief: Initialises the Vertex object with the coordinates passed in vector x.
   * @param x: Set of values representing the coordinates (x1,...,xN)
   */
  explicit Vertex(const vec_type& x);

  /**
   * @brief: Initialises the Vertex object with the coordinates passed in vector x (move semantics).
   * @param x: Set of values representing the coordinates (x1,...,xN)
   */
  explicit Vertex(vec_type&& x);

  // endregion

  // region Inspectors:
  /**
   * Returns the size of the coordinates_ vector
   * @return Integer number of coordinates
   */
  int dimension() const;

  /**
   * Set of coordinate values representing the spatial point
   * @return Vector of coordinates
   */
  const std::vector<x_type>& coordinates() const;
  // endregion

  // region Overloads:
  const x_type& operator[](int idx) const;

  friend std::ostream& operator<<(std::ostream& os, const Vertex& vertex);
  // endregion
};

}}

#endif