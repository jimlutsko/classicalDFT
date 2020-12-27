#ifndef CLASSICALDFT_ELEMENT_H
#define CLASSICALDFT_ELEMENT_H

#include "vertex.h"

#include <functional>
#include <unordered_map>
#include <vector>

namespace dft_core {
namespace geometry {

const std::vector<Vertex> DEFAULT_VERTICES_RAW = std::vector<Vertex>{};

// region Type definitions:

typedef std::reference_wrapper<Vertex> vertex_refwrap;
typedef std::vector<Vertex> vertex_vec;
typedef std::unordered_map<int, vertex_refwrap> vertex_map;

// endregion

/**
 * @brief Element is a convenient class to represent a spatial element. This class serves as base
 *         for more specific geometric elements, e.g. square-box, trapezoid, etc.
 */
class Element {

 protected:

  /**
   * @brief The spatial dimension of element
   */
  int dimension_ = DEFAULT_DIMENSION;

  /**
   * @brief The list of vertices (`std::vector<Vertex>`)
   */
  vertex_vec vertices_raw_ = DEFAULT_VERTICES_RAW;

  /**
   * @brief Map `int -> Vertex` which allows for index access to the vertices
   */
  vertex_map vertices_ = {};

  /**
   * @brief Initialisation of the `dimension_` property
   */
  void initialise_dimension();

  /**
   * @brief  Initialisation of the `int -> Vertex` mapping
   * @param v_vec Vector of vertices to be mapped
   * @param vertex_map Unordered-map which assigns an `int` to each Vertex
   */
  void initialise_vertex_map(vertex_vec& v_vec, vertex_map& vertex_map);

  /**
   * @brief Initialiser for the various constructors
   */
  void initialise_element();

 public:
  // region Cttors:

  /**
   * Default constructor
   */
  Element() = default;

  /**
   * Element constructor when given an arbitrarily large list of vertices
   * @param vertices Set of vertices in the form of a `Initializer_list<Vertex>`
   */
  Element(const std::initializer_list<Vertex>& vertices);

  /**
   * Element constructor when given a `std::vector` of vertices
   * @param vertices Set of vertices stored in a `std::vector<Vertex>` object
   */
  explicit Element(const std::vector<Vertex>& vertices);

  /**
   * Element constructor when given a `std::vector` of vertices (move semantics)
   * @param vertices Set of vertices stored in a `std::vector<Vertex>` object
   */
  explicit Element(std::vector<Vertex>&& vertices);

  // endregion

  // region Inspectors:
  /**
   * Gets the list of vertices stored in a std::vector object
   * @return std::vector<Vertex>
   */
  const vertex_vec& vertices_raw() const;

  /**
   * Gets the unordered map which indexes the vertices
   * @return std::unordered_map
   */
  const std::unordered_map<int, vertex_refwrap>& vertices() const;

  /**
   * Gets the number of vertices being part of the element
   * @return Integer number of vertices
   */
  int number_of_vertices() const;

  /**
   * Gets the dimension of the element
   * @return Integer dimension
   */
  virtual int dimension() const { return dimension_; };

  /**
   * Prototype function for when we need to measure the length in a given direction
   * @return Length in the requested direction
   */
  //virtual double length() const { return 0.0; };

  /**
   * Prototype function to compute the element volume
   * @return Total volume of the element
   */
  virtual double volume() const { return 0.0; };

  // endregion

  // region Overloads:
  /**
   * Indexer to get a given Vertex object via the inner index.
   * @param idx
   * @return Vertex object
   */
  const Vertex& operator[](int idx) const;

  /**
   * Overload of the operator `<<` to use with sdt::cout
   * @param os
   * @param element
   * @return
   */
  friend std::ostream& operator<<(std::ostream& os, const Element& element);
  // endregion
};

// region Abstract SquareBox:

const double DEFAULT_SQUAREBOX_LENGTH = 1;

/**
 * SquareBox is the abstract template for the creation of square-box elements in
 * 2D, 3D (and higher dimensions)
 */
class SquareBox : public Element {
 protected:
  /**
   * Width/height/length of the box
   */
  double length_ = DEFAULT_SQUAREBOX_LENGTH;

  /**
   * Element's origin coordinates
   */
  std::vector<double> origin_ = {};

  /**
   * Initialiser of the element's attributes which are general for all the possible dimensions
   * @param dx The length/width/height of the box
   * @param origin The origin's coordinates to build the box
   */
  void initialise(double dx, const std::vector<double>& origin);

 public:
  /**
   * Default constructor
   */
  SquareBox() = default;

  /**
   * Gets the volume of the element: pow(length_, dimension_)
   * @return The element's volume
   */
  double volume() const final;
};

// endregion

}}

#endif
