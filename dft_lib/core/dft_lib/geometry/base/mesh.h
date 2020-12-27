#ifndef CLASSICALDFT_MESH_H
#define CLASSICALDFT_MESH_H

#include "vertex.h"
#include "element.h"

#include <vector>

namespace dft_core {
namespace geometry {

enum class Direction {
  X = 0,
  Y = 1,
  Z = 2
};

enum class Plane {
  XY,
  XZ,
  YZ,
};

// region Pure-Abstract Mesh:

/**
 * Convenient abstraction of what the basic ingredients of a general mesh must have
 */
class Mesh {

 protected:
  /**
   * The shape of the mesh which will record the number of vertices in each direction
   */
  std::vector<long> shape_ = {};

  /**
   * The dimensions of the mesh in each possible direction
   */
  std::vector<double> dimensions_ = {};

  /**
   * Vector of vertices which comprise the entire mesh
   */
  vertex_vec vertices_raw_ = {};

  /**
   * Unordered map linking each of the raw vertices with a unique global index
   */
  vertex_map vertices_ = {};

  /**
   * Total number of vertices in the mesh
   */
  long number_vertices_ = 0;

  /**
   * Total number of elements in the mesh
   */
  long number_elements_ = 0;

  /**
   * Checks whether the indexes being passed to the index accessor [] can be found
   * @param idxs Index array
   * @param maxs The index limits
   * @throws std::runtime_error if any idxs[k] > maxs[k]
   */
  static void check_index_in_bounds(const std::vector<long>& idxs, const std::vector<long>& maxs);

  /**
   * Checks whether the indexes being passed to the index accessor are correct in dimension
   * @param idx The indexes to be looked for in the vertices_map
   * @param dimension The spatial dimension of the mesh
   * @throws std::runtime_error if idx.size() != dimension
   */
  static void check_correct_size_indexes(const std::vector<long>& idx, const int& dimension);

  /**
   * Checks whether the indexes passed to the index operator `[]' are negative.
   * In case they are, they are translated to positive. This makes it possible inverted
   * index access, e.g. mesh[-1] (a la python).
   * @param idxs Index array
   * @param maxs The index limits
   */
  static void correct_negative_indexes(std::vector<long>& idxs, std::vector<long> maxs);

  /**
   * Transforms an array of cartesian indexes into a single global index to look for the
   * corresponding vertex in the vertex_map attribute
   * @param idxs Index array
   * @param shape
   * @return
   */
  virtual long cartesian_to_global_index(const std::vector<long>& idxs, const std::vector<long>& shape) const = 0;

  /**
   * Transforms a single global index into an array of cartesian indexes associated with the
   * corresponding vertex in the vertex_map attribute
   * @param idx The global index
   * @param shape
   * @return
   */
  virtual std::vector<long> global_index_to_cartesian(long idx, const std::vector<long>& shape) const = 0;

 public:
  // region Cttors:
  /**
   * Default constructor
   */
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
  virtual void plot() const = 0;

  // endregion

  // region Overloads:

  virtual const Vertex& operator[](const std::vector<long>& idx) const = 0;
  friend std::ostream& operator<<(std::ostream& os, const Mesh& mesh);

  // endregion
};

// endregion

// region Abstract Lattice:

class SUQMesh : public Mesh {
 protected:
  std::vector<long> idx_max_ = {};
  std::vector<double> origin_ = {};

  void initialise_dimensions(double dx, std::vector<double>& dimensions);
  long cartesian_to_global_index(const std::vector<long>& idxs, const std::vector<long>& shape) const final;
  std::vector<long> global_index_to_cartesian(long idx, const std::vector<long>& shape) const final;

 public:
  // region Cttors:
  explicit SUQMesh(double dx, std::vector<double>& dimensions, std::vector<double>& origin);
  // endregion

  // region Methods:
  double volume() const final;
  virtual double element_volume() const = 0;
  // endregion

  // region Operators:
  const Vertex& operator[](const std::vector<long>& idx) const final;
  // endregion
};

// endregion

}}
#endif
