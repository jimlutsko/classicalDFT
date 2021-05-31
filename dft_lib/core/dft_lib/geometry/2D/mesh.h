#ifndef CLASSICALDFT_2D_MESH_H
#define CLASSICALDFT_2D_MESH_H

#include "dft_lib/geometry/base/mesh.h"
#include "dft_lib/geometry/2D/element.h"

namespace dft_core {
namespace geometry {
namespace two_dimensional {

typedef std::reference_wrapper<SquareBox> sqbox_refwrap;
typedef std::vector<SquareBox> sqbox_vec;
typedef std::unordered_map<int, sqbox_refwrap> sqbox_map;

 class SUQMesh : public dft_core::geometry::SUQMesh
{
 private:
  // region Attributes:

  /**
   * Vector of elements which constitutes the mesh
   */
  sqbox_vec elements_raw_ = {};

  /**
   * Dictionary which maps a global index with every mesh element
   */
  sqbox_map elements_ = {};

  // endregion

  void initialise(double dx);

 public:
  // region Cttors:

  explicit SUQMesh(double dx, std::vector<double>& dimensions, std::vector<double>& origin);

  // endregion

  // region Overloads:

  /**
   * Gets the graphical representation of the SUQMesh
   */
  void plot() const override ;

  /**
   * Gets a reference to the vector which contains the mesh elements
   * @return
   */
  const std::vector<SquareBox>& elements() const;

  /**
   * Gets the volume of the fundamental element of the mesh
   * @return
   */
  double element_volume() const final;

  // endregion
};

typedef SUQMesh Lattice;

}}}

#endif
