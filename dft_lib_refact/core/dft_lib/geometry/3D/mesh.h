#ifndef CLASSICALDFT_3D_MESH_H
#define CLASSICALDFT_3D_MESH_H

#include "dft_lib/geometry/base/mesh.h"
#include "dft_lib/geometry/3D/element.h"

namespace dft_core {
namespace geometry {
namespace three_dimensional {

typedef std::reference_wrapper<SquareBox> sqbox_refwrap;
typedef std::vector<SquareBox> sqbox_vec;
typedef std::unordered_map<int, sqbox_refwrap> sqbox_map;

class SUQMesh : public dft_core::geometry::SUQMesh
{
 private:
  // region Attributes

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

  // region Methods:

  void plot() const override;

  const std::vector<SquareBox>& elements() const;

  double element_volume() const final;

  // endregion
};

typedef SUQMesh Lattice;


}}}
#endif