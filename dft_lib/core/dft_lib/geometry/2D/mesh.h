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
  sqbox_vec elements_raw_ = {};
  sqbox_map elements_ = {};

  void initialise(double dx);

 public:
  // region Cttors:
  explicit SUQMesh(double dx, std::vector<double>& dimensions, std::vector<double>& origin);
  // endregion

  // region Overwrites:
  void plot() const override ;

  const std::vector<SquareBox>& elements() const;
  double element_volume() const final;
  // endregion
};

typedef SUQMesh Lattice;

}}}

#endif
