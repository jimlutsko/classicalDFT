#ifndef CLASSICALDFT_2D_MESH_H
#define CLASSICALDFT_2D_MESH_H

#include "dft_lib/geometry/base/mesh.h"
#include "dft_lib/geometry/2D/element.h"

namespace dft_core {
namespace geometry {
namespace two_dimensional {

const int DEFAULT_NUMBER_VERTICES = 1000;

typedef std::reference_wrapper<SquareBox> element_refwrap;
typedef std::vector<SquareBox> element_vec;
typedef std::unordered_map<int, element_refwrap> element_map;

class SUQMesh : public Mesh
{
 private:

  element_vec elements_raw_ = {};
  element_map elements_ = {};
  std::vector<long> idx_max_ = {0, 0};
  std::vector<double> origin_ = {0, 0};

 public:
  // region Cttors:
  explicit SUQMesh(double dx, double length, std::vector<double>& origin);
  // endregion

  // region Overwrites:
  double volume() const override ;
  const std::vector<SquareBox>& elements() const {return elements_raw_;}
  const Vertex& operator[](const std::vector<long>& idx) const override;
  // endregion
};

typedef SUQMesh Lattice;

}}}

#endif
