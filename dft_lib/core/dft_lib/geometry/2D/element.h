#ifndef CLASSICALDFT_2D_BOX_H
#define CLASSICALDFT_2D_BOX_H

#include "dft_lib/geometry/base/element.h"

namespace dft_core {
namespace geometry {
namespace two_dimensional {

const double DEFAULT_BOX_LENGTH = 1;
const std::vector<double> DEFAULT_2D_BOX_ORIGIN = std::vector<double>{0, 0};

class SquareBox : public Element {
 private:
  double length_ = DEFAULT_BOX_LENGTH;
  void _initialise(double dx, const std::vector<double>& origin);

 public:
  SquareBox();
  explicit SquareBox(double length, const std::vector<double>& origin);
  explicit SquareBox(std::vector<Vertex>&& vertices);

  double volume() const override;
};

}}}

#endif