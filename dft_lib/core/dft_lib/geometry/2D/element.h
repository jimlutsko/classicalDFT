#ifndef CLASSICALDFT_2DBOX_H
#define CLASSICALDFT_2DBOX_H

#include "dft_lib/geometry/base/element.h"

namespace dft_core {
namespace geometry {
namespace two_dimensional {

const int DEFAULT_DIMENSION = 2;
const double DEFAULT_LENGTH = 1;
const std::vector<double> DEFAULT_ORIGIN = std::vector<double>{0, 0};

class SquareBox : public Element {
 private:
  double length_ = DEFAULT_LENGTH;

 public:
  SquareBox();
  explicit SquareBox(double length, const std::vector<double>& origin);
};

}}}

#endif