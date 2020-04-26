//
// Created by migduroli on 19/04/2020.
//

#ifndef CLASSICALDFT_SALUTE_H
#define CLASSICALDFT_SALUTE_H

#include <string>

namespace dft_core {
  int MultiplyByTwo(int x);

  template <class T>
  static T& MultiplyByTwo(T& x)
  {
    return static_cast<T>(2.0*x);
  }

  std::string Salute(std::string& name);

  template <class T>
  static std::string GetType(T arg) {
    return typeid(arg).name();
  };
}

#endif
