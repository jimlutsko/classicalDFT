#ifndef CLASSICALDFT_VECTORIZE_H
#define CLASSICALDFT_VECTORIZE_H

#include <iostream>
#include <functional>
#include <memory>

namespace dft_core {
namespace utils {
namespace functions {

template<class T, typename return_type, typename input_type = return_type>
using class_method = std::function<return_type(const T&, input_type)>;

template<typename return_type, typename input_type = return_type>
using general_method = std::function<return_type(input_type)>;

template<typename x_type = double,
    typename = typename std::enable_if<std::is_arithmetic<x_type>::value, x_type>::type
>
std::vector<x_type> apply_vector_wise(
    general_method<x_type, x_type> func,
    const std::vector<x_type>& x
)
{
  auto y = std::vector<x_type>();
  for (const auto& k : x) { y.push_back(func(k)); }
  return y;
}

template<class T,
    typename x_type = double,
    typename = typename std::enable_if<std::is_arithmetic<x_type>::value, x_type>::type
>
std::vector<x_type> apply_vector_wise(
    const T& obj,
    class_method<T, x_type, x_type> method,
    const std::vector<x_type>& x
)
{
  auto y = std::vector<x_type>();
  for (const auto& k : x) { y.push_back(method(obj, k)); }
  return y;
}

}}}
#endif  // CLASSICALDFT_VECTORIZE_H
