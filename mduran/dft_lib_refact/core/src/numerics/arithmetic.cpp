#include "dft_lib/numerics/arithmetic.h"

#include <vector>
#include <cmath>
#include <numeric>
#include <iostream>

namespace dft_core
{
namespace numerics
{
namespace arithmetic
{
namespace summation {

return_type KahanBabuskaAdd(x_type x, x_type sum, error_type error)
{
  auto c_error = error.front();
  auto x_corrected = x - c_error;
  auto t = sum + x_corrected;
  c_error = (t - sum) - x_corrected;
  sum = t;
  return return_type{ sum, error_type{c_error} };
}

return_type KahanBabuskaSum(const std::vector<x_type>& x_series, const x_type& sum_ini, const error_type & error_ini)
{
  // The method is given in here:
  // wikipedia: https://en.wikipedia.org/wiki/Kahan_summation_algorithm#The_algorithm

  x_type sum = sum_ini;
  error_type error = error_ini;

  for (const auto& x : x_series)
  {
    std::tie(sum, error) = KahanBabuskaAdd(x, sum, error);
  }

  // When catching the result we need to use std::tie(x, y) [https://en.cppreference.com/w/cpp/utility/tuple]
  return return_type{ sum, error };
}

return_type KahanBabuskaNeumaierAdd(x_type x, x_type sum, error_type error)
{
  auto t = sum + x;
  auto c_error = error.front();
  if (std::abs(sum) >= std::abs(x))
  {
    c_error += (sum - t) + x;
  }
  else
  {
    c_error += (x - t) + sum;
  }
  sum = t;
  return return_type{ sum, error_type{c_error} };
}

return_type KahanBabuskaNeumaierSum(const std::vector<x_type>& x_series, const x_type& sum_ini, const error_type & error_ini)
{
  // The method is given in here:
  // wikipedia: https://en.wikipedia.org/wiki/Kahan_summation_algorithm#Further_enhancements

  auto sum = sum_ini;
  auto error = error_ini;

  for (const auto& x : x_series)
  {
    std::tie(sum, error) = KahanBabuskaNeumaierAdd(x, sum, error);
  }

  sum += std::accumulate(error.begin(), error.end(), 0.0);
  return return_type{ sum, error };
}

return_type KahanBabuskaKleinAdd(x_type x, x_type sum, error_type error)
{
  auto t = sum + x;
  auto err = error;
  x_type c;

  if (std::abs(sum) >= std::abs(x)) {
    c = (sum - t) + x;
  }
  else
  {
    c = (x - t) + sum;
  }

  double cc;
  sum = t;
  t = err[0] + c;

  if (std::abs(err[0]) >= std::abs(c))
  {
    cc = (err[0] - t) + c;
  }
  else
  {
    cc = (c - t) + err[0];
  }

  err[0] = t;
  err[1] += cc;

  return return_type{sum, err};
}

return_type KahanBabuskaKleinSum(const std::vector<x_type>& x_series, const x_type& sum_ini, const error_type & error_ini)
{
  // The method is given in here:
  // wikipedia: https://en.wikipedia.org/wiki/Kahan_summation_algorithm#Further_enhancements

  auto sum = sum_ini;
  auto error = error_ini;

  for (const auto& x : x_series)
  {
    std::tie(sum, error) = KahanBabuskaKleinAdd(x, sum, error);
  }

  sum += std::accumulate(error.begin(), error.end(), 0.0);

  return return_type{ sum, error };
}

const error_type& CompensatedSum::error() const { return this->error_; }
const summation::Type& CompensatedSum::type() const { return this->type_; }

CompensatedSum::CompensatedSum(summation::Type type)
: sum_(0.0), type_(type)
{
  if (type == summation::Type::KahanBabuskaKlein) {
    error_ = std::vector<double>{0.0, 0.0};
    method_ = &KahanBabuskaKleinAdd;
  }
  else if (type == summation::Type::KahanBabuska){
    error_ = std::vector<double>{0.0};
    method_ = &KahanBabuskaAdd;
  }
  else if (type == summation::Type::KahanBabuskaNeumaier){
    error_ = std::vector<double>{0.0};
    method_ = &KahanBabuskaNeumaierAdd;
  }
  else {
    console::Warning("Method not implemented yer, defaulting to Kahan-Babuska-Neumaier");
    error_ = std::vector<double>{0.0};
    method_ = &KahanBabuskaNeumaierAdd;
  }
}

void CompensatedSum::add(x_type value)
{
  std::tie(sum_, error_) = this->method_(value, sum_, error_);
}

x_type CompensatedSum::Sum() const
{
  if (this->type_ == summation::Type::KahanBabuska) { return this->sum_; }
  return sum_ + std::accumulate(error_.begin(),error_.end(), 0.0);
}

CompensatedSum::operator double() const
{
  return this->Sum();
}

CompensatedSum& CompensatedSum::operator+=(const x_type x)
{
  this->add(x);
  return *this;
}

CompensatedSum& CompensatedSum::operator-=(const x_type x)
{
  this->add(-x);
  return *this;
}

CompensatedSum& CompensatedSum::operator+=(const CompensatedSum& other)
{
  sum_ += other.sum_;
  for (auto k = 0; k < error_.size(); k++)
  {
    this->error_[k] += other.error_[k];
  }
  return *this;
}

CompensatedSum& CompensatedSum::operator-=(const CompensatedSum& other)
{
  sum_ -= other.sum_;
  for (auto k = 0; k < error_.size(); k++)
  {
    this->error_[k] -= other.error_[k];
  }
  return *this;
}

CompensatedSum& CompensatedSum::operator+(x_type x)
{
  this->add(x);
  return *this;
}

CompensatedSum& CompensatedSum::operator-(x_type x)
{
  this->add(-x);
  return *this;
}

CompensatedSum& CompensatedSum::operator+=(const std::vector<x_type>& x)
{
  for (const auto& e : x)
  {
    this->add(e);
  }
  return *this;
}

CompensatedSum& CompensatedSum::operator-=(const std::vector<x_type>& x)
{
  for (const auto& e : x)
  {
    this->add(-e);
  }
  return *this;
}

CompensatedSum& CompensatedSum::operator+(const CompensatedSum& other) {
  sum_ += other.sum_;
  for (auto k = 0; k < error_.size(); k++)
  {
    this->error_[k] += other.error_[k];
  }
  return *this;
}

CompensatedSum& CompensatedSum::operator-(const CompensatedSum& other) {
  sum_ -= other.sum_;
  for (auto k = 0; k < error_.size(); k++)
  {
    this->error_[k] -= other.error_[k];
  }
  return *this;
}

CompensatedSum& CompensatedSum::operator=(x_type x)
{
  sum_ = x;
  error_ = (type_ == summation::Type::KahanBabuskaKlein)
      ? error_type{0.0, 0.0}
      : error_type{0.0};
  return *this;
}

std::ostream& operator<<(std::ostream& output, const CompensatedSum& obj)
{
  output << obj.Sum();
  return output;
}

}}}}