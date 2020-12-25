#include "dft_lib/geometry/2D/mesh.h"

#include <boost/range/combine.hpp>
#include <utility>

#include "dft_lib/geometry/2D/element.h"

namespace dft_core {
namespace geometry {
namespace two_dimensional {

const static double _scaling_dx = 1E-8;

SUQMesh::SUQMesh(double dx, double length, std::vector<double>& origin)
    : origin_(std::move(origin))
{
  dimensions_ = {length, length};
  shape_ = {
      static_cast<long>((dimensions_[0] + _scaling_dx * dx) / dx) + 1,
      static_cast<long>((dimensions_[1] + _scaling_dx * dx) / dx) + 1,
  };
  number_vertices_ = shape_[0] * shape_[1];
  idx_max_ = {shape_[0] - 1, shape_[1] - 1};

  //region Vertices init:
  auto vertex_index = 0;
  auto element_index = 0;

  auto x = origin_[0]; auto y = origin_[1];
  vertices_raw_ = vertex_vec(number_vertices_);
  elements_raw_ = element_vec((shape_[0]-1)*(shape_[1]-1));

  for (auto y_idx = 0; y_idx < shape_[1]; y_idx++)
  {
    for (auto x_idx = 0; x_idx < shape_[0]; x_idx++)
    {
      vertices_raw_[vertex_index] = Vertex({x, y});
      vertices_.insert({vertex_index, std::ref(vertices_raw_[vertex_index])});

      if ((x_idx < idx_max_[0]) && (y_idx < idx_max_[1]))
      {
        elements_raw_[element_index] = SquareBox(dx, {x, y});
        element_index += 1;
      }

      vertex_index += 1;
      x += dx;
    }
    x = 0; y += dx;
  }
  //endregion
}

static long _transform_coordinates_to_index(const std::vector<long>& idxs, const std::vector<long>& shape)
{
  return (idxs[0] * shape[0] + idxs[1]);
}

// region Indexing functionality:
static void _check_index_in_bounds(const std::vector<long>& idxs, const std::vector<long>& maxs)
{
  auto zip_list = boost::combine(idxs, maxs);
  auto test = false;
  for (auto tup : zip_list)
  {
    double k, k_max; boost::tie(k,k_max) = tup;
    if (std::abs(k) > k_max + 1) { test = true; break;}
  }

  if (test) {throw std::runtime_error("[!] Indexes are out of bound in mesh-indexer");}
}

static void _check_correct_size_indexes(const std::vector<long>& idx, const int& dimension)
{
  if (idx.size() != dimension) {
    throw std::runtime_error("[!] The index array cannot have more than 2 components");
  }
}

static void _correct_negative_indexes(std::vector<long>& idxs, std::vector<long> maxs)
{
  auto zip_list = boost::combine(idxs, maxs);
  for (auto tup : zip_list)
  {
    double k, k_max; boost::tie(k,k_max) = tup;
    if ((k < 0) && (std::abs(k) <= k_max + 1)) {
      tup.get<0>() += tup.get<1>() + 1;
    }
  }
}

// endregion

const Vertex& SUQMesh::operator[](const std::vector<long>& idx) const
{
  auto idxs = idx;
  _correct_negative_indexes(idxs, idx_max_);
  _check_correct_size_indexes(idxs, shape_.size());
  _check_index_in_bounds(idxs, idx_max_);
  return vertices_.at(_transform_coordinates_to_index(idxs, shape_));
}

double SUQMesh::volume() const { return dimensions_[0] * dimensions_[1]; }

}}}