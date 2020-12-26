#include "dft_lib/geometry/2D/mesh.h"
#include "dft_lib/geometry/2D/element.h"
#include "dft_lib/graph/grace.h"

#include <boost/range/combine.hpp>
#include <numeric>
#include <utility>

namespace dft_core {
namespace geometry {
namespace two_dimensional {

const static double _scaling_dx = 1E-8;

SUQMesh::SUQMesh(double dx, std::vector<double>& dimensions, std::vector<double>& origin)
    : origin_(std::move(origin))
{
  dimensions_ = std::move(dimensions);

  for (auto l : dimensions_)
  {
    auto s = static_cast<long>((l + _scaling_dx * dx) / dx) + 1;
    shape_.push_back(s);
    idx_max_.push_back(s - 1);
  }

  number_vertices_ = std::accumulate(begin(shape_), end(shape_), 1, std::multiplies<>());
  number_elements_ = std::accumulate(begin(idx_max_), end(idx_max_), 1, std::multiplies<>());

  //region Vertices init:
  auto vertex_index = 0;
  auto element_index = 0;

  auto x = origin_[0]; auto y = origin_[1];
  vertices_raw_ = vertex_vec(number_vertices_);
  elements_raw_ = element_vec(number_elements_);

  for (auto j_idx = 0; j_idx < shape_[1]; j_idx++)
  {
    for (auto i_idx = 0; i_idx < shape_[0]; i_idx++)
    {
      vertices_raw_[vertex_index] = Vertex({x, y});
      vertices_.insert({vertex_index, std::ref(vertices_raw_[vertex_index])});

      if ((i_idx < idx_max_[0]) && (j_idx < idx_max_[1]))
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

void SUQMesh::plot() const
{
  auto g = dft_core::grace_plot::Grace();
  for (const auto& v : vertices_raw_)
  { g.AddPoint(v.coordinates()[0], v.coordinates()[1]); }

  auto dx = std::vector<double>{
    0.1 * dimensions_[0], 0.1 * dimensions_[1]
  };

  g.SetLimits(
      {origin_[0]-dx[0], (dimensions_[0]+origin_[0]) + dx[0]},
      {origin_[1]-dx[1], (dimensions_[1]+origin_[1]) + dx[1]}
      );

  g.SetLineType(dft_core::grace_plot::LineStyle::NO_LINE, 0);
  g.SetSymbol(dft_core::grace_plot::Symbol::SQUARE, 0);
  g.SetSymbolColor(dft_core::grace_plot::Color::BLUE, 0);
  g.SetSymbolFill(dft_core::grace_plot::Color::DARKGREEN, 0);

  g.RedrawAndWait();
}
}}}