#include <dft_lib/geometry/base/mesh.h>

#include <boost/range/combine.hpp>
#include <numeric>
#include <cmath> // std::abs

namespace dft_core {
namespace geometry {

// region Abstract Mesh:

const std::vector<long>& dft_core::geometry::Mesh::shape() const { return shape_; }
const std::vector<double>& dft_core::geometry::Mesh::dimensions() const { return dimensions_; }
long Mesh::number_vertices() const { return number_vertices_; }
const vertex_vec& Mesh::vertices() const { return vertices_raw_; }

std::ostream& operator<<(std::ostream& os, const Mesh& mesh)
{
  os << "Mesh object [" << std::addressof(mesh) << "]" << std::endl;
  os << "\u2b91 Volume: " << mesh.volume() << std::endl;
  os << "\u2b91 Number of vertices: " << mesh.number_vertices() << std::endl;

  os << "\u2b91 Shape: [" ;
  for (auto k = 0; k < mesh.shape().size(); k++)
  {
    os << mesh.shape()[k];
    if (k != (mesh.shape().size()-1)) { os << ", "; }
  }
  os << "]" << std::endl;

  os << "\u2b91 Dimensions: [";
  for (auto k = 0; k < mesh.dimensions().size(); k++)
  {
    os << mesh.dimensions()[k];
    if (k != (mesh.shape().size()-1)) { os << ", "; }
  }
  os << "]" << std::endl;

  return os;
}


// region Indexing functionality:

void Mesh::check_index_in_bounds(const std::vector<long>& idxs, const std::vector<long>& maxs)
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

void Mesh::check_correct_size_indexes(const std::vector<long>& idx, const int& dimension)
{
  if (idx.size() != dimension) {
    throw std::runtime_error("[!] The index array cannot have more than 2 components");
  }
}

void Mesh::correct_negative_indexes(std::vector<long>& idxs, std::vector<long> maxs)
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

// endregion

// region SUQ Mesh
const static double _scaling_dx = 1E-8;

void SUQMesh::initialise_dimensions(double dx, std::vector<double>& dimensions)
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
}

SUQMesh::SUQMesh(double dx, std::vector<double>& dimensions, std::vector<double>& origin)
    : origin_(std::move(origin))
{
  initialise_dimensions(dx, dimensions);
}

long SUQMesh::cartesian_to_global_index(const std::vector<long>& idxs, const std::vector<long>& shape) const
{
  long v = idxs[static_cast<unsigned long>(Direction::X)];

  for (auto k = 0; k < idxs.size()-1; k++)
  {
    v = (idxs[k+1] + shape[k+1] * v);
  }
  return v;
}

std::vector<long> SUQMesh::global_index_to_cartesian(long pos, const std::vector<long>& shape) const
{
  auto v = std::vector<long>(shape_.size());
  auto k_max = v.size()-1;
  for (auto k = k_max; k > 0; k--)
  {
    v[k] = pos % (shape_[k]);
    if (1 < k) { pos = (pos - v[k]) / shape_[k]; }
    else { v[k-1] = pos / shape_[k]; }
    //  k = pos%(Nz_); pos = (pos-k)/Nz_;
    //  j = pos%Ny_; i = pos/Ny_;}
  }
  return v;
}

double SUQMesh::volume() const
{
  return std::accumulate(begin(dimensions_), end(dimensions_), 1, std::multiplies<>());
}

const Vertex& SUQMesh::operator[](const std::vector<long>& idx) const
{
  auto idxs = idx;
  dft_core::geometry::Mesh::correct_negative_indexes(idxs, idx_max_);
  dft_core::geometry::Mesh::check_correct_size_indexes(idxs, shape_.size());
  dft_core::geometry::Mesh::check_index_in_bounds(idxs, idx_max_);
  return vertices_.at(cartesian_to_global_index(idxs, shape_));
}

// endregion
}}