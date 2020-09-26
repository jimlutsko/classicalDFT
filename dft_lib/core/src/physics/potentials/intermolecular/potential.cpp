#include "dft_lib/physics/potentials/intermolecular/potential.h"
#include "dft_lib/utils/functions.h"

#include <cmath>

#include "dft_lib/numerics/integration.h"
using namespace dft_core::numerics::integration;
using namespace dft_core::utils;

namespace dft_core {
namespace physics {
namespace potentials {
namespace intermolecular {

// region Potential (base class):

Potential::Potential(double sigma, double epsilon, double r_cutoff)
    : sigma_(sigma),
      epsilon_(epsilon),
      r_cutoff_(r_cutoff),
      bh_perturbation_(false),
      epsilon_shift_(DEFAULT_ZERO),
      r_min_(DEFAULT_ZERO),
      v_min_(DEFAULT_ZERO),
      r_attractive_min_(DEFAULT_ZERO),
      r_zero_(DEFAULT_LENGTH_SCALE),
      kT_(DEFAULT_ENERGY_SCALE) {}

double Potential::sigma() const { return sigma_; }
double Potential::epsilon() const { return epsilon_; }
double Potential::epsilon_shift() const { return epsilon_shift_; }
double Potential::r_cutoff() const { return r_cutoff_; }
double Potential::r_min() const { return r_min_; }
double Potential::v_min() const { return v_min_; }
double Potential::r_attractive_min() const { return r_attractive_min_; }
double Potential::r_zero() const { return r_zero_; }
bool Potential::bh_perturbation() const { return bh_perturbation_; }
double Potential::kT() const { return kT_; }
const PotentialName& Potential::id() const { return potential_id_; }

std::string Potential::identifier() const
{
  std::string name;
  switch (this->id()) {
    case PotentialName::LennardJones:
      name = "LennardJones";
      break;
    case PotentialName::tenWoldeFrenkel:
      name = "tenWoldeFrenkel";
      break;
    case PotentialName::WangRamirezDobnikarFrenkel:
      name = "WangRamirezDobnikarFrenkel";
      break;
    default:
      name = "Unknown";
  }

  return name + "_"
         + std::to_string(this->sigma()) + "_"
         + std::to_string(this->epsilon()) + "_"
         + std::to_string(this->r_cutoff()) + "_"
         + std::to_string(this->r_attractive_min()) + "_"
         + std::to_string(this->bh_perturbation());
}

double Potential::_v_potential(double r) const { return vr_(r) - epsilon_shift(); }

double Potential::v_potential(double r) const { return _v_potential(r); }

std::vector<double> Potential::v_potential(const std::vector<double>& r) const {
  return functions::apply_vector_wise<Potential, double>(*this, &Potential::_v_potential, r);
}

arma::vec Potential::v_potential(const arma::vec& r) const
{
  return arma::vec(this->v_potential(arma::conv_to<std::vector<double>>::from(r)));
}

double Potential::_v_potential_r2(double r_squared) const { return vr2_(r_squared) - epsilon_shift(); }

double Potential::v_potential_r2(double r_squared) const { return _v_potential_r2(r_squared); }

std::vector<double> Potential::v_potential_r2(const std::vector<double>& r_squared) const
{
  return functions::apply_vector_wise<Potential, double>(*this, &Potential::_v_potential_r2, r_squared);
}

arma::vec Potential::v_potential_r2(const arma::vec& r_squared) const
{
  return arma::vec(this->v_potential_r2(arma::conv_to<std::vector<double>>::from(r_squared)));
}

double Potential::_w_repulsive(double r) const
{
  if (this->bh_perturbation_) {
    return (r < r_zero() ? v_potential(r) : 0.0);
  }

  return (r < r_min_ ? v_potential(r) - v_min() : 0.0);
}

double Potential::w_repulsive(double r) const { return _w_repulsive(r); }

std::vector<double> Potential::w_repulsive(const std::vector<double>& r) const
{
  return functions::apply_vector_wise<Potential, double>(*this, &Potential::_w_repulsive, r);
}

arma::vec Potential::w_repulsive(const arma::vec& r) const
{
  return arma::vec(this->w_repulsive(arma::conv_to<std::vector<double>>::from(r)));
}

void Potential::SetWCALimit(double r) { r_attractive_min_ = r; }

void Potential::SetBHPerturbation()
{
  bh_perturbation_ = true;
  r_attractive_min_ = r_zero();
}

double Potential::bh_diameter_kernel(double r) const
{
  return (1.0 - std::exp(-w_repulsive(r) / kT_));
}

double Potential::_w_attractive(double r) const { return this->w_attractive_r2(r * r); }

double Potential::w_attractive(double r) const { return _w_attractive(r); }

std::vector<double> Potential::w_attractive(const std::vector<double>& r) const
{
  return functions::apply_vector_wise<Potential, double>(*this, &Potential::_w_attractive, r);
}

arma::vec Potential::w_attractive(const arma::vec& r) const
{
  return arma::vec(this->w_attractive(arma::conv_to<std::vector<double>>::from(r)));
}

double Potential::w_attractive_r2(double r_squared) const
{
  double ret = 0.0;

  // zero outside of cutoff
  if (r_squared < (r_cutoff() * r_cutoff())) {
    if (bh_perturbation_) {
      return (r_squared < r_zero() * r_zero() ? 0.0 : v_potential_r2(r_squared));
    } else if (r_squared < (r_attractive_min() * r_attractive_min())) {
      // Our "generalized" WCA
      return 0.0;
    } else if (r_squared < (r_min() * r_min())) {
      // WCA continuation inside r_min
      return v_min();
    } else {
      // Just the potential outsize of r_min
      return v_potential_r2(r_squared);
    }
  }
  return ret;
}

double Potential::vdw_kernel(double r) const { return r*r*w_attractive(r); }

double Potential::FindHardSphereDiameter(double kT)
{
  kT_ = kT;
  auto dHardCore = FindHardCoreDiameter();
  auto rLimit = bh_perturbation_ ? r_zero() : r_min();
  auto integrator = Integrator<Potential>(*this, &Potential::bh_diameter_kernel, 1e-4, 1e-6);
  auto dHardSphere = dHardCore + integrator.DefiniteIntegral(dHardCore, rLimit);
  return dHardSphere;
}

double Potential::ComputeVanDerWaalsIntegral(double kT)
{
  kT_ = kT;
  auto prefactor = (2*M_PI/kT);
  auto integrator = Integrator<Potential>(*this, &Potential::vdw_kernel, 1e-6, 1e-8);

  auto limit_superior = std::max(r_cutoff(), 0.0);
  auto integral = bh_perturbation() ?
                  integrator.DefiniteIntegral(r_zero(), limit_superior)
                                    : (integrator.DefiniteIntegral(r_attractive_min(), limit_superior)
                                       + integrator.DefiniteIntegral(r_min(), limit_superior));

  return prefactor * integral;
}

double Potential::operator()(double r) const { return v_potential(r); }

std::vector<double> Potential::operator()(const std::vector<double>& r) const { return v_potential(r); }

arma::vec Potential::operator()(const arma::vec& r) const { return v_potential(r); }

//endregion

//region LennardJones (12-6):

LennardJones::LennardJones(): Potential()
{
  potential_id_ = PotentialName::LennardJones;
  epsilon_shift_ = (r_cutoff_ < 0 ? 0.0 : this->vr_(r_cutoff_));
  r_min_ = this->FindRMin();
  v_min_ = this->v_potential(r_min_);
  r_zero_ = std::pow(0.5 * std::sqrt(1 + epsilon_shift_) + 0.5, -1.0/6.0);
}

LennardJones::LennardJones(double sigma, double epsilon, double r_cutoff)
    : Potential(sigma, epsilon, r_cutoff)
{
  potential_id_ = PotentialName::LennardJones;
  epsilon_shift_ = (r_cutoff_ < 0 ? 0.0 : this->vr_(r_cutoff_));
  r_min_ = this->FindRMin();
  v_min_ = this->v_potential(r_min_);
  r_zero_ = std::pow(0.5 * std::sqrt(1 + epsilon_shift_) + 0.5, -1.0/6.0);
}

double LennardJones::vr_(double r) const
{
  auto y = sigma_ / r;
  auto y6 = y * y * y * y * y * y;
  return 4 * epsilon_ * (y6 * y6 - y6);
}

double LennardJones::vr2_(double r2) const
{
  auto y2 = sigma_ * sigma_ / r2;
  auto y6 = y2 * y2 * y2;
  return 4 * epsilon_ * (y6 * y6 - y6);
}

double LennardJones::FindRMin() const { return std::pow(2.0, 1.0/6.0) * sigma_; }

double LennardJones::FindHardCoreDiameter() const { return 0; }

//endregion

//region tenWoldeFrenkel:

tenWoldeFrenkel::tenWoldeFrenkel(): Potential()
{
  potential_id_ = PotentialName::tenWoldeFrenkel;
  alpha_ = DEFAULT_ALPHA_PARAMETER;
  r_min_ = this->FindRMin();
  v_min_  = this->v_potential(r_min_);
  r_zero_ = std::sqrt(1 + std::pow(25 * sqrt(1+epsilon_shift_) + 25, -1.0/3.0));
  epsilon_shift_ = (r_cutoff_ <= 0.0 ? 0.0 : this->vr_(r_cutoff_));
}

tenWoldeFrenkel::tenWoldeFrenkel(double sigma, double epsilon, double r_cutoff, double alpha)
    : Potential(sigma, epsilon, r_cutoff), alpha_(alpha)
{
  potential_id_ = PotentialName::tenWoldeFrenkel;
  epsilon_shift_ = (r_cutoff_ <= 0.0 ? 0.0 : this->vr_(r_cutoff_));
  r_min_ = this->FindRMin();
  v_min_  = this->v_potential(r_min_);
  r_zero_ = std::sqrt(1 + std::pow(25 * sqrt(1+epsilon_shift_) + 25, -1.0/3.0));
}
double tenWoldeFrenkel::alpha() const { return alpha_; }

double tenWoldeFrenkel::vr_(double r) const
{
  if (r < sigma_) { return MAX_POTENTIAL_VALUE; }

  double s = r / sigma_;
  double y = 1.0/(s * s - 1);
  double y3 = y * y * y;

  return (4 * epsilon_ / (alpha_ * alpha_)) * (y3 * y3 - alpha_ * y3);
}

double tenWoldeFrenkel::vr2_(double r2) const
{
  if(r2 < sigma_*sigma_) return MAX_POTENTIAL_VALUE;

  double s2 = r2 / (sigma_ * sigma_);
  double y = 1.0 / (s2 - 1);
  double y3 = y * y * y;

  return (4 * epsilon_ / (alpha_ * alpha_)) * (y3 * y3 - alpha_ * y3);
}

double tenWoldeFrenkel::FindHardCoreDiameter() const { return sigma_; }

double tenWoldeFrenkel::FindRMin() const
{
  return sigma_ * std::sqrt(1 + std::pow(2.0 /alpha_, 1.0/3.0));
}

//endregion

//region WangRamirezDobnikarFrenkel:

double WangRamirezDobnikarFrenkel::vr_(double r) const
{
  double y = sigma_ / r;
  double z = r_cutoff_ / r;
  return (r < r_cutoff_) ? epsilon_ * (y * y - 1) * (z * z - 1) * (z * z - 1) : 0.0;
}
double WangRamirezDobnikarFrenkel::vr2_(double r2) const
{
  double y = sigma_ * sigma_ / r2;
  double z = r_cutoff_ * r_cutoff_ / r2;
  return (r2 < r_cutoff_ * r_cutoff_) ? epsilon_ * (y - 1) * (z - 1) * (z - 1) : 0.0;
}

WangRamirezDobnikarFrenkel::WangRamirezDobnikarFrenkel(): Potential()
{
  potential_id_ = PotentialName::WangRamirezDobnikarFrenkel;

  // The cutoff is built in:
  r_cutoff_ = 3;
  epsilon_ *= 2 * std::pow(r_cutoff_ / sigma_, 2)
              * std::pow(2 * ((r_cutoff_ / sigma_) * (r_cutoff_ / sigma_)-1)/3.0, -3.0);
  epsilon_shift_ = 0.0;

  r_min_ = this->FindRMin();
  v_min_ = this->v_potential(r_min_);
  r_zero_ = 1.0;
}

WangRamirezDobnikarFrenkel::WangRamirezDobnikarFrenkel(double sigma, double epsilon, double r_cutoff)
    : Potential(sigma, epsilon, r_cutoff)
{
  potential_id_ = PotentialName::WangRamirezDobnikarFrenkel;

  // The cutoff is built in:
  epsilon_ *= 2 * std::pow(r_cutoff / sigma, 2)
                * std::pow(2 * ((r_cutoff / sigma) * (r_cutoff / sigma)-1)/3.0, -3.0);
  epsilon_shift_ = 0.0;

  r_min_  = this->FindRMin();
  v_min_ = this->v_potential(r_min_);
  r_zero_ = 1.0;
}

double WangRamirezDobnikarFrenkel::FindHardCoreDiameter() const { return 0; }

double WangRamirezDobnikarFrenkel::FindRMin() const
{
  return r_cutoff_ * std::pow((1.0 + 2.0*(r_cutoff_/sigma_) * (r_cutoff_/sigma_))/3.0, -0.5);
}

//endregion
}}}}