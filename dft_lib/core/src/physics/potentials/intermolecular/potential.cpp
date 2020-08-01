#include "dft_lib/physics/potentials/intermolecular/potential.h"


namespace dft_core
{
namespace physics
{
namespace potentials
{
namespace intermolecular {

Potential::Potential(double sigma, double epsilon, double r_cutoff)
    : sigma_(sigma),
      epsilon_(epsilon),
      r_cutoff_(r_cutoff),
      bh_perturbation_(false),
      eps_shift_(DEFAULT_ZERO),
      r_min_(DEFAULT_ZERO),
      v_min_(DEFAULT_ZERO),
      r_attractive_min_(DEFAULT_ZERO),
      r_zero_(DEFAULT_LENGTH_SCALE),
      kT_(DEFAULT_ENERGY_SCALE) {}

double Potential::sigma() const { return sigma_; }
double Potential::epsilon() const { return epsilon_; }
double Potential::eps_shift() const { return eps_shift_; }
double Potential::r_cutoff() const { return r_cutoff_; }
double Potential::r_min() const { return r_min_; }
double Potential::v_min() const { return v_min_; }
double Potential::r_attractive_min() const { return r_attractive_min_; }
double Potential::r_zero() const { return r_zero_; }
bool Potential::bh_perturbation() const { return bh_perturbation_; }
double Potential::kT() const { return kT_; }
const PotentialName& Potential::id() const { return potential_id_; }

std::string Potential::identifier() const {
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
      name = "LennardJones";
  }

  return name + "_" + std::to_string(this->sigma()) + "_" + std::to_string(this->epsilon()) + "_" +
         std::to_string(this->r_cutoff()) + "_" + std::to_string(this->r_attractive_min()) + "_" +
         std::to_string(this->bh_perturbation());
}

}}}}