#ifndef CLASSICALDFT_POTENTIAL_H
#define CLASSICALDFT_POTENTIAL_H

#include <iostream>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>

namespace dft_core
{
namespace physics
{
namespace potentials
{

const double DEFAULT_LENGTH_SCALE = 1.0;
const double DEFAULT_ENERGY_SCALE = 1.0;
const double DEFAULT_ZERO = 0.0;

namespace intermolecular
{

enum class PotentialName
{
  LennardJones = 0,
  tenWoldeFrenkel,
  WangRamirezDobnikarFrenkel
};

/**
 * @brief Potential-energy general model: Base for all the interaction potentials.
 * The potential class brings the ability to split the potential contribution into two parts:
 * - hard sphere and,
 * - attractive,
 * and, also comes with the functionality to compute the hard-sphere diameter.
 */
class Potential
{
 private:
  //region Attributes:

  /// The potential length scale
  double sigma_ = DEFAULT_ENERGY_SCALE;
  /// The potential energy scale
  double epsilon_ = DEFAULT_ENERGY_SCALE;
  /// The amount of energy the potential will be shifted
  double epsilon_shift_ = DEFAULT_ZERO;
  /// The cut-off distance, as of which the potential will be truncated
  double r_cutoff_ = -DEFAULT_LENGTH_SCALE;
  /// The position where the potential reaches its minimum
  double r_min_ = DEFAULT_LENGTH_SCALE;
  /// The minimum of the potential energy
  double v_min_ = -DEFAULT_ENERGY_SCALE;
  /// Weeks-Chandler-Andersen continuation limit
  double r_attractive_min_ = DEFAULT_ZERO;
  /// The position at which potential goes to zero
  double r_zero_ = DEFAULT_LENGTH_SCALE;
  /// Boolean variable to indicate whether Barker-Henderson model (repulsive-attractive split) is applied
  bool bh_perturbation_ = false;
  /// The non-dimensional energy (used to pass the temperature during evaluation of numerical integrals)
  double kT_ = DEFAULT_ENERGY_SCALE;
  /// The abbreviation of the intermolecular potential name
  PotentialName potential_id_;

  /// The underlying potential evaluated at r
  virtual double vr_(double r) const = 0;
  /// The underlying potential evaluated at r, computed from r^2
  virtual double vr2_(double r2) const = 0;
  /// Integral kernel for calculating Barker-Henderson hard-sphere diameter
  virtual double bh_diameter_kernel(double r) const;
  /// Integral kernel for calculating van der Waals parameter
  virtual double vdw_kernel(double r) const;

  //endregion

 public:
  //region Cttor:

  Potential() = default;
  Potential(double sigma, double epsilon, double r_cutoff);

  //endregion

  //region Inspectors:

  double sigma() const;
  double epsilon() const;
  double epsilon_shift() const;
  double r_cutoff() const;
  double r_min() const;
  double v_min() const;
  double r_attractive_min() const;
  double r_zero() const;
  bool bh_perturbation() const;
  double kT() const;
  const PotentialName& id() const;

  // endregion

  //region Methods:

  /// The repulsive part of the potential
  double w_repulsive(double r) const;

  /// The attractive tail
  double w_attractive(double r) const;

  /// The attractive part calculated from r2
  double w_attractive_r2(double r_squared) const;

  /// The cut and shifted potential at point r
  double v_potential(double r) const;

  /// The cut and shifted potential at point r calculated from r^2
  double v_potential_r2(double r_squared) const;

  /// Identifier string with the name of the potential and some characteristic parameters
  std::string identifier() const;

  /// Distance from origin that attractive continuate extends
  void SetWCALimit(double r);

  /// Use the Barker-Henderson split of the potential
  void SetBHPerturbation();

  /// Computes (and sets) the hard-core diameter of the original potential (e.g. zero for LJ)
  virtual double FindHardCoreDiameter() const = 0;

  /// Computes the position where the
  virtual double FindRMin() = 0;

  /// Computes the hard-sphere diameter by numerical integration
  double FindHardSphereDiameter(double kT);

  /// Compute the van der Waals pair-correlation-integral contribution to the free energy
  double ComputeVanDerWaalsIntegral(double kT);
  //endregion

  friend class boost::serialization::access;
  /**
   * @brief Allows the serialization of an object for saving it
   * @details For the saving of an object we need a serialization method which tells the program
   *    how (and what) to save from the object at hand. This is precisely what we obtain with
   *    the boost::serialization library. When the class `Archive` corresponds to an output archive,
   *    the `&` operator is defined similar the iostream's operator `<<`.  Likewise, when the class
   *    `Archive` is a type of input archive the `&` operator is defined similar to the
   *    counterpart `>>`.
   *
   * @tparam Archive The input/output archive object class
   * @param archive The archive-object reference
   * @param version The serialization library stores a version number in the archive for each class
   *    serialized. By default this version number is 0. When the archive is loaded, the version number
   *    under which it was saved is read.
   *    More info: https://www.boost.org/doc/libs/1_56_0/libs/serialization/doc/tutorial.html
   */
  template<class Archive>
  void serialize(Archive& archive, const unsigned int version)
  {
    archive & sigma_;
    archive & epsilon_;
    archive & epsilon_shift_;
    archive & r_cutoff_;
    archive & r_min_;
    archive & v_min_;
    archive & r_attractive_min_;
    archive & r_zero_;
    archive & bh_perturbation_;
    archive & kT_;
  }
};

}}}}

#endif  // CLASSICALDFT_POTENTIAL_H
