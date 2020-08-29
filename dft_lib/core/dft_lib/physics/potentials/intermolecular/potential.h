#ifndef CLASSICALDFT_POTENTIAL_H
#define CLASSICALDFT_POTENTIAL_H

#include <iostream>
#include <vector>

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
const double DEFAULT_ALPHA_PARAMETER = 50.0;
const double DEFAULT_ZERO = 0.0;
const double MAX_POTENTIAL_VALUE = 1e50;

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
 *
 * @details The potential class brings the ability to split the potential contribution into two parts:
 *      - hard sphere (purely repulsive) and,
 *      - attractive,
 *      and, also comes with the functionality to compute the hard-sphere diameter.
 *
 *      By default construction, the potential is written as the sum: w_repulsive(r) + w_attractive(r)
 *      with w_repulsive(r) = 0, for r > r_minimum_ (where r_minimum_ is where the minimum of the
 *      potential vr_(r) is located at); and w_attractive(r) = v(r_minimum_), for r < r_minimum_,
 *      i.e. the minimum of the potential.
 *
 *      If r_attractive_min_ > 0, then w_attractive(r) is zero for r < r_attractive_min_,
 *      i.e. the continuation inside the core is truncated.
 *
 *      If the bh_perturbation_ (Barker-Henderson model) flag is true, the potential is split at the
 *      point it reaches zero and the attractive part is only non-zero for positions greater than
 *      this point.
 */
class Potential
{
 protected:
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

  //endregion

  //region Methods:

  /// The underlying potential evaluated at r
  virtual double vr_(double r) const = 0;
  /// The underlying potential evaluated at r, computed from r^2
  virtual double vr2_(double r2) const = 0;
  /// Integral kernel for calculating Barker-Henderson hard-sphere diameter
  virtual double bh_diameter_kernel(double r) const;
  /// Integral kernel for calculating van der Waals parameter
  virtual double vdw_kernel(double r) const;

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

  //endregion

 public:
  //region Cttor:
  /**
   * @brief Default constructor of the class. It comes with the private parameters initialised
   *    with the default energy or length scale: DEFAULT_ENERGY_SCALE and DEFAULT_LENGTH_SCALE,
   *    respectively
   */
  Potential() = default;
  /**
   * @brief Constructor used for the parameterization of a Potential object (sigma, epsilon, r_cutoff)
   * @param sigma The typical length scale defining the problem at hand
   * @param epsilon The typical energy scale defining the problem at hand
   * @param r_cutoff The distance at which the potential energy is considered negligible, hence used
   *        for truncation purposes. E.g., if `r_cutoff = 2.5` the potential will be set to zero
   *        from `r=r_cutoff` onwards, which is equivalent to shift the potential by
   *        `epsilon_shift = v(r_cutoff)`
   */
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
  std::vector<double> v_potential(const std::vector<double>& r) const;

  /// The cut and shifted potential at point r calculated from r^2
  double v_potential_r2(double r_squared) const;
  std::vector<double> v_potential_r2(const std::vector<double>& r_squared) const;

  /// Identifier string with the name of the potential and some characteristic parameters
  std::string identifier() const;

  /// Distance from origin that attractive continuate extends
  void SetWCALimit(double r);

  /// Use the Barker-Henderson split of the potential
  void SetBHPerturbation();

  /// Computes the hard-core diameter of the potential (e.g., zero for LJ). The hard-core diameter
  /// determines the integration lower limit when computing the hard-sphere diameter.
  virtual double FindHardCoreDiameter() const = 0;
  /// Computes the position where the potential reaches its minimum
  virtual double FindRMin() const = 0;

  /// Computes the hard-sphere diameter by numerical integration
  double FindHardSphereDiameter(double kT);

  /// Compute the van der Waals pair-correlation-integral contribution to the free energy
  double ComputeVanDerWaalsIntegral(double kT);
  //endregion
};

/**
 * @brief Lennard-Jones potential (also termed 12-6 potential)
 * @details The `LennardJones` class embodies the Lennard-Jones potential of interaction.
 *      This is a mathematically simple model that approximates the intermolecular potential energy
 *      between a pair of neutral atoms or molecules as metals or cyclic alkanes.
 *      More information about the LJ[12-6] potential can be found at: https://en.wikipedia.org/wiki/Lennard-Jones_potential
 */
class LennardJones final: public Potential
{
 protected:
  //region Methods:

  /// The underlying potential evaluated at r
  double vr_(double r) const override;
  /// The underlying potential evaluated at r, computed from r^2
  double vr2_(double r2) const override;

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive& archive, const unsigned int version)
  {
    archive & boost::serialization::base_object<Potential>(*this);
    boost::serialization::void_cast_register<LennardJones, Potential>(static_cast<LennardJones*>(nullptr),static_cast<Potential*>(nullptr));
  }

  //endregion

 public:
  //region Cttors:

  /**
   * @brief Default constructor of the class. It comes with the private parameters initialised
   *    with the default energy or length scale: DEFAULT_ENERGY_SCALE and DEFAULT_LENGTH_SCALE,
   *    respectively
   */
  LennardJones();
  /**
   * @brief Constructor used for the parameterization of a LennardJones object (sigma, epsilon, r_cutoff)
   * @param sigma The typical length scale defining the problem at hand
   * @param epsilon The typical energy scale defining the problem at hand
   * @param r_cutoff The distance at which the potential energy is considered negligible, hence used
   *        for truncation purposes. E.g., if `r_cutoff = 2.5` the potential will be set to zero
   *        from `r=r_cutoff` onwards, which is equivalent to shift the potential by
   *        `epsilon_shift = v(r_cutoff)`
   */
  LennardJones(double sigma, double epsilon, double r_cutoff);

  //endregion

  //region Methods:

  double FindHardCoreDiameter() const override;
  double FindRMin() const override;

  //endregion
};

/**
 * @brief ten Wolde-Frenkel potential
 * @details The `tenWoldeFrenkel` class embodies the potential of interaction introduced in the
 *          work of ten Wolde and Frenkel (https://science.sciencemag.org/content/277/5334/1975)
 *          for the simulation of globular proteins with short-range attractive interactions.
 */
class tenWoldeFrenkel final: public Potential
{
 private:
  //region Attributes:

  /// The parameter alpha modulating the potential shape
  double alpha_ = DEFAULT_ALPHA_PARAMETER;

  //endregion

 protected:
  //region Methods:

  /// The underlying potential evaluated at r
  double vr_(double r) const override;
  /// The underlying potential evaluated at r, computed from r^2
  double vr2_(double r2) const override;

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive& archive, const unsigned int version)
  {
    archive & boost::serialization::base_object<Potential>(*this);
    archive & alpha_;
    boost::serialization::void_cast_register<tenWoldeFrenkel, Potential>(static_cast<tenWoldeFrenkel*>(nullptr),static_cast<Potential*>(nullptr));
  }

  //endregion

 public:
  //region Cttors:

  /**
   * @brief Default constructor of the class. It comes with the private parameters initialised
   *    with the default energy or length scale: DEFAULT_ENERGY_SCALE and DEFAULT_LENGTH_SCALE,
   *    respectively
   */
  tenWoldeFrenkel();
  /**
   * @brief Constructor used for the parameterization of a LennardJones object (sigma, epsilon, r_cutoff)
   * @param sigma The typical length scale defining the problem at hand
   * @param epsilon The typical energy scale defining the problem at hand
   * @param r_cutoff The distance at which the potential energy is considered negligible, hence used
   *        for truncation purposes. E.g., if `r_cutoff = 2.5` the potential will be set to zero
   *        from `r=r_cutoff` onwards, which is equivalent to shift the potential by
   *        `epsilon_shift = v(r_cutoff)`
   * @param alpha The `alpha` parameter in the tWF potential (default = 50)
   */
  tenWoldeFrenkel(double sigma, double epsilon, double r_cutoff, double alpha = DEFAULT_ALPHA_PARAMETER);

  //endregion

  //region Inspectors:

  /// Gets the value
  double alpha() const;

  //endregion

  //region Methods:

  double FindHardCoreDiameter() const override;
  double FindRMin() const override;

  //endregion
};

}}}}

#endif  // CLASSICALDFT_POTENTIAL_H
