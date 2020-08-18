#include <gtest/gtest.h>
#include <cmath>

#include "dft_lib/physics/potentials/intermolecular/potential.h"

using namespace dft_core::physics::potentials;
using namespace dft_core::physics::potentials::intermolecular;

//region Cttors:

/**
 * The class fake potential is required to test the abstract class Potential.
 * This class is a trivial inheritor of the Potential class so that we can test the
 * underlying functionality of the mother class.
 */
class fake_potential: public intermolecular::Potential
{
 public:
  fake_potential(): intermolecular::Potential() {}
  fake_potential(double sigma, double epsilon, double r_cutoff): Potential(sigma, epsilon, r_cutoff) {}
  double HardCoreDiameter() const { return 0; }
};

TEST(intermolecular_potential, potential_cttor_works_ok)
{
  auto v_test = fake_potential();
  EXPECT_DOUBLE_EQ(v_test.sigma(), DEFAULT_LENGTH_SCALE);
  EXPECT_DOUBLE_EQ(v_test.epsilon(), DEFAULT_ENERGY_SCALE);
  EXPECT_DOUBLE_EQ(v_test.epsilon_shift(), DEFAULT_ZERO);
  EXPECT_DOUBLE_EQ(v_test.r_cutoff(), -DEFAULT_LENGTH_SCALE);
  EXPECT_DOUBLE_EQ(v_test.r_min(), DEFAULT_LENGTH_SCALE);
  EXPECT_DOUBLE_EQ(v_test.v_min(), -DEFAULT_ENERGY_SCALE);
  EXPECT_DOUBLE_EQ(v_test.r_attractive_min(), -DEFAULT_ZERO);
  EXPECT_DOUBLE_EQ(v_test.r_zero(), DEFAULT_LENGTH_SCALE);
  EXPECT_EQ(v_test.bh_perturbation(), false);
  EXPECT_DOUBLE_EQ(v_test.kT(), DEFAULT_ENERGY_SCALE);
}

//endregion