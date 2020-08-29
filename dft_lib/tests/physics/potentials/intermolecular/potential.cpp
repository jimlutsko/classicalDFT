#include <gtest/gtest.h>
#include <cmath>

#include "dft_lib/physics/potentials/intermolecular/potential.h"

using namespace dft_core::physics::potentials;
using namespace dft_core::physics::potentials::intermolecular;

//region Cttors:

/**
 * The class FakePotential is required to test the abstract class Potential.
 * This class is a trivial inheritor of the Potential class so that we can test the
 * underlying functionality of the mother class.
 */
class FakePotential : public intermolecular::Potential
{
 private:
  double vr_(double r) const override { return vr2_(r*r); }
  double vr2_(double r2) const override { return r2; }

 public:
  FakePotential(): intermolecular::Potential() {}
  FakePotential(double sigma, double epsilon, double r_cutoff): Potential(sigma, epsilon, r_cutoff) {}
  double FindHardCoreDiameter() const override { return 0; }
  double FindRMin() const override { return 0; }
};

TEST(intermolecular_potential, potential_cttor_works_ok)
{
  auto v_test = FakePotential();
  EXPECT_DOUBLE_EQ(v_test.sigma(), DEFAULT_LENGTH_SCALE);
  EXPECT_DOUBLE_EQ(v_test.epsilon(), DEFAULT_ENERGY_SCALE);
  EXPECT_DOUBLE_EQ(v_test.epsilon_shift(), DEFAULT_ZERO);
  EXPECT_DOUBLE_EQ(v_test.r_cutoff(), -DEFAULT_LENGTH_SCALE);
  EXPECT_DOUBLE_EQ(v_test.r_min(), DEFAULT_LENGTH_SCALE);
  EXPECT_DOUBLE_EQ(v_test.v_min(), -DEFAULT_ENERGY_SCALE);
  EXPECT_DOUBLE_EQ(v_test.r_attractive_min(), -DEFAULT_ZERO);
  EXPECT_DOUBLE_EQ(v_test.r_zero(), DEFAULT_LENGTH_SCALE);
  EXPECT_DOUBLE_EQ(v_test.bh_perturbation(), false);
  EXPECT_DOUBLE_EQ(v_test.kT(), DEFAULT_ENERGY_SCALE);

  EXPECT_DOUBLE_EQ(v_test.w_repulsive(0), DEFAULT_ENERGY_SCALE);
  EXPECT_DOUBLE_EQ(v_test.w_repulsive(0.5), 1.25);
  EXPECT_DOUBLE_EQ(v_test.w_attractive(0), -DEFAULT_ENERGY_SCALE);
  EXPECT_DOUBLE_EQ(v_test.w_attractive(0.5), -DEFAULT_ENERGY_SCALE);

  EXPECT_DOUBLE_EQ(v_test.v_potential(0.5), 0.25);

  v_test.SetBHPerturbation();
  EXPECT_DOUBLE_EQ(v_test.r_attractive_min(), DEFAULT_LENGTH_SCALE);
  EXPECT_DOUBLE_EQ(v_test.w_repulsive(0), DEFAULT_ZERO);
  EXPECT_DOUBLE_EQ(v_test.w_repulsive(0.5), 0.25);

  v_test.SetWCALimit(0.5);
  EXPECT_DOUBLE_EQ(v_test.r_attractive_min(), 0.5);

  auto d_hs = v_test.FindHardSphereDiameter(0);
  EXPECT_DOUBLE_EQ(d_hs, DEFAULT_LENGTH_SCALE);
  d_hs = v_test.FindHardSphereDiameter(0.1);
  EXPECT_DOUBLE_EQ(d_hs, 0.719752609493357);

  auto vdw = v_test.ComputeVanDerWaalsIntegral(0.5);
  EXPECT_DOUBLE_EQ(vdw, 0);
}

//endregion