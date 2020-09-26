#include <gtest/gtest.h>
#include <cmath>

#include "dft_lib/physics/potentials/intermolecular/potential.h"

#include <boost/range/combine.hpp>

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
  EXPECT_DOUBLE_EQ(v_test.r_cutoff(), DEFAULT_CUTOFF);
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

TEST(intermolecular_potential, potential_brackets_works)
{
  auto lj = LennardJones();
  auto d_expected = lj.v_potential(1);
  auto d_actual = lj(1);
  EXPECT_DOUBLE_EQ(d_expected, d_actual);

  auto x = std::vector<double>{ 1, 2 };
  auto vec_expected = lj.v_potential(x);
  auto vec_actual = lj(x);

  for (auto tup : boost::combine(vec_expected, vec_actual))
  {
    double x, y; boost::tie(x,y) = tup;
    EXPECT_DOUBLE_EQ(x, y);
  }

  auto r = arma::linspace(1, 5, 10);
  auto arma_expected = lj.v_potential(r);
  auto arma_actual = lj(r);

  for (auto tup : boost::combine(arma_expected, arma_actual))
  {
    double x, y; boost::tie(x,y) = tup;
    EXPECT_DOUBLE_EQ(x, y);
  }
}

TEST(intermolecular_potential, potential_attractive_part_ok)
{
  auto twf = tenWoldeFrenkel();
  auto v_actual = twf.w_attractive(0);
  EXPECT_DOUBLE_EQ(twf.v_min(), v_actual);

  v_actual = twf.w_attractive(twf.r_min());
  EXPECT_DOUBLE_EQ(twf.v_min(), v_actual);

  v_actual = twf.w_attractive(2 * twf.r_min());
  EXPECT_DOUBLE_EQ(twf(2 * twf.r_min()), v_actual);
}

TEST(intermolecular_potential, potential_attractive_part_bh_ok)
{
  auto twf = tenWoldeFrenkel();
  twf.SetBHPerturbation();

  auto v_actual = twf.w_attractive(0.0);
  EXPECT_DOUBLE_EQ(0.0, v_actual);

  v_actual = twf.w_attractive(0.999 * twf.r_zero());
  EXPECT_DOUBLE_EQ(0.0, v_actual);

  v_actual = twf.w_attractive(1.001 * twf.r_zero());
  EXPECT_DOUBLE_EQ(twf(1.001 * twf.r_zero()), v_actual);

  v_actual = twf.w_attractive(twf.r_min());
  EXPECT_DOUBLE_EQ(twf.v_min(), v_actual);

  v_actual = twf.w_attractive(2 * twf.r_min());
  EXPECT_DOUBLE_EQ(twf(2 * twf.r_min()), v_actual);
}

TEST(intermolecular_potential, potential_repulsive_part_ok)
{
  auto lj = LennardJones();
  auto r = 0.9 * lj.r_min();

  auto actual = lj.w_repulsive(r);
  auto expected = lj(r) - lj.v_min();
  EXPECT_DOUBLE_EQ(expected, actual);

  r = 1.1 * lj.r_min();
  actual = lj.w_repulsive(r);
  expected = 0.0;
  EXPECT_DOUBLE_EQ(expected, actual);

  lj.SetBHPerturbation();

  r = 0.999 * lj.r_zero();
  actual = lj.w_repulsive(r);
  expected = lj(r);
  EXPECT_DOUBLE_EQ(expected, actual);

  r = 1.001 * lj.r_zero();
  actual = lj.w_repulsive(r);
  expected = 0.0;
  EXPECT_DOUBLE_EQ(expected, actual);
}

TEST(intermolecular_potential, WCA_split_works_ok)
{
  auto lj = LennardJones();
  auto x = arma::linspace(0.5, 5, 100);
  auto expected = lj(x);
  arma::vec actual = lj.w_attractive(x) + lj.w_repulsive(x);

  for (auto tup : boost::combine(expected, actual))
  {
    double x, y; boost::tie(x,y) = tup;
    EXPECT_NEAR(x, y, 1e-6);
  }
}

TEST(intermolecular_potential, WCA_split_works_with_BH_ok)
{
  auto lj = LennardJones();
  lj.SetBHPerturbation();

  auto x = arma::linspace(0.5, 5, 100);
  auto expected = lj(x);
  arma::vec actual = lj.w_attractive(x) + lj.w_repulsive(x);

  for (auto tup : boost::combine(expected, actual))
  {
    double x, y; boost::tie(x,y) = tup;
    EXPECT_NEAR(x, y, 1e-6);
  }
}
//endregion