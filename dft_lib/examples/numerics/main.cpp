#include "classical_dft"

class LocalProblem
{
 private:
  double param_ = 0.1;

 public:
  explicit LocalProblem(double param = 0.1): param_(param) {}
  double Function(double x) const { return param_ * exp(-x); }
};

int main()
{
  using namespace dft_core;

  auto problem = LocalProblem(1.0);
  auto integrator = numerics::Integrator<LocalProblem>(problem, &LocalProblem::Function);

  /*
  std::cout << integrator.function(M_PI_2) << std::endl;
  std::cout << numerics::Integrator<LocalProblem>::integrand_function(M_PI_2, &integrator) << std::endl;
  */

  auto result = integrator.DefiniteIntegral(0, -log(0.5));
  std::cout << "result [0, -log(0.5)]: " << integrator.numerical_result() << std::endl;
  std::cout << "error: " << integrator.numerical_error() << std::endl;

  auto result_fast = integrator.DefiniteIntegralFast(0, -log(0.5));
  std::cout << "result (fast) [0, -log(0.5)]: " << integrator.numerical_result() << std::endl;
  std::cout << "error: " <<  integrator.numerical_error() << std::endl;

  auto result_semi_inf = integrator.UpperSemiInfiniteIntegral(0);
  std::cout << "result [0, +inf)]: " << integrator.numerical_result() << std::endl;
  std::cout << "error: " <<  integrator.numerical_error() << std::endl;
}