#include "classical_dft"

namespace example
{
class TestProblem
{
 private:
  double param_ = 0.1;

 public:
  explicit TestProblem(double param = 0.1) : param_(param) {}
  double NegativeExp(double x) const { return param_ * exp(-x); }
  double PositiveExp(double x) const { return param_ * exp(x); }
  double NormalDist(double x) const { return param_ * exp(-x * x * 0.5) / sqrt(2 * M_PI); }
};

double NegativeExpFunction(double x, const std::vector<double>& param) { return param[0] * param[1] * exp(-x); }
}

int main()
{
  using namespace dft_core::numerics::integration;

  auto problem = example::TestProblem(1.0);
  auto integrator = Integrator<example::TestProblem>(problem, &example::TestProblem::NegativeExp);

  console::WriteLine(console::format::Blink("Testing passing object->method:"));
  console::WriteLine("integrator.f(2*M_PI) = " + std::to_string(integrator.function(M_PI_2)));
  console::WriteLine("Integrator<LocalProblem>::integrand_function(2*M_PI, &integrator) = "
                     + std::to_string(Integrator<example::TestProblem>::integrand_function(M_PI_2, &integrator)));

  auto result = integrator.DefiniteIntegral(0, -log(0.5));
  console::NewLine();
  console::WriteLine(console::format::Blink("Testing integration methods: QAGS"));
  console::WriteLine("int[0, -log(0.5)] exp(-x) dx = " + std::to_string(integrator.numerical_result()));
  console::WriteLine("numerical_error =  " + std::to_string(integrator.numerical_error()));

  auto result_fast = integrator.DefiniteIntegralFast(0, -log(0.5));
  console::NewLine();
  console::WriteLine(console::format::Blink("Testing integration methods: QNG"));
  console::WriteLine("int[0, -log(0.5)] exp(-x) dx = " + std::to_string(integrator.numerical_result()));
  console::WriteLine("numerical_error =  " + std::to_string(integrator.numerical_error()));

  auto result_semi_inf = integrator.UpperSemiInfiniteIntegral(0);
  console::NewLine();
  console::WriteLine(console::format::Blink("Testing integration methods: QAGIU"));
  console::WriteLine("int[0, +inf] exp(-x) dx = " + std::to_string(integrator.numerical_result()));
  console::WriteLine("numerical_error =  " + std::to_string(integrator.numerical_error()));

  auto integrator_neg = Integrator<example::TestProblem>(problem, &example::TestProblem::PositiveExp);
  auto result_semi_neg = integrator_neg.LowerSemiInfiniteIntegral(0);
  console::NewLine();
  console::WriteLine(console::format::Blink("Testing integration methods: QAGIL"));
  console::WriteLine("int[-inf, 0] exp(x) dx = " + std::to_string(integrator.numerical_result()));
  console::WriteLine("numerical_error =  " + std::to_string(integrator.numerical_error()));

  auto integrator_gauss = Integrator<example::TestProblem>(problem, &example::TestProblem::NormalDist);
  auto result_gauss = integrator_gauss.FullInfiniteIntegral();
  console::NewLine();
  console::WriteLine(console::format::Blink("Testing integration methods: QAGI"));
  console::WriteLine("int[-inf, +inf] normal(x) dx = " + std::to_string(integrator.numerical_result()));
  console::WriteLine("numerical_error =  " + std::to_string(integrator.numerical_error()));

  auto parameters = std::vector<double>{ 0.5, 2 };
  auto func_integrator = FunctionIntegrator<std::vector<double>>(&example::NegativeExpFunction, parameters);

  result = func_integrator.DefiniteIntegral(0, -log(0.5));
  console::NewLine();
  console::WriteLine(console::format::Blink("Testing function-integration vector parameters: p = [0.5, 2]"));
  console::WriteLine("int[0, -log(0.5)] p[0] * p[1] * exp(-x) dx = " + std::to_string(func_integrator.numerical_result()));
  console::WriteLine("numerical_error =  " + std::to_string(func_integrator.numerical_error()));

  console::Wait();
}