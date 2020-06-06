#include "classical_dft"

int main()
{
  using namespace dft_core::numerics::arithmetic;

  auto x_input_1 = std::vector<double>{1.0 + 1E-14, 2.5 + 1E-14, 3.0 + 1E-14, 4.0 + 1E-14};
  auto x_input_2 = std::vector<double>{1.00100001, 2.50010002, 3.00020001, 4.00010003};

  auto trivial_sum = summation::StandardVectorSum<>(x_input_1);
  double kb_sum, err;
  std::tie(kb_sum, err) = summation::KahanBabuskaSum(x_input_1);

  console::WriteLine("Trivial sum [1]:");
  console::WriteLine(trivial_sum);
  console::WriteLine("Kahan-Babuska sum [1]:");
  console::WriteLine(kb_sum);
  console::WriteLine("Kahan-Babuska err [1] = ");
  console::WriteLine(err);

  console::NewLine();
  trivial_sum += summation::StandardVectorSum(x_input_2);
  std::tie(kb_sum, err) = summation::KahanBabuskaSum(x_input_2, kb_sum, err);
  console::WriteLine("Trivial sum [2]: ");
  console::WriteLine(trivial_sum);
  console::WriteLine("Kahan-Babuska sum [2]: ");
  console::WriteLine(kb_sum);
  console::WriteLine("Kahan-Babuska err [2]: ");
  console::WriteLine(err);

  console::Wait();
}