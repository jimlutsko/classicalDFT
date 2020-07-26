#include "classical_dft"

int main()
{
  using namespace dft_core::numerics::arithmetic;

  // region set-up:
  auto x_input_1 = std::vector<double>{1.0 + 1E-14, 2.5 + 1E-14, 3.0 + 1E-14, 4.0 + 1E-14};
  auto x_input_2 = std::vector<double>{1.00100001, 2.50010002, 3.00020001, 4.00010003};

  auto trivial_sum = summation::StandardVectorSum(x_input_1);
  console::WriteLine("Test 1: Sum from scratch");
  console::Write("Trivial sum [1]: "); console::WriteLine(trivial_sum);
  // endregion

  // region free-methods:

  double kb_sum; std::vector<double> kb_err;
  std::tie(kb_sum, kb_err) = summation::KahanBabuskaSum(x_input_1);
  console::Write("Kahan-Babuska sum [1]:"); console::WriteLine(kb_sum);
  console::Write("Kahan-Babuska err [1] = "); console::WriteLine(kb_err.front());

  double neum_sum; std::vector<double> neum_err;
  std::tie(neum_sum, neum_err) = summation::KahanBabuskaNeumaierSum(x_input_1);
  console::Write("Kahan-Babuska-Neumaier sum [1]: "); console::WriteLine(neum_sum);
  console::Write("Kahan-Babuska-Neumaier err [1]: "); console::WriteLine(neum_err.front());

  double klein_sum; std::vector<double> klein_err;
  std::tie(klein_sum, klein_err) = summation::KahanBabuskaKleinSum(x_input_1);
  console::Write("Kahan-Babuska-Klein sum [1]: "); console::Write(klein_sum);
  console::Write("Kahan-Babuska-Klein err [1/a]: "); console::WriteLine(klein_err[0]);
  console::Write("Kahan-Babuska-Klein err [1/b]: "); console::WriteLine(klein_err[1]);

  console::NewLine();
  trivial_sum += summation::StandardVectorSum(x_input_2);
  std::tie(kb_sum, kb_err) = summation::KahanBabuskaSum(x_input_2, kb_sum, kb_err);
  std::tie(neum_sum, neum_err) = summation::KahanBabuskaNeumaierSum(x_input_2, neum_sum, neum_err);
  std::tie(klein_sum, klein_err) = summation::KahanBabuskaKleinSum(x_input_2, klein_sum, klein_err);

  console::WriteLine("Test 2: Continuation of previous sum");
  console::Write("Trivial sum [2]: "); console::WriteLine(trivial_sum);

  console::Write("Kahan-Babuska sum [2]: "); console::WriteLine(kb_sum);
  console::Write("Kahan-Babuska err [2]: "); console::WriteLine(kb_err.front());

  console::Write("Kahan-Babuska-Neumaier sum [2]: "); console::WriteLine(neum_sum);
  console::Write("Kahan-Babuska-Neumaier err [2]: "); console::WriteLine(neum_err.front());

  console::Write("Kahan-Babuska-Klein sum [1]: "); console::WriteLine(klein_sum);
  console::Write("Kahan-Babuska-Klein err [1/a]: "); console::WriteLine(klein_err[0]);
  console::Write("Kahan-Babuska-Klein err [1/b]: "); console::WriteLine(klein_err[1]);

  console::NewLine();
  // endregion

  // region class:
  console::WriteLine("Test 3: Using CompensatedSum class");
  summation::CompensatedSum x_cs;

  x_cs += x_input_1;
  console::Write("Testing the operator '+=' [1] = "); console::WriteLine(x_cs);

  x_cs += x_input_2;
  console::Write("Testing the operator '+=' [2] = "); console::WriteLine(x_cs);

  console::Wait();
}