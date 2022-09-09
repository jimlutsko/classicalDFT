# Compensated Summation

### Introduction

Very likely, one of the most common operations we need to carry out whilst developing a numerical project is the `add` method, typically used via the arithmetic operators `+` or `+=`  (equivalently `-` or `-=`). However, adding a [sequence](https://en.wikipedia.org/wiki/Sequence) of finite-[precision](https://en.wikipedia.org/wiki/Decimal_precision) [floating-point numbers](https://en.wikipedia.org/wiki/Floating-point_number) is not as trivial as it might seem due to the [numerical error](https://en.wikipedia.org/wiki/Numerical_error) induced by every member of the sequence under consideration. Indeed, the trivial summation of a sequence of *n* numbers yields (worst-case scenario) an error which grows linearly with the number of elements of the sequence. This obviously becomes a major problem when we need to sum a very large sequence of floating-point numbers and the problem at hand is very sensitive to errors. The main goal (and characteristic) of **compensated summation** is to achieve a worst-case error bound effectively independent of *n*, so a large number of values can be summed with an error that only depends on the floating-point precision. The `classicalDFT` library comes with an implementation of the **Kahan summation algorithm** (and its most important variants).  Such functionality is found under the namespace `dft_core::numerics::arithmetic::sumation`, which provides us with three free methods:

* `KahanBabuskaSum`: The standard *compensated-summation* algorithm attributed to Kahan and Babuska
* `KahanBabuskaNeumaierSum`: A improvement of the compensated-summation algorithm by Neumaier
* `KahanBabuskaKleinSum`: A further improvement of the algorithm by Klein

Besides, we can finde an interface class `CompensatedSum` which allows us to use them effortlessly.

### Examples

The best way of showing the convenience offered by `arithmetic::summation` is by example. Thus, we are going to proceed by inserting the code in [`main.cpp`](main.cpp):

```c++
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
```

After compilation and running we will get the following results:

```bash
Test 1: Sum from scratch
Trivial sum [1]: 10.5
Kahan-Babuska sum [1]:10.5
Kahan-Babuska err [1] = -8.88178e-16
Kahan-Babuska-Neumaier sum [1]: 10.5
Kahan-Babuska-Neumaier err [1]: 1.11022e-15
Kahan-Babuska-Klein sum [1]: 10.5
Kahan-Babuska-Klein err [1/a]: 1.11022e-15
Kahan-Babuska-Klein err [1/b]: 0

Test 2: Continuation of previous sum
Trivial sum [2]: 21.0014
Kahan-Babuska sum [2]: 21.0014
Kahan-Babuska err [2]: -1.77636e-15
Kahan-Babuska-Neumaier sum [2]: 21.0014
Kahan-Babuska-Neumaier err [2]: 3.55271e-15
Kahan-Babuska-Klein sum [1]: 21.0014
Kahan-Babuska-Klein err [1/a]: 3.55271e-15
Kahan-Babuska-Klein err [1/b]: 0

Test 3: Using CompensatedSum class
Testing the operator '+=' [1] = 10.5
Testing the operator '+=' [2] = 21.0014
Press enter to continue...
```

