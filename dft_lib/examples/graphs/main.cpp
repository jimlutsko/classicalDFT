#include "classical_dft"

int main(int argc, char **argv)
{
  console::info("Initialising Grace...");
  console::wait();

  auto g = dft_core::grace_plot::Grace();
  g.AddPoint(1,1);
  g.AddPoint(2,1);
  g.AddPoint(3,2);

  console::wait();
}