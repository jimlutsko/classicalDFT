#include <math.h>
#include <armadillo>
#include "classical_dft"

int main(int argc, char **argv)
{
  console::info("Initialising Grace...");

  //region Cttor
  auto g = dft_core::grace_plot::Grace();
  const int N_POINTS = 100;
  //endregion

  //region Example of adding points to the default dataset 0
  console::info("Plotting sin(x): x in [0, 2*PI]");

  auto x_vector = arma::linspace(0, 2*M_PI, N_POINTS);
  auto y_vector = arma::sin(x_vector);

  for (auto k = 0; k < x_vector.size(); k++)
  { g.AddPoint(x_vector[k], y_vector[k]); }

  g.SetXLimits(x_vector.min(), x_vector.max());
  g.SetYLimits(y_vector.min(), y_vector.max());

  g.Redraw();
  console::wait();
  //endregion

  //region Example of adding dataset
  console::info("Adding new dataset");
  arma::vec z_vector = arma::cos(x_vector);

  auto graph_id = g.AddDataset(
      arma::conv_to<std::vector<double>>::from(x_vector),
      arma::conv_to<std::vector<double>>::from(z_vector)
      );

  g.Redraw();
  console::debug("The graph id = " + std::to_string(graph_id));
  console::wait();
  //endregion
}