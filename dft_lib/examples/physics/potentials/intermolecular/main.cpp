#include <cmath>
#include <armadillo>
#include "classical_dft"

auto conv_arma_to_vec(const arma::vec& x)
{
  auto y = arma::conv_to<std::vector<double>>::from(x);
  return y;
}

int main(int argc, char **argv)
{
  console::Info("Initialising Grace...");

  //region Cttor
  auto g = dft_core::grace_plot::Grace();
  const int N_POINTS = 200;
  //endregion

  using namespace dft_core::physics::potentials::intermolecular;
  auto lj = LennardJones();
  auto twf = tenWoldeFrenkel();

  auto x = arma::linspace(0.75, 1.5, N_POINTS);
  auto x_vector = conv_arma_to_vec(x);
  auto lj_vector = lj.v_potential(x_vector);
  auto twf_vector = twf.v_potential(x_vector);

  auto lj_ds = g.AddDataset(x_vector, lj_vector);
  auto twf_ds = g.AddDataset(x_vector, twf_vector);

  g.SetXLimits(x.min(), x.max());
  g.SetYLimits(-2, 10);

  auto lj_min = g.AddDataset(std::vector<double>{lj.r_min()}, std::vector<double>{lj.v_min()});
  g.SetSymbol(dft_core::grace_plot::Symbol::SQUARE, lj_min);
  g.SetSymbolFill(dft_core::grace_plot::Color::RED, lj_min);

  auto twf_min = g.AddDataset(std::vector<double>{twf.r_min()}, std::vector<double>{twf.v_min()});
  g.SetSymbol(dft_core::grace_plot::Symbol::DIAMOND, twf_min);
  g.SetSymbolFill(dft_core::grace_plot::Color::BLUE, twf_min);

  g.RedrawAndWait();

  /*
  //region Example of adding dataset
  console::Info("Adding new dataset");
  arma::vec z_vector = arma::cos(x_vector);

  auto dataset_id = g.AddDataset(
      arma::conv_to<std::vector<double>>::from(x_vector),
      arma::conv_to<std::vector<double>>::from(z_vector)
  );
  console::Debug("The new dataset id = " + std::to_string(dataset_id));
  g.SetColor(dft_core::grace_plot::Color::RED, dataset_id);
  g.RedrawAndWait();
  //endregion

  //region Example of replacing dataset
  console::Info("Example: Replacing existing dataset");
  arma::vec w_vector = arma::tan(x_vector);

  auto existing_dataset_id = g.AddDataset(
      arma::conv_to<std::vector<double>>::from(x_vector),
      arma::conv_to<std::vector<double>>::from(w_vector)
  );
  console::Debug("The new dataset id = " + std::to_string(dataset_id));
  g.SetColor(dft_core::grace_plot::Color::BLUE, existing_dataset_id);
  g.RedrawAndWait();

  console::Warning("Replacing dataset id: " + std::to_string(existing_dataset_id));
  w_vector = arma::sqrt(x_vector)/sqrt(2*M_PI);
  g.ReplaceDataset(
      arma::conv_to<std::vector<double>>::from(x_vector),
      arma::conv_to<std::vector<double>>::from(w_vector),
      existing_dataset_id
  );
  console::Debug("The replaced graph id = " + std::to_string(existing_dataset_id));

  g.SetLegend("x\\S1/2\\N", existing_dataset_id);
  //g.SetColor(dft_core::grace_plot::MAGENTA, existing_dataset_id);
  g.RedrawAndWait();
  //endregion

  //region Example of axis labels
  console::Info("Example: Setting the X and Y labels");
  g.SetLabel("This is X", dft_core::grace_plot::Axis::X);
  g.SetLabel("This is Y", dft_core::grace_plot::Axis::Y);
  g.RedrawAndWait();
  //endregion

  //region Example of graph title
  console::Info("Example: Setting the Title and Subtitle");
  g.SetTitle("This is the title");
  g.SetSubtitle("And this is the subtitle");
  g.RedrawAndWait();
  //endregion

  //region Example of setting limits
  console::Info("Example: Setting the limits");
  g.SetLimits(std::vector<double>{ -0.1, 2*M_PI+0.1 }, std::vector<double>{-1.2, 1.2});
  g.RedrawAndWait();
  //endregion

  //region Example of setting ticks
  console::Info("Example: Setting the tick spacing");
  g.SetTicks(0.5, 0.1);
  g.RedrawAndWait(false, false);
  //endregion

  //region Example of setting line type
  console::Info("Example: Setting the line type");
  g.SetLineType(dft_core::grace_plot::LineType::NO_LINE, 0);
  g.RedrawAndWait(false, false);
  g.SetLineType(dft_core::grace_plot::LineType::LINE, 0);
  g.RedrawAndWait(false, false);
  g.SetLineType(dft_core::grace_plot::LineType::DOTTEDLINE, 0);
  g.RedrawAndWait(false, false);
  g.SetLineType(dft_core::grace_plot::LineType::DASHEDLINE_EN, 0);
  g.RedrawAndWait(false, false);
  //endregion

  //region Example of setting symbols and fills
  console::Info("Example: Setting symbol and symbol color");
  g.SetSymbol(dft_core::grace_plot::Symbol::TRIANGLE_DOWN, 0);
  g.SetSymbolColor(dft_core::grace_plot::Color::BLUE, 0);
  g.SetSymbolFill(dft_core::grace_plot::Color::BLUE, 0, 0, 4);

  g.SetSymbol(dft_core::grace_plot::Symbol::TRIANGLE_LEFT, 1);
  g.SetSymbolColor(dft_core::grace_plot::Color::DARKGREEN, 1);

  g.SetSymbol(dft_core::grace_plot::Symbol::DIAMOND, 2);
  g.SetSymbolFill(dft_core::grace_plot::Color::RED, 2);
  g.RedrawAndWait(false, false);
  //endregion

  //region Example of setting symbol size
  console::Info("Example: Setting symbol size");
  g.SetSymbolSize(0.5,0);
  g.SetSymbolSize(1.5,1);
  g.RedrawAndWait(false, false);
  //endregion

  //region Example of Export
  console::Info("Example: Saving the result as PNG");
  g.PrintToFile("test_graph.png", dft_core::grace_plot::ExportFormat::PNG);
  //g.PrintToFile("test_graph.pdf", dft_core::grace_plot::ExportFormat::PDF);
  //g.PrintToFile("test_graph.jpg", dft_core::grace_plot::ExportFormat::JPG);
  //g.PrintToFile("test_graph.eps", dft_core::grace_plot::ExportFormat::EPS1);
  //g.PrintToFile("test_graph.ps", dft_core::grace_plot::ExportFormat::PS);
  g.RedrawAndWait(false, false);
  //endregion
  */
}