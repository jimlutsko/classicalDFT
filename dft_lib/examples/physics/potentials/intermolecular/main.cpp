#include "classical_dft"
#include <armadillo>

/// A convenient wrapper to convert arma::vec -> std::vector
auto conv_arma_to_vec(const arma::vec& x)
{
  auto y = arma::conv_to<std::vector<double>>::from(x);
  return y;
}

int main(int argc, char **argv)
{
  console::Info("Initialising Grace...");

  //region Grace set up:
  auto g = dft_core::grace_plot::Grace();
  const int N_POINTS = 80;

  //region Grid set up:
  auto x_vector = arma::linspace(0.75, 1.5, N_POINTS);
  auto y_lims = std::vector<double>{-2, 10};
  g.SetXLimits(x_vector.min(), x_vector.max());
  g.SetYLimits(y_lims[0], y_lims[1]);
  //endregion

  //endregion

  //region Instantiation of the potentials:
  using namespace dft_core::physics::potentials::intermolecular;
  auto lj = LennardJones();
  auto twf = tenWoldeFrenkel();
  auto wrdf = WangRamirezDobnikarFrenkel();
  //endregion

  //region Lennard-Jones:
  //region Potential:
  auto lj_vector = lj.v_potential(x_vector); // equivalent to lj(x_vector);
  auto lj_ds = g.AddDataset(conv_arma_to_vec(x_vector), conv_arma_to_vec(lj_vector));
  //endregion

  //region Minimum:
  auto lj_min = g.AddDataset(std::vector<double>{lj.r_min()}, std::vector<double>{lj.v_min()});
  g.SetLineType(dft_core::grace_plot::LineStyle::NO_LINE, lj_min);
  g.SetSymbol(dft_core::grace_plot::Symbol::SQUARE, lj_min);
  g.SetSymbolFill(dft_core::grace_plot::Color::RED, lj_min);
  //endregion

  //region Hard-sphere diameter:
  auto y_vec = arma::linspace(y_lims[0], y_lims[1], 10);
  auto hs_diameter = lj.FindHardSphereDiameter(1.0);
  auto hs_x = arma::vec(10, arma::fill::ones); hs_x *= hs_diameter;
  auto lj_hs = g.AddDataset(conv_arma_to_vec(hs_x), conv_arma_to_vec(y_vec));
  g.SetColor(dft_core::grace_plot::Color::RED, lj_hs);
  g.SetLineType(dft_core::grace_plot::LineStyle::DASHEDLINE_EN, lj_hs);
  console::Info("LJ hard-sphere diameter (kT = 1.0) = " + std::to_string(hs_diameter));
  //endregion

  //region Pertubation theory:
  auto lj_att = lj.w_attractive(x_vector);
  auto lj_att_ds = g.AddDataset(conv_arma_to_vec(x_vector), conv_arma_to_vec(lj_att));
  g.SetColor(dft_core::grace_plot::Color::RED, lj_att_ds);
  g.SetLineType(dft_core::grace_plot::LineStyle::D_DOTTEDDASHEDLINE_EM, lj_att_ds);
  g.SetSymbol(dft_core::grace_plot::Symbol::CIRCLE, lj_att_ds);
  g.SetSymbolSize(0.25, lj_att_ds);

  auto lj_rep = lj.w_repulsive(x_vector);
  auto lj_rep_ds = g.AddDataset(conv_arma_to_vec(x_vector), conv_arma_to_vec(lj_rep));
  g.SetColor(dft_core::grace_plot::Color::RED, lj_rep_ds);
  g.SetLineType(dft_core::grace_plot::LineStyle::D_DOTTEDDASHEDLINE_EN, lj_rep_ds);
  g.SetSymbol(dft_core::grace_plot::Symbol::STAR, lj_rep_ds);
  g.SetSymbolSize(0.25, lj_rep_ds);
  //endregion
  //endregion

  //region ten Wolde-Frenkel:
  //region Potential:
  auto twf_vector = twf(x_vector);
  auto twf_ds = g.AddDataset(conv_arma_to_vec(x_vector), conv_arma_to_vec(twf_vector));
  g.SetColor(dft_core::grace_plot::Color::BLUE,twf_ds);
  //endregion

  //region Minimum:
  auto twf_min = g.AddDataset(std::vector<double>{twf.r_min()}, std::vector<double>{twf.v_min()});
  g.SetLineType(dft_core::grace_plot::LineStyle::NO_LINE, twf_min);
  g.SetSymbol(dft_core::grace_plot::Symbol::DIAMOND, twf_min);
  g.SetSymbolFill(dft_core::grace_plot::Color::BLUE, twf_min);
  //endregion

  //region Hard-sphere diameter:
  hs_x /= hs_diameter;
  hs_diameter = twf.FindHardSphereDiameter(1.0); hs_x *= hs_diameter;
  auto twf_hs = g.AddDataset(conv_arma_to_vec(hs_x), conv_arma_to_vec(y_vec));
  g.SetColor(dft_core::grace_plot::Color::BLUE, twf_hs);
  g.SetLineType(dft_core::grace_plot::LineStyle::DASHEDLINE_EN, twf_hs);
  console::Info("tWF hard-sphere diameter (kT = 1.0) = " + std::to_string(hs_diameter));
  //endregion

  //region Pertubation theory:
  auto twf_att = twf.w_attractive(x_vector);
  auto twf_att_ds = g.AddDataset(conv_arma_to_vec(x_vector), conv_arma_to_vec(twf_att));
  g.SetColor(dft_core::grace_plot::Color::BLUE, twf_att_ds);
  g.SetLineType(dft_core::grace_plot::LineStyle::DOTTEDLINE, twf_att_ds);
  g.SetSymbol(dft_core::grace_plot::Symbol::STAR, twf_att_ds);
  g.SetSymbolSize(0.25, twf_att_ds);

  auto twf_rep = twf.w_repulsive(x_vector);
  auto twf_rep_ds = g.AddDataset(conv_arma_to_vec(x_vector), conv_arma_to_vec(twf_rep));
  g.SetColor(dft_core::grace_plot::Color::BLUE, twf_rep_ds);
  g.SetLineType(dft_core::grace_plot::LineStyle::D_DOTTEDDASHEDLINE_EN, twf_rep_ds);
  g.SetSymbol(dft_core::grace_plot::Symbol::STAR, twf_rep_ds);
  g.SetSymbolSize(0.25, twf_rep_ds);
  //endregion
  //endregion

  //region Wang-Ramirez-Dobnikar-Frenkel:
  //region Potential:
  auto wrdf_vector = wrdf(x_vector);
  auto wrdf_ds = g.AddDataset(conv_arma_to_vec(x_vector), conv_arma_to_vec(wrdf_vector));
  g.SetColor(dft_core::grace_plot::Color::ORANGE, wrdf_ds);
  //endregion

  //region Minimum:
  auto wrdf_min = g.AddDataset(std::vector<double>{wrdf.r_min()}, std::vector<double>{wrdf.v_min()});
  g.SetLineType(dft_core::grace_plot::LineStyle::NO_LINE, wrdf_min);
  g.SetSymbol(dft_core::grace_plot::Symbol::SQUARE, wrdf_min);
  g.SetSymbolFill(dft_core::grace_plot::Color::ORANGE, wrdf_min);
  //endregion

  //region Hard-sphere diameter:
  hs_diameter = wrdf.FindHardSphereDiameter(1.0);
  hs_x = arma::vec(10, arma::fill::ones); hs_x *= hs_diameter;
  auto wrdf_hs = g.AddDataset(conv_arma_to_vec(hs_x), conv_arma_to_vec(y_vec));
  g.SetColor(dft_core::grace_plot::Color::ORANGE, wrdf_hs);
  g.SetLineType(dft_core::grace_plot::LineStyle::DASHEDLINE_EN, wrdf_hs);
  console::Info("WRDF hard-sphere diameter (kT = 1.0) = " + std::to_string(hs_diameter));
  //endregion

  //region Pertubation theory:
  auto wrdf_att = wrdf.w_attractive(x_vector);
  auto wrdf_att_ds = g.AddDataset(conv_arma_to_vec(x_vector), conv_arma_to_vec(wrdf_att));
  g.SetColor(dft_core::grace_plot::Color::ORANGE, wrdf_att_ds);
  g.SetLineType(dft_core::grace_plot::LineStyle::D_DOTTEDDASHEDLINE_EM, wrdf_att_ds);
  g.SetSymbol(dft_core::grace_plot::Symbol::PLUS, wrdf_att_ds);
  g.SetSymbolSize(0.25, wrdf_att_ds);

  auto wrdf_rep = wrdf.w_repulsive(x_vector);
  auto wrdf_rep_ds = g.AddDataset(conv_arma_to_vec(x_vector), conv_arma_to_vec(wrdf_rep));
  g.SetColor(dft_core::grace_plot::Color::ORANGE, wrdf_rep_ds);
  g.SetLineType(dft_core::grace_plot::LineStyle::D_DOTTEDDASHEDLINE_EN, wrdf_rep_ds);
  g.SetSymbol(dft_core::grace_plot::Symbol::STAR, wrdf_rep_ds);
  g.SetSymbolSize(0.25, wrdf_rep_ds);
  //endregion
  //endregion

  //region Legend:
  g.SetLegend("Lennard-Jones (LJ)", lj_ds);
  g.SetLegend("ten Wolde-Frenkel (tWF)", twf_ds);
  g.SetLegend("Wang-Ramirez-Dobnikar-Frenkel (WRDF)", wrdf_ds);

  g.SetLegend("LJ min", lj_min);
  g.SetLegend("LJ attractive", lj_att_ds);
  g.SetLegend("LJ repulsive", lj_rep_ds);
  g.SetLegend("LJ HS diameter", lj_hs);

  g.SetLegend("tWF min", twf_min);
  g.SetLegend("tWF HS diameter", twf_hs);
  g.SetLegend("tWF attractive", twf_att_ds);
  g.SetLegend("tWF repulsive", twf_rep_ds);

  g.SetLegend("WRDF min", wrdf_min);
  g.SetLegend("WRDF HS diameter", wrdf_hs);
  g.SetLegend("WRDF attractive", wrdf_att_ds);
  g.SetLegend("WRDF repulsive", wrdf_rep_ds);
  //endregion

  g.RedrawAndWait();

  //g.PrintToFile("potentials.png", dft_core::grace_plot::ExportFormat::PNG);
  //g.RedrawAndWait();
}