# Intermolecular potentials

### Introduction

One of the few and key elements required when doing research on statistical physics is to have a good estimate (or at least a justifiable model) for the intermolecular potential energy between a pair of constituent particles, e.g. atoms or molecules. The choice of the *intermolecular potential* model describing the intermolecular forces is not something we should come up with hastily, since this **fundamental quantity** governs the thermodynamics of the system at hand. I.e., the same simulation with the same boundary and initial conditions will very likely result in quite different outcomes if we choose to model our system with a [Lennard-Jones](https://en.wikipedia.org/wiki/Lennard-Jones_potential) (LJ) or a [ten Wolde-Frenkel](https://science.sciencemag.org/content/277/5334/1975) (tWF) potential. For this reason, a myriad of potentials have been developed ever since the introduction of the widely-known **LJ(6-12) potential** in 1924 by [John Lennard-Jones](https://en.wikipedia.org/wiki/John_Lennard-Jones). However, most of them share common elements which can be encapsulated in a `base class` to be reutilised for the implementation of every particular one of them. 

The `classicalDFT` library comes with an implementation of such an abstract class which can be found in the namespace `dft_core::physics::potentials::intermolecular`, class `Potential`. Besides this abstract class, `classicalDFT` offers the implementation of the two intermolecular potentials mentioned above, namely LJ(6-12) and the tWF potential, `LennardJones` and `tenWoldeFrenkel` classes under the same namespace. However, the implementations offered here not only come with the basic functionality, but also are enriched with some convenient methods and properties which are required for thermodynamic perturbation theory (TPT) analysis. Particularly, any potential which derives from the abstract class `Potential` will automatically know how to split the potential contribution into two parts (as suggested by [Weeks-Chandler-Andersen TPT](http://www.sklogwiki.org/SklogWiki/index.php/Weeks-Chandler-Andersen_perturbation_theory)): a) hard-sphere (purely repulsive) contribution; and b) the attractive part. E.g., for the particular case of the LJ(6-12) potential:
$$
V_{\text{LJ}}(\mathbf{r}; \varepsilon, \sigma) = 4\varepsilon\left[\left(\frac{\sigma}{|\mathbf{r}|}\right)^{12}-\left(\frac{\sigma}{|\mathbf{r}|}\right)^{6}\right]\doteq\Phi_{\text{LJ}}(\mathbf{r}; \varepsilon, \sigma)
$$
The WCA-PT assumes the following decomposition of the inter-particle potential:
$$
V_{\text{LJ}}(\mathbf{r};\varepsilon,\sigma)=\Phi_{\text{repulsive}}(\mathbf{r}; \varepsilon, \sigma) + \Phi_{\text{attractive}}(\mathbf{r}; \varepsilon, \sigma)
$$
where:
$$
\Phi_{\text{repulsive}}(\mathbf{r};\varepsilon, \sigma) = \Theta(r_*-|\mathbf{r}|)\times(\Phi_{\text{LJ}}
(\mathbf{r}; \varepsilon, \sigma) - \Phi_{\text{LJ}}(r_*))
$$
where $r_*$ is the distance at which the potential reaches its minimum: $\nabla\Phi|_{r_*}=0$, and
$$
\Phi_{\text{attractive}}(\mathbf{r};\varepsilon, \sigma) = \Theta(r_*-|\mathbf{r}|)\times\Phi_{\text{LJ}}(r_*)
+ \Theta(|\mathbf{r}|-r_*)\times\Phi_{\text{LJ}}(\mathbf{r};\varepsilon,\sigma)
$$
The `Potential` class comes with a generalisation of this decomposition which allows for a different position where to split the interaction, and the underlying potential $\Phi$ was renamed as `v_potential`. Using the built-in method `Potential::SetBHPerturbation()` the split occurs at the position where the potential is zero, i.e.  $V(r_*) = 0$.

The potential-like object also comes with the functionality to compute the hard-sphere, by using the approximation:
$$
d_{\text{HS}} (T) = \int_{r_{\text{HC}}}^{r_*}\left(1-e^{-V(r)/k_BT}\right) dr
$$
Where, $T$ is the absolute temperature of the system (which typically is given through $\beta=(k_BT)^{-1}$); $r_{\text{HC}}$ is the hard-core diameter related to the interaction potential at hand.

### Examples

The best way of showing the convenience offered by the `physics::potentials::intermolecular::Potential` is by example. Thus, we are going to proceed by inserting the code in [`main.cpp`](main.cpp):

```c++
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
```

After compilation and running we will get the following results:

<img src="figures/potentials.png" alt="potentials" style="zoom:50%;" />

The figure's quality might seem poor because it has been taken from a screenshot. For a better resolution and quality, you can always use one of the export methods provided by `grace` (as is shown at the bottom of the above code block).