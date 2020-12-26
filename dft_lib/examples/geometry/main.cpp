#include <cmath>
#include <armadillo>
#include "classical_dft"

int main(int argc, char **argv)
{
  using namespace dft_core::geometry;

  // region Vertex: Initializer cttor
  console::Info("Vertex | Constructor");
  auto v1 = Vertex({0, 1, 2, 3});
  auto v2 = Vertex({3, 4, 5, 6});

  std::cout << "v1: " << v1 << std::endl;
  std::cout << "v2: " << v2 << std::endl << std::endl;
  // endregion

  // region Element: Cttor with pass-by-ref std::vector<Vertex> (copy)
  auto v_list = std::vector<Vertex>{ v1, v2 };

  console::Info("Element | Copy constructor");
  auto e1 = Element(v_list);

  std::cout << "v_list[0]: " << v_list[0] << std::endl
            << "v_list[1]: " << v_list[1] << std::endl << std::endl;

  std::cout << "Element 1: " << std::endl
            << "\u2b91 raw[0] = " << e1.vertices_raw()[0] << std::endl
            << "\u2b91 map[0]  = " << e1.vertices().at(0).get() << std::endl
            << "\u2b91 raw[1] = " << e1.vertices_raw()[1] << std::endl
            << "\u2b91 map[1]  = " << e1.vertices().at(1).get() << std::endl << std::endl;

  // endregion

  // region Element: Cttor with std::move of std::vector<Vertex> (move)
  console::Info("Element | Move semantics");

  std::cout << "v_list[0]: " << v_list[0] << std::endl
            << "v_list[1]: " << v_list[1] << std::endl << std::endl;

  auto e2 = Element(std::move(v_list));

  std::cout << "Element 2: " << std::endl
            << "\u2b91 raw[0] = " << e2.vertices_raw()[0] << std::endl
            << "\u2b91 map[0]  = " << e2[0] << std::endl
            << "\u2b91 raw[1] = " << e2.vertices_raw()[1] << std::endl
            << "\u2b91 map[1]  = " << e2[1] << std::endl << std::endl;

  std::cout << e2 << std::endl;
  // endregion

  // region 2D: SquareBox
  console::Info("Square boxes:");
  auto default_box = two_dimensional::SquareBox();
  std::cout << "Default square-box:" << std::endl
            << default_box << std::endl;

  auto custom_box = two_dimensional::SquareBox(0.25, {0,0});
  std::cout << "Customised square-box:" << std::endl
            << custom_box << std::endl;

  auto v_list_b = std::vector<Vertex>{ {0,0}, {0,1}, {1,1}, {1,0} };
  std::cout << "Move-semantics square-box:" << std::endl
            << "v_list[3]: " << v_list_b[3] << std::endl
            << "v_list[2]: " << v_list_b[2] << std::endl
            << "v_list[1]: " << v_list_b[1] << std::endl
            << "v_list[0]: " << v_list_b[0] << std::endl << std::endl;

  auto move_box = two_dimensional::SquareBox(std::move(v_list_b));
  std::cout << move_box << std::endl;
  // endregion

  auto origin = std::vector<double>{0,0};
  auto lengths = std::vector<double>{1.0,1.0};
  auto suq_mesh = two_dimensional::SUQMesh(0.25, lengths, origin);

  std::cout << suq_mesh << std::endl;
  std::cout << "Vertex[-1,-1]: " << suq_mesh[{-1,-1}] << std::endl;
  std::cout << suq_mesh.elements()[10] << std::endl;

  suq_mesh.plot();
}
