# Geometry wrapper

### Introduction

One of the main ingredients required when solving a `classicalDFT` problem is the mesh on which the partial differential equations will be integrated over time. The mesh is nothing but an abstraction of the physical space, hence it can be 1D, 2D or even 3D, depending on the type of problem we are aiming at solving. This abstraction comes already implemented under the namespace `dft_core::geometry`. Under such namespace, we can find both `abstract classes` which serve the purpose of base templates, and classes which inherit from such abstractions to create the specific mesh we want in the dimension we desire. 

The base classes are:

* **Vertex**: A class which aims at representing a vertex on a mesh of any given dimension
* **Element**: A set of vertices which represent a distinctive element of the mesh (e.g, a square in a 2D mesh)
* **Mesh**: A set of vertices and elements which will represent a physical geometry

The instantiatable classes which are inherited from these three are categorised according to the dimensionality of the physical domain where they are going be living, i.e. [2D](core/dft_lib/geometry/2D) and [3D](core/dft_lib/geometry/3D).

### Example

The best way to show the convenience offered by `geometry` is by example. The following simple piece of code shows some of the functionality mentioned above:

```c++
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
  console::Info("Square boxes 2D: default");
  auto default_box_2D = two_dimensional::SquareBox();
  std::cout << "Default square-box:" << std::endl
            << default_box_2D << std::endl;

  console::Info("Square boxes 2D: customized");
  auto custom_box_2D = two_dimensional::SquareBox(0.25, {0,0});
  std::cout << "Customised square-box:" << std::endl
            << custom_box_2D << std::endl;

  auto v_list_b = std::vector<Vertex>{ {0,0}, {0,1}, {1,1}, {1,0} };
  std::cout << "Move-semantics square-box:" << std::endl
            << "v_list[3]: " << v_list_b[3] << std::endl
            << "v_list[2]: " << v_list_b[2] << std::endl
            << "v_list[1]: " << v_list_b[1] << std::endl
            << "v_list[0]: " << v_list_b[0] << std::endl << std::endl;

  auto move_box_2D = two_dimensional::SquareBox(std::move(v_list_b));
  std::cout << move_box_2D << std::endl;
  // endregion

  // region 2D::SUQMesh
  console::Info("Mesh 2D: lattice");

  auto origin_2D = std::vector<double>{0,0};
  auto lengths_2D = std::vector<double>{1.0,1.0};
  auto lattice_2D = two_dimensional::Lattice(0.25, lengths_2D, origin_2D);

  std::cout << lattice_2D << std::endl;

  std::cout << "Vertex[-1,-1]: " << lattice_2D[{-1,-1}] << std::endl;
  std::cout << "Vertex[2,3]: " << lattice_2D[{2,3}] << std::endl;

  std::cout << "Element[10]: " << std::endl;
  std::cout << lattice_2D.elements()[10] << std::endl;

  std::cout << "Element volume: " << std::endl
            << lattice_2D.element_volume() << std::endl;

  console::Info("Mesh 2D: lattice plot");
  lattice_2D.plot();

  // endregion

  // region 3D::SquareBox
  console::Info("Square boxes 3D: default");
  auto default_box_3D = three_dimensional::SquareBox();
  std::cout << "Default square-box (3D):" << std::endl
            << default_box_3D << std::endl;

  console::Info("Mesh 3D: lattice");
  auto origin_3D = std::vector<double>{0,0,0};
  auto lengths_3D = std::vector<double>{1.0,1.0,1.0};
  auto lattice_3D = three_dimensional::Lattice(0.25, lengths_3D, origin_3D);

  std::cout << lattice_3D << std::endl;
  std::cout << "Vertex[-1,-1,-1]: " << lattice_3D[{-1,-1,-1}] << std::endl;
  std::cout << "Vertex[2,3,1]: " << lattice_3D[{2,3,1}] << std::endl;

  std::cout << "Element[10]: " << std::endl;
  std::cout << lattice_3D.elements()[10] << std::endl;

  std::cout << "Element volume: " << std::endl
            << lattice_3D.element_volume() << std::endl;
  // endregion
}
```

