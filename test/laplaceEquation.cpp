#include "../library/DualDifferentiation.h"
#include "../library/Matrix.h"
#include "../library/Mesh2.h"
#include "../library/Vector.h"
#include "../library/interpolate.h"
#include <iostream>

int main() {
  constexpr Index nx{101}, ny{101};
  double l{1.0}, h{1.0};
  Mesh2D<double> domain{l, h, nx, ny};
  const auto &nodes = domain.getAllNodes();
  
}
