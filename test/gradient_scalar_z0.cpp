#include "../library/DualDiffrentiation.h"

#include <cmath>
#include <fstream>
#include <iostream>

int main() {
  // Scalar field: f(x,y,z) = 1.5*x*y + 2*z
  // We evaluate on plane z = 0 and compute gradient in (x,y) using Dual AD.

  const int nx = 13;
  const int ny = 13;
  const double xmin = 0.0, xmax = 1.0;
  const double ymin = 0.0, ymax = 1.0;
  const double z = 0.0;

  auto f_xy = [z](Dual x, Dual y) -> Dual { return 1.5 * x * y + 2.0 * z; };

  std::ofstream out("./gradient_scalar_z0.dat");
  if (!out) {
    std::cerr << "Cannot open output file report/gradient_scalar_z0.dat\n";
    return 1;
  }

  out << "# columns: x y gx gy grad_xy_norm f\n";

  // Columns: x y gx gy |grad_xy| f(x,y,0)
  for (int i = 0; i < nx; ++i) {
    const double x = xmin + (xmax - xmin) * static_cast<double>(i) / (nx - 1);
    for (int j = 0; j < ny; ++j) {
      const double y = ymin + (ymax - ymin) * static_cast<double>(j) / (ny - 1);

      const auto [gx, gy] = grad2D(f_xy, x, y);
      const double gxy = std::sqrt(gx * gx + gy * gy);
      const double f = f_xy(Dual{x, 0.0}, Dual{y, 0.0}).getVal();

      out << x << ' ' << y << ' ' << gx << ' ' << gy << ' ' << gxy << ' ' << f
          << '\n';
    }
    out << '\n';
  }

  std::cout << "Wrote ./gradient_scalar_z0.dat\n";
  return 0;
}
