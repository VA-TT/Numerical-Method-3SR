#include "../library/interpolate.h"
#include <fstream>
#include <iomanip>
#include <iostream>

int main() {
  constexpr Index nPts = 8;
  StaticVector<double, nPts> data_x{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
  StaticVector<double, nPts> data_y{0.0, 1.0, 1.0, 0.5, 0.2, -0.3, 0.1, 0.8};

  std::ofstream out("interpolation_three_methods.dat");
  if (!out) {
    std::cerr << "Cannot open output file: interpolation_three_methods.dat\n";
    return 1;
  }

  out << "# x polynomial linear cubic_spline\n";
  out << std::fixed << std::setprecision(8);

  const int nSamples = 200;
  const double xMin = data_x[0];
  const double xMax = data_x[nPts - 1];
  for (int i = 0; i <= nSamples; ++i) {
    double x = xMin + (xMax - xMin) * static_cast<double>(i) /
                          static_cast<double>(nSamples);
    double y_poly = interpolatePolynomial(data_x, data_y, x);
    double y_lin = interpolateLinear(data_x, data_y, x);
    double y_spline = interpolateCubicSpline(data_x, data_y, x);
    out << x << " " << y_poly << " " << y_lin << " " << y_spline << "\n";
  }

  std::cout << "Wrote interpolation_three_methods.dat\n";
  return 0;
}
