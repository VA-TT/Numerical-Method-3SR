#include "../library/coordinates.h"

#include <cassert>
#include <cmath>
#include <iostream>

int main() {
  // Test at a regular point (r != 0).
  const double r0 = 2.0;
  const double p0_deg = 30.0;
  const double z0 = 3.0;
  CylinCoord<double> c{r0, p0_deg, z0};

  constexpr double tol = 1e-8;

  // Scalar field: f(r,p,z) = r^2 + z^2
  auto f = [](Dual r, Dual /*p*/, Dual z) -> Dual { return r * r + z * z; };

  const auto g = c.gradient(f);
  assert(std::abs(g[0] - 2.0 * r0) < tol);
  assert(std::abs(g[1] - 0.0) < tol);
  assert(std::abs(g[2] - 2.0 * z0) < tol);

  // Vector field: u = (r^2, r*z, z^2)
  auto u = [](Dual r, Dual /*p*/, Dual z) -> StaticVector<Dual, 3> {
    return {r * r, r * z, z * z};
  };

  const double div_u = c.div(u);
  const double div_expected = 3.0 * r0 + 2.0 * z0;
  assert(std::abs(div_u - div_expected) < tol);

  const auto curl_u = c.curl(u);
  assert(std::abs(curl_u[0] - (-r0)) < tol);
  assert(std::abs(curl_u[1] - 0.0) < tol);
  assert(std::abs(curl_u[2] - (2.0 * z0)) < tol);

  // Scalar Laplacian: f = r^2 + z^2 => 6
  const double lap_f = c.laplacian(f, 1e-5, 1e-5, 1e-5);
  assert(std::abs(lap_f - 6.0) < 1e-4);

  // Vector Laplacian component-wise test:
  // v = (r^2, z^2, r*z) => (4, 2, z/r)
  auto v = [](Dual r, Dual /*p*/, Dual z) -> StaticVector<Dual, 3> {
    return {r * r, z * z, r * z};
  };
  const auto lap_v = c.laplacian(v, 1e-5, 1e-5, 1e-5);
  assert(std::abs(lap_v[0] - 4.0) < 1e-4);
  assert(std::abs(lap_v[1] - 2.0) < 1e-4);
  assert(std::abs(lap_v[2] - (z0 / r0)) < 1e-3);

  // Jacobian check
  assert(std::abs(c.jacobian() - r0) < tol);

  std::cout << "cylindrical_coordinate_test passed\n";
  return 0;
}
