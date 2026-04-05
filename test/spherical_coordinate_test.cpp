#include "../library/coordinates.h"

#include <cassert>
#include <cmath>
#include <iostream>

int main() {
  const double r0 = 2.0;
  const double t0 = 60.0; // degree
  const double p0 = 40.0; // degree

  SphereCoor<double> s{r0, t0, p0};

  constexpr double tol = 1e-6;

  // Jacobian: J = r^2 sin(theta)
  const double jac_expected = r0 * r0 * std::sin(t0 * constants::pi / 180.0);
  assert(std::abs(s.jacobian() - jac_expected) < tol);

  // Scalar field f(r,t,p) = r^2
  auto f = [](Dual r, Dual /*t*/, Dual /*p*/) -> Dual { return r * r; };

  const auto g = s.gradient(f);
  assert(std::abs(g[0] - 2.0 * r0) < tol);
  assert(std::abs(g[1]) < tol);
  assert(std::abs(g[2]) < tol);

  const double lap_f = s.laplacian(f, 1e-5, 1e-5, 1e-5);
  assert(std::abs(lap_f - 6.0) < 1e-4);

  // Vector field u = (r^2, 0, 0)
  auto u = [](Dual r, Dual /*t*/, Dual /*p*/) -> StaticVector<Dual, 3> {
    return {r * r, Dual{0.0}, Dual{0.0}};
  };

  const double div_u = s.div(u);
  assert(std::abs(div_u - 4.0 * r0) < tol);

  const auto curl_u = s.curl(u);
  assert(std::abs(curl_u[0]) < tol);
  assert(std::abs(curl_u[1]) < tol);
  assert(std::abs(curl_u[2]) < tol);

  const auto lap_u = s.laplacian(u, 1e-5, 1e-5, 1e-5);
  assert(std::abs(lap_u[0] - 6.0) < 1e-4);
  assert(std::abs(lap_u[1]) < 1e-4);
  assert(std::abs(lap_u[2]) < 1e-4);

  std::cout << "spherical_coordinate_test passed\n";
  return 0;
}
