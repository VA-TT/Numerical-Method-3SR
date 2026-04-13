#include "../library/coordinates.h"
#include "../library/gaussQuadrature.h"

#include <cassert>
#include <cmath>
#include <iostream>

int main() {
  const double k1 = 2.5;
  const double k2 = -1.75;

  auto u_rt = [k1, k2](Dual r, Dual /*t*/) -> StaticVector<Dual, 2> {
    return {k1 / r, k2 / r};
  };

  // Evaluate at several points r != 0.
  struct Sample {
    double r;
    double thetaDeg;
  };

  const Sample samples[] = {
      {0.5, 10.0},
      {1.0, 45.0},
      {2.0, 120.0},
      {3.5, 270.0},
  };

  constexpr double tol = 1e-10;

  for (const auto &s : samples) {
    PolarCoord<double> p{s.r, s.thetaDeg};

    const double div_u = p.div(u_rt);
    const double curl_u = p.curl(u_rt);

    assert(std::abs(div_u) < tol && "Expected div(u) = 0 for r != 0");
    assert(std::abs(curl_u) < tol && "Expected curl(u) = 0 for r != 0");

    // f(r,theta) = r^2 => Laplacian in 2D polar is 4.
    auto f_rt = [](Dual r, Dual /*t*/) -> Dual { return r * r; };
    const double lap_f = p.laplacian(f_rt, 1e-5, 1e-5);
    assert(std::abs(lap_f - 4.0) < 1e-5 && "Expected laplacian(r^2) = 4");
  }

  // Jacobian test via Gauss quadrature:
  // A = integral_0^{2pi} integral_0^R (J dr dtheta), with J = r.
  const double R = 2.0;
  const auto jacobian_integrand = [](double r, double theta_rad) -> double {
    const double theta_deg = theta_rad * 180.0 / constants::pi;
    PolarCoord<double> p{r, theta_deg};
    return p.jacobian();
  };
  const double area_num = integrationGauss2D(0.0, R, 0.0, 2.0 * constants::pi,
                                             jacobian_integrand, 4);

  const double area_exact = constants::pi * R * R;
  const double area_err = std::abs(area_num - area_exact);
  assert(area_err < 1e-3 &&
         "Expected area from Jacobian integral to match pi R^2");

  std::cout << "polar_coordinate_test passed: div(u)=0, curl(u)=0, "
               "laplacian(r^2)=4, Jacobian area test passed\n";
  return 0;
}
