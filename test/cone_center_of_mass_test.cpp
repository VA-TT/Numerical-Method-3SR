#include "../library/coordinates.h"
#include "../library/gaussQuadrature.h"
#include "../library/physicConstants.h"

#include <cassert>
#include <cmath>
#include <iostream>

int main() {
  // Uniform right circular cone:
  // apex at z = 0, base at z = b, base radius = a.
  // r_max(z) = (a/b) z.

  const double a = 2.0;
  const double b = 3.0;

  // Use cylindrical Jacobian J=r from cylinCoor directly.
  auto dV_without_phi = [a, b](double z) -> double {
    const double rMax = (a / b) * z;
    auto inner = [z](double r) -> double {
      cylinCoor<double> c{r, 0.0, z};
      return c.jacobian();
    };
    return integrationGauss1D(0.0, rMax, inner, 4);
  };

  const double volume_num =
      2.0 * constants::pi * integrationGauss1D(0.0, b, dV_without_phi, 4);
  const double volume_exact = (1.0 / 3.0) * constants::pi * a * a * b;

  auto zMoment_without_phi = [a, b](double z) -> double {
    const double rMax = (a / b) * z;
    auto inner = [z](double r) -> double {
      cylinCoor<double> c{r, 0.0, z};
      return z * c.jacobian();
    };
    return integrationGauss1D(0.0, rMax, inner, 4);
  };

  const double mz_num =
      2.0 * constants::pi * integrationGauss1D(0.0, b, zMoment_without_phi, 4);

  const double z_cm_num = mz_num / volume_num;
  const double z_cm_exact = 3.0 * b / 4.0;

  constexpr double tol = 1e-10;
  assert(std::abs(volume_num - volume_exact) < tol && "Cone volume mismatch");
  assert(std::abs(z_cm_num - z_cm_exact) < tol && "Cone z_cm mismatch");

  std::cout << "cone_center_of_mass_test passed: V=" << volume_num
            << ", z_cm=" << z_cm_num << '\n';
  return 0;
}
