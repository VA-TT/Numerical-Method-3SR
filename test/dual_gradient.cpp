#include "../library/DualDiffrentiation.h"
#include <cmath>
#include <iostream>

int main() {
  // test point
  double x0 = 1.2345;
  double y0 = 0.7654;

  // example scalar function f(x,y) = sin(x)*y^2 + x*y
  auto f = [](Dual x, Dual y) -> Dual { return sin(x) * y * y + x * y; };

  auto grad = gradient2(f, x0, y0);

  // finite-difference check (central differences)
  double h = 1e-7;
  auto ff = [&](double x, double y) { return std::sin(x) * y * y + x * y; };
  double dfdx_fd = (ff(x0 + h, y0) - ff(x0 - h, y0)) / (2.0 * h);
  double dfdy_fd = (ff(x0, y0 + h) - ff(x0, y0 - h)) / (2.0 * h);

  std::cout.setf(std::ios::fixed);
  std::cout.precision(10);
  std::cout << "AD partials: df/dx = " << grad.first
            << ", df/dy = " << grad.second << "\n";
  std::cout << "FD approx  : df/dx = " << dfdx_fd << ", df/dy = " << dfdy_fd
            << "\n";
  std::cout << "Errors     : |dx|=" << std::fabs(grad.first - dfdx_fd)
            << ", |dy|=" << std::fabs(grad.second - dfdy_fd) << "\n";

  return 0;
}
