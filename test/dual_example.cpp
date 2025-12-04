#include "../library/DualDiffrentiation.h"
#include <cmath>
#include <functional>
#include <iostream>

int main() {
  // point where we evaluate derivative
  double x0 = 1.234567;

  // Example function: f(x) = sin(x) * x^2
  auto f_dual = [](Dual x) -> Dual {
    return sin(x) * x * x; // uses Dual overloads
  };

  // Use automaticDiff template (callable taking Dual) to get derivative
  double df = automaticDiff<decltype(f_dual)>(f_dual, x0);

  // To get the function value and derivative together, evaluate on a Dual with
  // der=1
  Dual xd{x0, 1.0};
  Dual y = f_dual(xd);

  // Numeric finite-difference check (central difference)
  double h = 1e-6;
  auto f_double = [](double x) -> double { return std::sin(x) * x * x; };
  double df_fd = (f_double(x0 + h) - f_double(x0 - h)) / (2.0 * h);

  std::cout.setf(std::ios::fixed);
  std::cout.precision(10);
  std::cout << "Point x0 = " << x0 << "\n";
  std::cout << "AD derivative df/dx = " << df << "\n";
  std::cout << "Dual evaluation: value = " << y.getVal()
            << ", der = " << y.getDer() << "\n";
  std::cout << "Finite-diff approx = " << df_fd << "\n";
  std::cout << "Absolute error (AD vs FD) = " << std::fabs(df - df_fd) << "\n";

  return 0;
}
