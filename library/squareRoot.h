#ifndef MY_CONSTEXPR_SQUARE_ROOT
#define MY_CONSTEXPR_SQUARE_ROOT

#include <cmath>
#include <iomanip>
#include <iostream>

// constexpr double square_root(double a, double x0, double epsilon = 1e-10) {
//   double xn = 0.5 * (x0 + a / x0); // x1 >= sqrt(a) due to AM-GM
//   while ((xn * xn - a) > epsilon) {
//     xn = 0.5 * (xn + a / xn);
//   }
//   return xn;
// }

// constexpr absolute
constexpr double constexpr_fabs(double x) { return x < 0.0 ? -x : x; }

// constexpr sqrt using Newton with fixed cap iterations and tolerance
constexpr double constexpr_sqrt(double a, double x0 = 1.0, double eps = 1e-12,
                                int maxIter = 100) {
  if (a <= 0.0)
    return a == 0.0 ? 0.0 : (0.0 / 0.0); // NaN for negative "a"
  // Chose the beneficial first iterative value r
  double r = x0 > 0.0 ? x0 : (a >= 1.0 ? a : 1.0);
  // AM-GM
  for (int i = 0; i < maxIter; ++i) {
    double x_n = 0.5 * (r + a / r);
    if (constexpr_fabs(x_n - r) <= eps) {
      return x_n;
    }
    r = x_n;
  }
  // if not converged within maxIter, return current approx
  return r;
}

#endif

// int main() {
//   // Error eps = 1e-12
//   double x = square_root(2, 1.2, 1e-12);

//   // Print with precisions
//   std::cout << std::fixed << std::setprecision(15);
//   std::cout << "Result: " << x << std::endl;

//   // Compare with standatd sqrt
//   std::cout << "std::sqrt(2): " << std::sqrt(2) << std::endl;

//   // Errors
//   double error = std::abs(x - std::sqrt(2));
//   std::cout << std::scientific << std::setprecision(3);
//   std::cout << "Error: " << error << std::endl;

//   return 0;
// }