#include "../library/newtonRaphson.h"
#include <iostream>

// Define your custom function here (not hardcoded in header!)
// Example: 8x³ - 4x² - 3x + 1
auto myFunc(auto x) {
  static const std::vector<double> coef = {1.0, -3.0, -4.0,
                                           8.0}; // a0, a1, a2, a3
  return polyFunc(coef, x);
}

// Derivative using automatic differentiation
double myDerivative(double x) {
  auto func = [](Dual d) { return myFunc(d); };
  return automaticDiff(func, x);
}

int main() {
  std::cout << "=== Newton-Raphson Root Finding ===\n\n";

  // Method 1: Using generic Newton-Raphson with custom function
  std::cout << "Method 1: Custom function with newtonRaphsonGeneric\n";
  auto f = [](double x) { return myFunc(x); };
  auto df = [](double x) { return myDerivative(x); };
  double root1 = newtonRaphsonGeneric(f, df, 0.5, 1e-10);
  std::cout << "  Root at x = " << root1 << ", f(x) = " << myFunc(root1)
            << "\n\n";

  // Method 2: Using polynomial coefficients directly
  std::cout << "Method 2: Using polynomial with coefficients\n";
  std::vector<double> coef = {1.0, -3.0, -4.0, 8.0}; // 1 - 3x - 4x² + 8x³
  double root2 = newtonRaphsonPoly(coef, 0.5, 1e-10);
  std::cout << "  Root at x = " << root2 << ", f(x) = " << polyFunc(coef, root2)
            << "\n\n";

  // Method 3: Find ALL roots using polynomial deflation
  std::cout << "Method 3: Finding all roots with deflation\n";
  std::vector<double> guesses = {-10.0, -1.0, 0.0, 0.5, 1.0, 10.0};
  auto all_roots = findAllRoots(coef, guesses, 1e-10);

  std::cout << "  Found " << all_roots.size() << " roots:\n";
  for (const auto &r : all_roots) {
    std::cout << "    x = " << r << ", f(x) = " << polyFunc(coef, r) << '\n';
  }

  return 0;
}