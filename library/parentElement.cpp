#ifndef PARENT_ELEMENT_SHAPE_FUNCTION
#define PARENT_ELEMENT_SHAPE_FUNCTION
#include "gaussQuadrature.h"
#include <functional>
// Range of 1D parent element Sr: [-1,1] -> mapping to [x1,x2] physical element
// First order polynomial shape function
auto N1_r = [](auto psy) { return (1 + psy) / 2; };
auto N2_r = [](auto psy) { return (1 - psy) / 2; };

// Discretizing x = phi(psy) = x1 N1(psy) + x2 N2(psy) = J*psy + M = psy *
// (x2-x1)/2 + (x2+x1)/2

// Jacobian: J = dx/d(psy) = (x2-x1)/2
template <typename T> T Jacobian(T x1, T x2) { return (x2 - x1) / 2; }
template <typename T> T midPoint(T x1, T x2) { return (x2 + x1) / 2; }

// Mapping deriviation: dN/dx = dN/d(psy) * d(psy)/dx = dN/d(psy) * 1/J
template <typename T> T dN1(T x1, T x2) { return 0.5 / Jacobian(x1, x2); }
template <typename T> T dN2(T x1, T x2) { return -0.5 / Jacobian(x1, x2); }

// Function Phi(psy) used to map
template <typename T> T physicCoor(T psy, T x1, T x2) {
  return N1_r(psy) * x1 + N2_r(psy) * x2;
}
template <typename T> T parentCoor(T x, T x1, T x2) {
  return (x - midPoint(x1, x2)) / Jacobian(x1, x2);
}

// Interpolation: int^(x2)_(x1) f(x)dx = int^1_(-1) f(phi(psy)) J d(psy)
double integrationGauss1D_ref(double x1, double x2,
                              std::function<double(double)> f_psy, int n = 2) {
  using namespace gaussQuadrature;

  const double *xi_ptr = nullptr;
  const double *w_ptr = nullptr;
  valueXiW(n, xi_ptr, w_ptr);
  double I = 0.0;
  double J = Jacobian(x1, x2);
  for (int i = 0; i < n; ++i) {
    I += J * w_ptr[i] * f_psy(xi_ptr[i]);
  }
  return I;
}

#endif

#include <iomanip>
#include <iostream>
int main() {
  double x = 4.0, x1 = 2.0, x2 = 5.0, psy = 0.0;
  std::cout << "=== Mapping Tests ===" << '\n';
  std::cout << "Ref coordinate of x=4: psy = " << parentCoor(x, x1, x2) << '\n';
  std::cout << "Physical coordinate of psy=0: x = " << physicCoor(psy, x1, x2)
            << '\n';
  std::cout << "Jacobian: " << Jacobian(x1, x2) << '\n';
  std::cout << "Derivatives of N1: dN1 = " << dN1(x1, x2) << '\n';
  std::cout << "Derivatives of N2: dN2 = " << dN2(x1, x2) << '\n';

  std::cout << "\n=== Integration Tests ===" << '\n';
  // Test 1: Integrate N1 over element [x1,x2]
  // Analytical: int N1 dpsy = int (1+psy)/2 dpsy from -1 to 1 = [psy/2 +
  // psy^2/4] = 1 With Jacobian J: Result = J * 1 = 1.5
  auto f1 = [](double psy) { return N1_r(psy); };
  double I1 = integrationGauss1D_ref(x1, x2, f1, 2);
  std::cout << std::fixed << std::setprecision(10);
  std::cout << "Integral of N1: " << I1 << " (expected: 1.5)" << '\n';

  // Test 2: Integrate N1*N2
  // Analytical: N1*N2 = (1+psy)/2 * (1-psy)/2 = (1-psy^2)/4
  // int (1-psy^2)/4 dpsy from -1 to 1 = 1/4 * [psy - psy^3/3]|_{-1}^{1}
  // = 1/4 * [(1 - 1/3) - (-1 + 1/3)] = 1/4 * [2/3 + 2/3] = 1/4 * 4/3 = 1/3
  // With Jacobian J: Result = J * 1/3 = 1.5/3 = 0.5
  auto f2 = [](double psy) { return N1_r(psy) * N2_r(psy); };
  double I2 = integrationGauss1D_ref(x1, x2, f2, 2);
  std::cout << "Integral of N1*N2: " << I2 << " (expected: 0.5)" << '\n';
}