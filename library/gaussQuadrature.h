#ifndef GAUSS_LEGENDRE_QUADRATURE_H
#define GAUSS_LEGENDRE_QUADRATURE_H

#include "squareRoot.h"
#include <array>
#include <cassert>
#include <functional>
#include <iostream>

// Note: In this method, we should use (n) Gauss points to approximate a
// polynomials with degree (2n - 1),

namespace gaussQuadrature {
constexpr std::array<int, 4> pointsNumber = {1, 2, 3, 4};

// Gauss Points
constexpr std::array<double, 1> xiOne = {0};
constexpr std::array<double, 2> xiTwo = {-1.0 / constexpr_sqrt(3.0),
                                         1.0 / constexpr_sqrt(3.0)};
constexpr std::array<double, 3> xiThree = {-constexpr_sqrt(3.0 / 5.0), 0.0,
                                           constexpr_sqrt(3.0 / 5.0)};
constexpr std::array<double, 4> xiFour = {
    -constexpr_sqrt(3.0 / 7.0 + 2.0 / 7.0 * constexpr_sqrt(6.0 / 5.0)),
    -constexpr_sqrt(3.0 / 7.0 - 2.0 / 7.0 * constexpr_sqrt(6.0 / 5.0)),
    constexpr_sqrt(3.0 / 7.0 - 2.0 / 7.0 * constexpr_sqrt(6.0 / 5.0)),
    constexpr_sqrt(3.0 / 7.0 + 2.0 / 7.0 * constexpr_sqrt(6.0 / 5.0))};

// Weights coefficients
constexpr std::array<double, 1> wOne = {2.0};
constexpr std::array<double, 2> wTwo = {1.0, 1.0};
constexpr std::array<double, 3> wThree = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
constexpr std::array<double, 4> wFour = {
    (18.0 - constexpr_sqrt(30.0)) / 36.0, (18.0 + constexpr_sqrt(30.0)) / 36.0,
    (18.0 + constexpr_sqrt(30.0)) / 36.0, (18.0 - constexpr_sqrt(30.0)) / 36.0};

// Function to set up the weight coefficient and points of integration based on
// number of integration points.
inline void valueXiW(int n, const double *&xi_ptr, const double *&w_ptr) {
  assert(n >= 1 && n <= 4 && "Number of Gauss points must be in [1,4]");

  // Uses pointers and member function .data() to the xi/w arrays to avoid
  // copying; selects arrays by n.
  switch (n) {
  case 1:
    xi_ptr = xiOne.data();
    w_ptr = wOne.data();
    break;
  case 2:
    xi_ptr = xiTwo.data();
    w_ptr = wTwo.data();
    break;
  case 3:
    xi_ptr = xiThree.data();
    w_ptr = wThree.data();
    break;
  case 4:
    xi_ptr = xiFour.data();
    w_ptr = wFour.data();
    break;
  }
}

}; // namespace gaussQuadrature

// Generic 1-D Gauss integration on [a,b] using n points (1 <= n <= 4).
// By default n = 2 for light calculation.
template <typename T, typename Func>
T integrationGauss1D(T x1, T x2, Func f, int n = 2) {
  const double *xi_ptr = nullptr;
  const double *w_ptr = nullptr;
  gaussQuadrature::valueXiW(n, xi_ptr, w_ptr);
  T I{0.0};
  T c1{(x2 - x1) / T{2}}; // mapping from [-1,1] to [a,b]
  T c2{(x2 + x1) / T{2}};
  for (int i = 0; i < n; ++i) {
    I += c1 * static_cast<T>(w_ptr[i]) * f(c1 * static_cast<T>(xi_ptr[i]) + c2);
  }
  return I;
}

// Implement 2D: Intergration in both x and y direction for a rectangular shape
// Assuming that number of Gauss points taken in both direction is the same

template <typename T, typename Func>
T integrationGauss2D(T x1, T x2, T y1, T y2, Func f, int n = 2) {
  const double *xi_ptr = nullptr;
  const double *w_ptr = nullptr;
  gaussQuadrature::valueXiW(n, xi_ptr, w_ptr);
  // use the same xi/w in both directions (tensor product)
  // Not actually necessary but choose to do like this for future implement
  // For example : triangular element:  m_x != m_y
  const double *xi_x = xi_ptr;
  const double *xi_y = xi_ptr;
  const double *w_x = w_ptr;
  const double *w_y = w_ptr;
  int m_x = n;
  int m_y = n;

  T I{0.0};

  T c1{(x2 - x1) / T{2}};
  T c2{(x1 + x2) / T{2}};
  T d1{(y2 - y1) / T{2}};
  T d2{(y1 + y2) / T{2}};
  for (int i = 0; i < m_x; ++i) {
    for (int j = 0; j < m_y; ++j) {
      I += c1 * d1 * static_cast<T>(w_x[i]) * static_cast<T>(w_y[j]) *
           f(c1 * static_cast<T>(xi_x[i]) + c2,
             d1 * static_cast<T>(xi_y[j]) + d2);
    }
  }
  return I;
}
#endif

// // Sample function
// double test(double x) { return x * x * x + 1; }
// double testfunctionUnit(double x, double y) { return (x * x + y * y); }
// int main() {
//   double expX{integrationGauss1D(-1, 1, test, 4)};
//   double expXY{integrationGauss2D(-2, 3, -4, 5, testfunctionUnit, 4)};
//   double expXY_default{integrationGauss2D(-2, 3, -4, 5, testfunctionUnit)};

//   // Print with 15 decimal precisions
//   std::cout << std::fixed << std::setprecision(15);
//   std::cout << "Integration [-1,1] of function x^3: I = " << expX <<
//   std::endl; std::cout << "Integration [-2,3]x[-4,5] of function (x^2+y^2): I
//   = " << expXY
//             << ' ' << expXY_default << std::endl;

//   return 0;
// }