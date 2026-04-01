#ifndef LAGRANGE_POLYNOMIAL_INTERPOLATE_H
#define LAGRANGE_POLYNOMIAL_INTERPOLATE_H

#include "Vector.h"
#include <array>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

// Note : To interpolate a polynomal of degree n, we need n+1 points!

// Lagrange basis at x : l_i(x)
// N_i(x) = l_i(x) = Product((x - x_i) / (x_i - x_j))
template <typename T>
T basisLagrange(Index i, const DynamicVector<double> &data_x, T x) {
  T l_i{1.0};
  for (Index j{0}; j < data_x.size(); ++j) {
    if (j != i)
      l_i = l_i * ((x - data_x[j]) / (data_x[i] - data_x[j]));
  }
  return l_i;
}

// y(x) = l_i(x) y_i
// Make interpolation templated in evaluation type T so it works with Dual
// numbers produced by automatic differentiation.
template <typename T>
T interpolatePolynomial(const DynamicVector<double> &data_x,
                        const DynamicVector<double> &data_y, T x) {
  assert(data_x.size() == data_y.size() &&
         "Size of data x and y must matched!");
  T y{0.0};
  for (Index i{0}; i < data_y.size(); ++i) {
    y += data_y[i] * basisLagrange<T>(i, data_x, x);
  }
  return y;
}

// Chebyshev nodes: Coordinate of node k from n nodes in range [a,b]
double chebyShevNode(int k, int n, double a, double b) {
  double x_k = std::cos(std::numbers::pi * (2.0 * k - 1.0) / (2.0 * n));
  return 0.5 * (a + b) + 0.5 * (b - a) * x_k;
}

#endif

// int main() {

//   DynamicVector<double> data_x{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

//   DynamicVector<double> data_y{0.0, 0.8415, 0.9093, 0.1411, -0.7568, -0.9589,
//   -0.2794};

//   std::cout << std::fixed << std::setprecision(6);

//   for (double x = 0.0; x <= 6.0; x += 0.1) {
//     double fx = interpolatePolynomial(data_x, data_y, x);
//     std::cout << "P(" << x << ") = " << fx << '\n';
//   }

//   return 0;
// }
