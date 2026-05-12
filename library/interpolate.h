#ifndef LAGRANGE_POLYNOMIAL_INTERPOLATE_H
#define LAGRANGE_POLYNOMIAL_INTERPOLATE_H

#include "Vector.h"
#include <array>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <type_traits>

// Note : To interpolate a polynomal of degree n, we need n+1 points!

namespace interpolation_detail {

template <typename T>
void assertInsideDomain(const DynamicVector<double> &data_x, const T &x,
                        const char *funcName) {
  const double xMin = data_x[0];
  const double xMax = data_x[data_x.size() - 1];
  assert(x >= xMin && x <= xMax);
  if (x < xMin || x > xMax) {
    throw std::out_of_range(std::string(funcName) +
                            ": x is out of data_x domain");
  }
}

template <Index n, typename T>
void assertInsideDomain(const StaticVector<double, n> &data_x, const T &x,
                        const char *funcName) {
  const double xMin = data_x[0];
  const double xMax = data_x[n - 1];
  assert(x >= xMin && x <= xMax);
  if (x < xMin || x > xMax) {
    throw std::out_of_range(std::string(funcName) +
                            ": x is out of data_x domain");
  }
}

} // namespace interpolation_detail

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

template <Index n, typename T>
T basisLagrange(Index i, const StaticVector<double, n> &data_x, T x) {
  T l_i{1.0};
  for (Index j{0}; j < data_x.size(); ++j) {
    if (j != i)
      l_i = l_i * ((x - data_x[j]) / (data_x[i] - data_x[j]));
  }
  return l_i;
}

// Polynomial interpolation: y(x) = l_i(x) y_i
// Make interpolation templated in evaluation type T so it works with Dual
// numbers produced by automatic differentiation.
template <typename T>
T interpolatePolynomial(const DynamicVector<double> &data_x,
                        const DynamicVector<double> &data_y, T x) {
  assert(data_x.size() == data_y.size() &&
         "Size of data x and y must matched!");
  assert(data_x.size() >= 2 && "Must have at least 2 points!");
  interpolation_detail::assertInsideDomain(data_x, x, "interpolatePolynomial");
  T y{0.0};
  for (Index i{0}; i < data_y.size(); ++i) {
    y += data_y[i] * basisLagrange<T>(i, data_x, x);
  }
  return y;
}

template <Index n, typename T>
T interpolatePolynomial(const StaticVector<double, n> &data_x,
                        const StaticVector<double, n> &data_y, T x) {
  interpolation_detail::assertInsideDomain(data_x, x, "interpolatePolynomial");
  T y{0.0};
  for (Index i{0}; i < data_y.size(); ++i) {
    y += data_y[i] * basisLagrange(i, data_x, x);
  }
  return y;
}

// Chebyshev nodes: Coordinate of node k from n nodes in range [a,b]
double chebyShevNode(int k, int n, double a, double b) {
  double x_k = std::cos(constants::pi * (2.0 * k - 1.0) / (2.0 * n));
  return 0.5 * (a + b) + 0.5 * (b - a) * x_k;
}

// Linear piecewise interpolation:
template <typename T>
T interpolateLinear(const DynamicVector<double> &data_x,
                    const DynamicVector<double> &data_y, T x) {
  assert(data_x.size() == data_y.size() &&
         "Size of data x and y must matched!");
  assert(data_x.size() >= 2 && "Must have at least 2 points!");
  interpolation_detail::assertInsideDomain(data_x, x, "interpolateLinear");

  for (Index i{0}; i < data_y.size() - 1; ++i) {
    if (data_x[i] <= x && x <= data_x[i + 1]) {
      T a_i = (data_y[i + 1] - data_y[i]) / (data_x[i + 1] - data_x[i]);
      return a_i * (x - data_x[i]) + data_y[i];
    }
  }
  throw std::out_of_range("interpolateLinear: x is out of data_x domain");
}

template <Index n, typename T>
T interpolateLinear(const StaticVector<double, n> &data_x,
                    const StaticVector<double, n> &data_y, T x) {
  interpolation_detail::assertInsideDomain(data_x, x, "interpolateLinear");

  for (Index i{0}; i < data_y.size() - 1; ++i) {
    if (data_x[i] <= x && x <= data_x[i + 1]) {
      T a_i = (data_y[i + 1] - data_y[i]) / (data_x[i + 1] - data_x[i]);
      return a_i * (x - data_x[i]) + data_y[i];
    }
  }
  throw std::out_of_range("interpolateLinear: x is out of data_x domain");
}

// Cubic spline interpolation: ai(x-x_i)^3 + b_i(x-x_i)^2 + c_i(x-x_i) + d_i
// Continuous in both 1st and 2nd derivative
// Used for fitting curve
template <Index nPoints, typename T>
T interpolateCubicSpline(const StaticVector<double, nPoints> &data_x,
                         const StaticVector<double, nPoints> &data_y, T x) {
  static_assert(nPoints >= 3, "Must have at least 3 points!");
  interpolation_detail::assertInsideDomain(data_x, x, "interpolateCubicSpline");

  if constexpr (nPoints <= 4) {
    // Knot-a-not constraints become singular for only 3 points.
    return interpolatePolynomial(data_x, data_y, x);
  }
  constexpr Index n = nPoints - 1; //(n+1) points -> n gap
  StaticVector<double, n> h{}, eta{};
  for (Index i{0}; i < n; ++i) {
    h[i] = data_x[i + 1] - data_x[i];
    eta[i] = data_y[i + 1] - data_y[i];
  }
  Matrix<double, nPoints, 1> rhs{};
  for (Index i{1}; i < n; ++i) {
    rhs(i, 0) = eta[i] / h[i] - eta[i - 1] / h[i - 1];
  }
  Matrix<double, nPoints, nPoints> H{};
  // knot-a-not condition:
  // First row:
  H(0, 0) = h[1];
  H(0, 1) = -(h[0] + h[1]);
  H(0, 2) = h[0];
  // Last row:
  H(n, n - 2) = h[n - 1];
  H(n, n - 1) = -(h[n - 2] + h[n - 1]);
  H(n, n) = h[n - 2];
  for (Index r{1}; r < n; ++r) {
    Index i = r - 1;
    H(r, i) = 1.0 / 3.0 * h[i];
    H(r, i + 1) = 2.0 / 3.0 * (h[i] + h[i + 1]);
    H(r, i + 2) = 1.0 / 3.0 * h[i + 1];
  }
  Matrix<double, nPoints, 1> b = solveLinearSystem(H, rhs);
  // Note : component b_n is just a technical term
  StaticVector<double, n> a{}, c{};
  for (Index i{0}; i < n; ++i) {
    a[i] = 1.0 / (3.0 * h[i]) * (b[i + 1] - b[i]);
    c[i] = eta[i] / h[i] - 1.0 / 3.0 * h[i] * (b[i + 1] + 2.0 * b[i]);
  }
  for (Index i{0}; i < n; ++i) {
    if (data_x[i] <= x && x <= data_x[i + 1]) {
      T result = a[i] * std::pow(x - data_x[i], 3) +
                 b[i] * std::pow(x - data_x[i], 2) + c[i] * (x - data_x[i]) +
                 data_y[i];
      return result;
    }
  }
  throw std::out_of_range("interpolateCubicSpline: x is out of data_x domain");
}

#endif

#if 0
int main() {

  DynamicVector<double> data_x{0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

  DynamicVector<double> data_y{0.0, 0.8415, 0.9093, 0.1411, -0.7568, -0.9589,
  -0.2794};

  std::cout << std::fixed << std::setprecision(6);

  for (double x = 0.0; x <= 6.0; x += 0.1) {
    double fx = interpolatePolynomial(data_x, data_y, x);
    std::cout << "P(" << x << ") = " << fx << '\n';
  }

  return 0;
}
#endif
