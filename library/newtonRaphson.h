#ifndef MY_NEWTON_RAPHSON_ROOT_H
#define MY_NEWTON_RAPHSON_ROOT_H

#include "DualDifferentiation.h"
#include "Matrix.h"
#include <algorithm>
#include <cmath>
#include <type_traits>
#include <vector>

// General polynomial function: coef[0] + coef[1]*x + coef[2]*x^2 + ...
// Works with both double and Dual types
template <typename T> auto polyFunc(const std::vector<double> &coef, T x) {

  T result{};
  if constexpr (std::is_same_v<T, Dual>) {
    result = Dual{0.0, 0.0};
  } else {
    result = T{0.0};
  }

  T x_power{};
  if constexpr (std::is_same_v<T, Dual>) {
    x_power = Dual{1.0, 0.0};
  } else {
    x_power = T{1.0};
  }

  for (const auto &c : coef) {
    result = result + c * x_power;
    x_power = x_power * x;
  }
  return result;
}

// Generic derivative for any polynomial using Dual automatic differentiation
inline double polyDerivative(const std::vector<double> &coef, double x) {
  auto func = [&coef](Dual d) { return polyFunc(coef, d); };
  return automaticDiff(func, x);
}

// Generic Newton-Raphson for any function with its derivative
// func: function to find root of
// dfunc: derivative of the function
template <typename Func, typename DFunc>
inline double newtonRaphsonGeneric(Func func, DFunc dfunc, double x0,
                                   double epsilon = 1e-10,
                                   int max_iter = 1000) {
  double x_n = x0;
  const double tiny = 1e-16;

  for (int i = 0; i < max_iter; ++i) {
    double f_n = func(x_n);
    double f_prime = dfunc(x_n);

    if (std::abs(f_prime) < tiny)
      break;

    double x_next = x_n - f_n / f_prime;

    if (std::abs(x_next - x_n) <= epsilon)
      return x_next;

    x_n = x_next;
  }
  return x_n;
}

// Generic Newton-Raphson for any polynomial
inline double newtonRaphsonPoly(const std::vector<double> &coef, double x0,
                                double epsilon = 1e-10, int max_iter = 1000) {
  double x_n = x0;
  const double tiny = 1e-16;

  for (int i = 0; i < max_iter; ++i) {
    double f_n = polyFunc(coef, x_n);
    double f_prime = polyDerivative(coef, x_n);

    if (std::abs(f_prime) < tiny)
      break;

    double x_next = x_n - f_n / f_prime;

    if (std::abs(x_next - x_n) <= epsilon)
      return x_next;

    x_n = x_next;
  }
  return x_n;
}

// ============= POLYNOMIAL DEFLATION =============

// Synthetic division: divide polynomial by (x - root)
// Returns deflated coefficients (degree reduced by 1)
inline std::vector<double> syntheticDivision(const std::vector<double> &coef,
                                             double root) {
  if (coef.size() <= 1)
    return {};

  std::vector<double> deflated(coef.size() - 1);
  deflated[coef.size() - 2] = coef[coef.size() - 1]; // highest degree term

  // Horner's method for synthetic division
  for (int i = static_cast<int>(coef.size()) - 3; i >= 0; --i) {
    deflated[i] = coef[i + 1] + root * deflated[i + 1];
  }

  return deflated;
}

// Find all real roots using polynomial deflation
// coef: coefficients [a0, a1, a2, ...] for a0 + a1*x + a2*x^2 + ...
// initial_guesses: starting points for Newton-Raphson
inline std::vector<double>
findAllRoots(const std::vector<double> &coef,
             const std::vector<double> &initial_guesses,
             double epsilon = 1e-10) {
  std::vector<double> roots;
  std::vector<double> current_coef = coef;

  for (const auto &x0 : initial_guesses) {
    if (current_coef.size() <= 1)
      break; // no more roots to find

    // Find one root
    double root = newtonRaphsonPoly(current_coef, x0, epsilon);

    // Verify it's actually a root
    double f_root = polyFunc(current_coef, root);
    if (std::abs(f_root) < 1e-6) {
      roots.push_back(root);

      // Deflate: divide by (x - root)
      current_coef = syntheticDivision(current_coef, root);
    }
  }

  // Remove duplicates (roots found multiple times)
  std::sort(roots.begin(), roots.end());
  roots.erase(
      std::unique(roots.begin(), roots.end(),
                  [](double a, double b) { return std::abs(a - b) < 1e-8; }),
      roots.end());

  return roots;
}

// Generic Newton-Raphson for 2D nonlinear systems:
template <typename T, typename ResidualFunc, typename JacobianFunc>
inline StaticVector<T, 2>
newtonRaphsonSystemEquations(ResidualFunc residual, JacobianFunc jacobian,
                             const StaticVector<T, 2> &initValue,
                             int maxIter = 20, T tol = static_cast<T>(1e-10)) {
  auto absValue = [](const T &v) {
    if constexpr (requires(const T &a) { a.getVal(); }) {
      return std::abs(v.getVal());
    } else {
      return std::abs(v);
    }
  };

  StaticVector<T, 2> X{initValue};
  const T tiny = static_cast<T>(1e-16);

  for (int iter = 0; iter < maxIter; ++iter) {
    const StaticVector<T, 2> R = residual(X);

    // Converged in residual norm
    if (absValue(R.x()) < absValue(tol) && absValue(R.y()) < absValue(tol)) {
      return X;
    }

    const Matrix<T, 2, 2> J = jacobian(X);
    const T detJ = det(J);
    if (absValue(detJ) <= absValue(tiny)) {
      break;
    }

    // Newton step: J * dX = -R
    const StaticVector<T, 2> rhs{-R.x(), -R.y()};
    const StaticVector<T, 2> dX = solveLinearSystem(J, rhs);
    X += dX;

    // Converged in increment norm
    if (absValue(dX.x()) < absValue(tol) && absValue(dX.y()) < absValue(tol)) {
      return X;
    }
  }

  return X;
}

StaticVector<double, 2>
newtonRaphsonSystemEquations(const VectorFunction2D &u,
                             const StaticVector<double, 2> &initValue) {
  StaticVector<double, 2> X{initValue};

  StaticVector<double, 2> dX{};
  constexpr int maxIter = 100;

  for (int iter = 0; iter < maxIter; ++iter) {
    Matrix<double, 2, 2> J = Jacobian(u, X);

    auto uVal = u(X[0], X[1]);
    StaticVector<double, 2> rhs{-uVal.x().getVal(), -uVal.y().getVal()};

    dX = solveLinearSystem(J, rhs);
    X += dX;

    bool isSolved = true;
    for (const auto e : dX) {
      if (!approximatelyEqualAbsRel(e, 0.0)) {
        isSolved = false;
        break;
      }
    }

    if (isSolved)
      return X;
  }

  return X;
}

StaticVector<double, 3>
newtonRaphsonSystemEquations(const VectorFunction3D &u,
                             const StaticVector<double, 3> &initValue) {
  StaticVector<double, 3> X{initValue};

  StaticVector<double, 3> dX{};
  constexpr int maxIter = 100;

  for (int iter = 0; iter < maxIter; ++iter) {
    Matrix<double, 3, 3> J = Jacobian(u, X);

    auto uVal = u(X[0], X[1], X[2]);
    StaticVector<double, 3> rhs{-uVal.x().getVal(), -uVal.y().getVal(),
                                -uVal.z().getVal()};

    dX = solveLinearSystem(J, rhs);
    X += dX;

    bool isSolved = true;
    for (const auto e : dX) {
      if (!approximatelyEqualAbsRel(e, 0.0)) {
        isSolved = false;
        break;
      }
    }

    if (isSolved)
      return X;
  }

  return X;
}

#endif
