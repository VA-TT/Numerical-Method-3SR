#ifndef PARENT_ELEMENT_SHAPE_FUNCTION_H
#define PARENT_ELEMENT_SHAPE_FUNCTION_H

#include "DualDifferentiation.h"
#include "Matrix.h"
#include "Vector.h"
#include "gaussQuadrature.h"
#include "interpolate.h"
#include <cassert>
#include <functional>
#include <iomanip>
#include <iostream>

// Range of 1D parent element Sr: [-1,1] -> mapping to [x1,x2] physical element
// First order polynomial shape function using Lagrange's basis in range [-1,1]
// Element node ordering (left to right):
//   (-1) n1 ------ n2 (1)
template <typename T> struct ElementL2 {
  // Variables
  T x1{}, x2{};     // real coordinates
  Index n1{}, n2 {} // nodes' ID in Mesh

  // Jacobian: J = dx/d(xi) = (x2-x1)/2
  T jacobian{}, midPoint{};

  // Mapping deriviation: dN/dx = dN/d(xi) * d(xi)/dx = dN/d(xi) * 1/J
  T dN1dx{}, dN2dx{};

  // Ctor
  ElementL2(T p1, T p2)
      : x1{p1}, x2{p2}, jacobian{(x2 - x1) / T{2}}, midPoint{(x2 + x1) / T{2}},
        dN1dx{-T{0.5} / jacobian}, dN2dx{T{0.5} / jacobian} {
    assert(x2 > x1 && "Element coordinates must satisfy x2 > x1");
  }

  // Functions
  //   Discretizing x = phi(xi) = x1 N1(xi) + x2 N2(xi) = J*xi + M = xi *
  //   (x2-x1)/2 + (x2+x1)/2
  static T N1_ref(const T &xi) const { return (1 - xi) / 2; }
  static T N2_ref(const T &xi) const { return (1 + xi) / 2; }

  // Function Phi(xi) used to map and its inverse function
  static T physicCoor(const T &xi) const {
    return N1_ref(xi) * x1 + N2_ref(xi) * x2;
  }
  static T parentCoor(const T &x) const { return (x - midPoint) / jacobian; }

  // Interpolation: int^(x2)_(x1) f(x)dx = int^1_(-1) f(phi(xi)) J d(xi)
  using namespace gaussQuadrature;
  inline double integrationGauss1D_ref(std::function<T(T)> f_xi, int n = 2) {

    const double *xi_ptr = nullptr;
    const double *w_ptr = nullptr;
    valueXiW(n, xi_ptr, w_ptr);
    double I = 0.0;
    for (int i = 0; i < n; ++i) {
      I += jacobian * w_ptr[i] * f_xi(xi_ptr[i]);
    }
    return I;
  }
};

//  Q4 2D rectangle parent element: [-1,-1] to [1,1]
// Element node ordering (counterclockwise from bottom-left):
//   (-1,1) n4 ------ n3 (1,1)
//          |          |
//          |          |
//  (-1,-1) n1 ------ n2 (1,-1)
// Use auto parameters to accept both double and Dual types

template <typename T> struct ElementQ4 {
  // Variables

  StaticVector<T, 4> x_nodes{}, y_nodes{}; // real coordinates
  T xi{}, T eta{};                         // reference coordinate
  Index n1{}, n2{}, n3{}, n4{};            // nodes' ID in Mesh

  // Jacobian: J = dx/d(xi) = (x2-x1)/2
  T jacobian{}, midPoint{};

  // Mapping deriviation: dN/dx = dN/d(xi) * d(xi)/dx = dN/d(xi) * 1/J
  T dN1dx{}, dN2dx{};

  // Ctor

  ElementL2(T p1, T p1)
      : x1{p1}, x2{p2}, jacobian{(x2 - x1) / T{2}}, midPoint{(x2 + x1) / T{2}},
        dN1dx{-T{0.5} / jacobian}, dN2dx{T{0.5} / jacobian} {
    assert(x2 > x1 && "Element coordinates must satisfy x2 > x1");
  }

  // Functions

  //   Discretizing x = phi(xi) = x1 N1(xi) + x2 N2(xi) = J*xi + M = xi *
  //   (x2-x1)/2 + (x2+x1)/2
  static T N1_ref(const T &xi, const T &eta) const {
    return T{0.25} * (T{1} - xi) * (T{1} - eta);
  }
  T static N2_ref(const T &xi, const T &eta) const {
    return T{0.25} * (T{1} + xi) * (T{1} - eta);
  }
  T static N3_ref(const T &xi, const T &eta) const {
    return T{0.25} * (T{1} + xi) * (T{1} + eta);
  }
  T static N4_ref(const T &xi, const T &eta) const {
    return T{0.25} * (T{1} - xi) * (T{1} + eta);
  }
  static StaticVector<T, 4> N_Q4(const T &xi, const T &eta) const {
    return {N1_ref(xi, eta), N2_ref(xi, eta), N3_ref(xi, eta), N4_ref(xi, eta)};
  }

  // Derivatives of shape functions with respect to xi
  static T dN1_dxi(auto xi, auto eta) const { return T{-0.25} * (T{1} - eta); }
  static T dN2_dxi(auto xi, auto eta) const { return T{0.25} * (T{1} - eta); }
  static T dN3_dxi(auto xi, auto eta) const { return T{0.25} * (T{1} + eta); }
  static T dN4_dxi(auto xi, auto eta) const { return T{-0.25} * (T{1} + eta); }

  // Derivatives of shape functions with respect to eta
  static T dN1_deta(auto xi, auto eta) const { return T{-0.25} * (T{1} - xi); }
  static T dN2_deta(auto xi, auto eta) const { return T{-0.25} * (T{1} + xi); }
  static T dN3_deta(auto xi, auto eta) const { return T{0.25} * (T{1} + xi); }
  static T dN4_deta(auto xi, auto eta) const { return T{0.25} * (T{1} - xi); }

  // Mapping from parent [-1,1]x[-1,1] to physical coordinates
  // x = sum(N_i * x_i)
  StaticVector<T, 2> physicCoor(T xi, T eta) {
    // Evaluate all shape functions at (xi, eta)
    StaticVector<T, 4> N{};
    N[0] = N1_ref(xi, eta);
    N[1] = N2_ref(xi, eta);
    N[2] = N3_ref(xi, eta);
    N[3] = N4_ref(xi, eta);

    T x = dotProduct(N, x_nodes); // dot product
    T y = dotProduct(N, y_nodes); // dot product
    return {x, y};
  }

  // Inverse mapping: (x,y) -> (xi,eta) - Optimized for rectangular elements
  // For general quadrilaterals, falls back to Newton-Raphson
  template <typename T>
  std::pair<T, T> parentCoor(T x, T y, int maxIter = 20, T tol = 1e-10) {
    // General case: Newton-Raphson iteration for non-rectangular elements
    T xi = 0.0, eta = 0.0; // Initial guess: center
    for (int iter = 0; iter < maxIter; ++iter) {
      // Compute residual: R = [x - x(xi,eta), y - y(xi,eta)]
      auto [x_curr, y_curr] = physicCoor(xi, eta);
      T Rx = x - x_curr;
      T Ry = y - y_curr;

      // Check convergence
      if (std::abs(Rx) < tol && std::abs(Ry) < tol) {
        return {xi, eta};
      }

      // Compute Jacobian
      Matrix<T, 2, 2> J = jacobian(xi, eta);

      // Solve J * delta = R for delta = [dxi, deta]
      T detJ = det(J);
      T dxi = (J(1, 1) * Rx - J(0, 1) * Ry) / detJ;
      T deta = (-J(1, 0) * Rx + J(0, 0) * Ry) / detJ;

      // Update
      xi += dxi;
      eta += deta;
    }

    // If not converged, throw warning or return current estimate
    std::cerr << "Warning: parentCoor2D did not converge!" << '\n';
    return {xi, eta};
  }

  // Jacobian matrix for 2D element: J = [dx/dxi  dx/deta]
  //                                     [dy/dxi  dy/deta]
  template <typename T> Matrix<T, 2, 2> jacobian(T xi, T eta) {

    Matrix<T, 2, 2> J;
    // Evaluate all derivatives at (xi, eta)
    StaticVector<T, 4> dN_dxi{}, dN_deta{};
    dN_dxi[0] = dN1_dxi(xi, eta);
    dN_dxi[1] = dN2_dxi(xi, eta);
    dN_dxi[2] = dN3_dxi(xi, eta);
    dN_dxi[3] = dN4_dxi(xi, eta);

    dN_deta[0] = dN1_deta(xi, eta);
    dN_deta[1] = dN2_deta(xi, eta);
    dN_deta[2] = dN3_deta(xi, eta);
    dN_deta[3] = dN4_deta(xi, eta);

    J(0, 0) = dotProduct(dN_dxi, x_nodes);  // dx/dxi
    J(0, 1) = dotProduct(dN_deta, x_nodes); // dx/deta
    J(1, 0) = dotProduct(dN_dxi, y_nodes);  // dy/dxi
    J(1, 1) = dotProduct(dN_deta, y_nodes); // dy/deta

    // Check for non-singular Jacobian
    assert(
        det(J) > 0.0 &&
        "Negative Jacobian - element is inverted or has wrong node ordering!");
    return J;
  }

  // 2D Gauss quadrature integration on reference element [-1,1]x[-1,1]
  // Computes: int_{-1}^{1} int_{-1}^{1} f(xi,eta) * det(J) dxi deta
  inline double integrationGauss2D_ref(std::function<T(T, T)> f_xi_eta,
                                       int n = 2) {
    const double *xi_ptr = nullptr;
    const double *w_ptr = nullptr;
    valueXiW(n, xi_ptr, w_ptr);

    double I = 0.0;

    // Double loop over Gauss points
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        double xi = xi_ptr[i];
        double eta = xi_ptr[j];

        // Compute Jacobian determinant at (xi, eta)
        Matrix<double, 2, 2> J = jacobian(xi, eta);

        // Accumulate weighted sum
        I += w_ptr[i] * w_ptr[j] * f_xi_eta(xi, eta) * det(J);
      }
    }

    return I;
  }

  // Compute derivatives of shape functions in physical coordinates
  std::pair<StaticVector<T, 4>, StaticVector<T, 4>> gradientN(T xi, T eta) {

    // Compute Jacobian and its inverse
    Matrix<T, 2, 2> J_inv = Jacobian(xi, eta).inverse();

    // Evaluate derivatives in parent coordinates (2x4 matrix)
    // [dN/dxi ] = [dN1/dxi  dN2/dxi  dN3/dxi  dN4/dxi ]
    // [dN/deta]   [dN1/deta dN2/deta dN3/deta dN4/deta]
    Matrix<T, 2, 4> dN_parent;
    dN_parent(0, 0) = dN1_dxi(xi, eta);
    dN_parent(0, 1) = dN2_dxi(xi, eta);
    dN_parent(0, 2) = dN3_dxi(xi, eta);
    dN_parent(0, 3) = dN4_dxi(xi, eta);

    dN_parent(1, 0) = dN1_deta(xi, eta);
    dN_parent(1, 1) = dN2_deta(xi, eta);
    dN_parent(1, 2) = dN3_deta(xi, eta);
    dN_parent(1, 3) = dN4_deta(xi, eta);

    // Transform to physical coordinates: J_inv * dN_parent = (2x2) * (2x4) =
    // (2x4) [dN/dx] = J_inv * [dN/dxi ] [dN/dy]           [dN/deta]
    Matrix<T, 2, 4> dN_physical = J_inv * dN_parent;

    // Extract results
    DynamicVector<T> dN_dx(4), dN_dy(4);
    for (int i = 0; i < 4; ++i) {
      dN_dx[i] = dN_physical(0, i);
      dN_dy[i] = dN_physical(1, i);
    }

    return {dN_dx, dN_dy};
  }
};

// void shapeFunction() {
//   // Equally divied
//   int local_i = 0;
//   int local_j = 1;
//   for (Index e{0}; e < nElements; ++e) {
//     int i = eleOrigin[e];
//     int j = eleEnd[e];
//     double x_i = nodes[i];
//     double x_j = nodes[j];
//     element[e] = {x_i, x_j};
//     length[e] = constexpr_fabs(x_j - x_i);
//     k[e] = EA / length[e];

//     // MUST set parameter "x" here to type <<auto>> in order to accept Dual
//     // class as input for derivative calculating
//     auto shapefunction_i = [=](auto x) {
//       return basisLagrange(local_i, element[e], x);
//     };
//     auto dShape_i = [=](double x) { return automaticDiff(shapefunction_i, x);
//     }; auto shapefunction_j = [=](auto x) {
//       return basisLagrange(local_j, element[e], x);
//     };
//     auto dShape_j = [=](double x) { return automaticDiff(shapefunction_j, x);
//     }; for (Index kk{0}; kk < nNodes; ++kk) {
//       if (kk == i) {
//         N[kk].push_back(shapefunction_i);
//         N_x[kk].push_back(dShape_i);
//       } else if (kk == j) {
//         N[kk].push_back(shapefunction_j);
//         N_x[kk].push_back(dShape_j);
//       } else {
//         N[kk].push_back([](double x) { return 0.0; });
//         N_x[kk].push_back([](double x) { return 0.0; });
//       }
//     }
//   }
// }

inline auto cubicBSpline = [](auto x, auto h) {
  if (-2.0 * h <= x <= -h)
    return (x * x * x) / (6.0 * h * h * h) + x * x / (h * h) + 2.0 * x / h +
           4.0 / 3.0;
  if (-h <= x <= 0.0)
    return (-x * x * x) / (2.0 * h * h * h) - x * x / (h * h) + 4.0 / 3.0;
  if (-2.0 * h <= x <= -h)
    return (x * x * x) / (2.0 * h * h * h) - x * x / (h * h) + 4.0 / 3.0;
  if (-2.0 * h <= x <= -h)
    return (-x * x * x) / (6.0 * h * h * h) + x * x / (h * h) - 2.0 * x / h +
           4.0 / 3.0;
  else
    return 0;
};

#endif
