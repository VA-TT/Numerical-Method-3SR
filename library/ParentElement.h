#ifndef PARENT_ELEMENT_SHAPE_FUNCTION_H
#define PARENT_ELEMENT_SHAPE_FUNCTION_H

#include "DualDifferentiation.h"
#include "Matrix.h"
#include "Particle-Node.h"
#include "Vector.h"
#include "gaussQuadrature.h"
#include "interpolate.h"
#include <cassert>
#include <functional>
#include <iostream>
#include <span>

// Linear mapping of 1D and 2D line elements, interpomation.h contains more type
// of interpolation

//  Range of 1D parent element Sr: [-1,1] -> mapping to [x1,x2] physical element
//  First order polynomial shape function using Lagrange's basis in range [-1,1]
//  Element node ordering (left to right):
//    (-1) n1 ------ n2 (1)
template <typename T> struct ElementL2 {
  // Variables
  Index n1{}, n2{};                   // nodes' ID in Mesh
  std::span<const Node1D<T>> nodes{}; // reference to nodal array
  bool isActive{false};
  // Ctor
  ElementL2(Index id1, Index id2, std::span<const Node1D<T>> ptn)
      : n1{id1}, n2{id2}, nodes{ptn} {
    assert(nodes[n1].pos < nodes[n2].pos &&
           "Element coordinates must satisfy x2 > x1");
  }
  // Other defaults
  ElementL2() = default;
  ElementL2(const ElementL2 &) = default;
  ElementL2(ElementL2 &&) = default;
  ElementL2 &operator=(const ElementL2 &) = default;
  ElementL2 &operator=(ElementL2 &&) = default;
  ~ElementL2() = default;

  T x1() const { return nodes[n1].pos; }
  T x2() const { return nodes[n2].pos; }
  // Jacobian: J = dx/d(xi) = (x2-x1)/2 = const
  T jacobian() const { return (x2() - x1()) / T{2}; }

  // Midpoint
  T midPoint() const { return (x2() + x1()) / T{2}; }

  // Shape function derivatives w.r.t x
  T dN1_dx() const { return -T{0.5} / jacobian(); }
  T dN2_dx() const { return T{0.5} / jacobian(); }

  T N1_ref(T xi) const { return (T{1} - xi) / T{2}; }
  T N2_ref(T xi) const { return (T{1} + xi) / T{2}; }

  // Function Phi(xi) used to map and its inverse function
  //   Discretizing x = phi(xi) = x1 N1(xi) + x2 N2(xi) = J*xi + M = xi *
  //   (x2-x1)/2 + (x2+x1)/2
  T physicCoord(T xi) const { return N1_ref(xi) * x1() + N2_ref(xi) * x2(); }
  T parentCoord(T x) const { return (x - midPoint()) / jacobian(); }

  // Interpolation: int^(x2)_(x1) f(x)dx = int^1_(-1) f(phi(xi)) J d(xi)

  template <typename Func>
  T integrationGauss1D_ref(Func f_xi, int n = 2) const {
    const double *xi_ptr = nullptr;
    const double *w_ptr = nullptr;
    gaussQuadrature::valueXiW(n, xi_ptr, w_ptr);
    T I{0.0};
    T J = jacobian(); // compute once per integration
    for (int i = 0; i < n; ++i) {
      I += J * static_cast<T>(w_ptr[i]) * f_xi(static_cast<T>(xi_ptr[i]));
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

template <typename T> struct ElementQ4 {

  // Variables
  Index n1{}, n2{}, n3{}, n4{};       // nodes' ID in Mesh
  std::span<const Node2D<T>> nodes{}; // reference to particle array
  bool isActive{false};

  ElementQ4(Index id1, Index id2, Index id3, Index id4,
            std::span<const Node2D<T>> ptn)
      : n1{id1}, n2{id2}, n3{id3}, n4{id4}, nodes{ptn} {}
  // Other defaults
  ElementQ4() = default;
  ElementQ4(const ElementQ4 &) = default;
  ElementQ4(ElementQ4 &&) = default;
  ElementQ4 &operator=(const ElementQ4 &) = default;
  ElementQ4 &operator=(ElementQ4 &&) = default;
  ~ElementQ4() = default;

  // Get x and y coordinate vectors from particle positions
  StaticVector<T, 4> getX_nodes() const {
    return {nodes[n1].pos.x(), nodes[n2].pos.x(), nodes[n3].pos.x(),
            nodes[n4].pos.x()};
  }

  StaticVector<T, 4> getY_nodes() const {
    return {nodes[n1].pos.y(), nodes[n2].pos.y(), nodes[n3].pos.y(),
            nodes[n4].pos.y()};
  }

  // Functions
  T N1_ref(T xi, T eta) const { return T{0.25} * (T{1} - xi) * (T{1} - eta); }
  T N2_ref(T xi, T eta) const { return T{0.25} * (T{1} + xi) * (T{1} - eta); }
  T N3_ref(T xi, T eta) const { return T{0.25} * (T{1} + xi) * (T{1} + eta); }
  T N4_ref(T xi, T eta) const { return T{0.25} * (T{1} - xi) * (T{1} + eta); }
  StaticVector<T, 4> N_Q4(T xi, T eta) const {
    return {N1_ref(xi, eta), N2_ref(xi, eta), N3_ref(xi, eta), N4_ref(xi, eta)};
  }

  // Derivatives of shape functions with respect to xi
  T dN1_dxi(auto xi, auto eta) const { return T{-0.25} * (T{1} - eta); }
  T dN2_dxi(auto xi, auto eta) const { return T{0.25} * (T{1} - eta); }
  T dN3_dxi(auto xi, auto eta) const { return T{0.25} * (T{1} + eta); }
  T dN4_dxi(auto xi, auto eta) const { return T{-0.25} * (T{1} + eta); }

  // Derivatives of shape functions with respect to eta
  T dN1_deta(auto xi, auto eta) const { return T{-0.25} * (T{1} - xi); }
  T dN2_deta(auto xi, auto eta) const { return T{-0.25} * (T{1} + xi); }
  T dN3_deta(auto xi, auto eta) const { return T{0.25} * (T{1} + xi); }
  T dN4_deta(auto xi, auto eta) const { return T{0.25} * (T{1} - xi); }

  StaticVector<T, 4> dN_dxi(T xi, T eta) const {
    return {dN1_dxi(xi, eta), dN2_dxi(xi, eta), dN3_dxi(xi, eta),
            dN4_dxi(xi, eta)};
  }
  StaticVector<T, 4> dN_deta(T xi, T eta) const {
    return {dN1_deta(xi, eta), dN2_deta(xi, eta), dN3_deta(xi, eta),
            dN4_deta(xi, eta)};
  }

  // Jacobian matrix for 2D element: J = [dx/dxi  dx/deta]
  //                                     [dy/dxi  dy/deta]
  Matrix<T, 2, 2> jacobian(T xi, T eta) const {

    StaticVector<T, 4> x_nodes = getX_nodes();
    StaticVector<T, 4> y_nodes = getY_nodes();

    Matrix<T, 2, 2> J;
    J(0, 0) = dotProduct(dN_dxi(xi, eta), x_nodes);  // dx/dxi
    J(0, 1) = dotProduct(dN_deta(xi, eta), x_nodes); // dx/deta
    J(1, 0) = dotProduct(dN_dxi(xi, eta), y_nodes);  // dy/dxi
    J(1, 1) = dotProduct(dN_deta(xi, eta), y_nodes); // dy/deta

    // Check for non-singular Jacobian
    double detTol = 1e-12;
    if (std::abs(det(J)) <= detTol) {
      throw std::runtime_error("Negative Jacobian - element is inverted or has "
                               "wrong node ordering!");
    }
    return J;
  }

  // Mapping from parent [-1,1]x[-1,1] to physical coordinates
  // x = sum(N_i * x_i)
  StaticVector<T, 2> physicCoord(T xi, T eta) const {
    StaticVector<T, 4> x_nodes = getX_nodes();
    StaticVector<T, 4> y_nodes = getY_nodes();
    StaticVector<T, 4> N{N_Q4(xi, eta)};
    StaticVector<T, 2> physic{};
    physic.x() = dotProduct(N, x_nodes);
    physic.y() = dotProduct(N, y_nodes);
    return physic;
  }

  // Inverse mapping: (x,y) -> (xi,eta)
  // Newton-Raphson method for general quadrilaterals
  StaticVector<T, 2> parentCoord(T x, T y, int maxIter = 20,
                                 T tol = 1e-10) const {
    auto absValue = [](const T &v) {
      if constexpr (requires(const T &a) { a.getVal(); }) {
        return std::abs(v.getVal());
      } else {
        return std::abs(v);
      }
    };

    // General case: Newton-Raphson iteration for non-rectangular elements
    T xi = 0.0, eta = 0.0; // Initial guess: center
    for (int iter = 0; iter < maxIter; ++iter) {
      // Compute residual: R = [x - x(xi,eta), y - y(xi,eta)]
      auto curr = physicCoord(xi, eta);
      T x_curr = curr.x();
      T y_curr = curr.y();
      T Rx = x - x_curr;
      T Ry = y - y_curr;

      // Check convergence
      if (absValue(Rx) < absValue(tol) && absValue(Ry) < absValue(tol)) {
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
    std::cerr << "Warning: parentCoord2D did not converge!" << '\n';
    return {xi, eta};
  }

  // 2D Gauss quadrature integration on reference element [-1,1]x[-1,1]
  // Computes: int_{-1}^{1} int_{-1}^{1} f(xi,eta) * det(J) dxi deta
  template <typename Func>
  T integrationGauss2D_ref(Func f_xi_eta, int n = 2) const {
    const double *xi_ptr = nullptr;
    const double *w_ptr = nullptr;
    gaussQuadrature::valueXiW(n, xi_ptr, w_ptr);

    T I{0.0};

    // Double loop over Gauss points
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        T xi = static_cast<T>(xi_ptr[i]);
        T eta = static_cast<T>(xi_ptr[j]);

        // Compute Jacobian determinant at (xi, eta)
        Matrix<T, 2, 2> J = jacobian(xi, eta);

        // Accumulate weighted sum
        I += static_cast<T>(w_ptr[i]) * static_cast<T>(w_ptr[j]) *
             f_xi_eta(xi, eta) * det(J);
      }
    }

    return I;
  }

  // Compute derivatives of shape functions in physical coordinates
  StaticVector<StaticVector<T, 4>, 2> gradientN(T xi, T eta) const {

    // Compute Jacobian and its inverse
    Matrix<T, 2, 2> J_inv = jacobian(xi, eta).inverse();

    // Evaluate derivatives in parent coordinates (2x4 matrix)
    // [dN/dxi ] = [dN1/dxi  dN2/dxi  dN3/dxi  dN4/dxi ]
    // [dN/deta]   [dN1/deta dN2/deta dN3/deta dN4/deta]
    Matrix<T, 2, 4> dN_parent;
    StaticVector<T, 4> dNxi = dN_dxi(xi, eta);
    StaticVector<T, 4> dNeta = dN_deta(xi, eta);
    for (int i = 0; i < 4; ++i) {
      dN_parent(0, i) = dNxi[i];
      dN_parent(1, i) = dNeta[i];
    }

    // Transform to physical coordinates:
    Matrix<T, 2, 4> dN_physical = J_inv * dN_parent;

    // Extract results
    StaticVector<T, 4> dN_dx{}, dN_dy{};
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

#endif
