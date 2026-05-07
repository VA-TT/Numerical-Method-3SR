#ifndef SHAPE_FUNCTION_H
#define SHAPE_FUNCTION_H

#include "interpolate.h"
#include <cmath>

enum class BasisFaminly {
  hat,
  quadBSpline,
  cubicBSpline,
  bernstein,
  typeCount,
};

// ============================================================================
// HatShapeFunction1D<T>: Linear basis function (piecewise linear)
// Support: [-1, 1]
// N(xi) = max(0, 1 - |xi|)
// dN/dxi = -sign(xi) for |xi| < 1, 0 otherwise
// ============================================================================
template <typename T> struct HatShapeFunction1D {
  static constexpr int supportRadius = 1;

  // Basis function value
  static T N(T xi) {
    T abs_xi = std::abs(xi);
    if (abs_xi >= T{1})
      return T{0};
    return T{1} - abs_xi;
  }

  // Derivative dN/dxi
  static T dN(T xi) {
    T abs_xi = std::abs(xi);
    if (abs_xi >= T{1})
      return T{0}; // zero derivative outside support
    if (xi < T{0})
      return T{1}; // slope +1 on left side [-1, 0)
    else if (xi > T{0})
      return -T{1}; // slope -1 on right side (0, 1]
    else
      return T{0}; // undefined at xi=0, use 0 or average
  }
};

// ============================================================================
// BSplineShapeFunction1D<T>: Cubic B-spline basis function
// Support: [-2, 2]
// Piecewise cubic polynomial with C2 continuity
// ============================================================================
template <typename T> struct BSplineShapeFunction1D {
  static constexpr int supportRadius = 2;

  // Basis function value using standard cubic B-spline
  // N(r) = {
  //   2/3 - r^2 + |r|^3/2        for |r| < 1
  //   (2 - |r|)^3 / 6            for 1 <= |r| < 2
  //   0                            for |r| >= 2
  // }
  static T N(T r) {
    T abs_r = std::abs(r);

    if (abs_r < T{1}) {
      // Region: |r| < 1
      T r2 = abs_r * abs_r;
      T r3 = r2 * abs_r;
      return (T{2} / T{3}) - r2 + r3 / T{2};
    } else if (abs_r < T{2}) {
      // Region: 1 <= |r| < 2
      T temp = T{2} - abs_r;
      return temp * temp * temp / T{6};
    } else {
      // Region: |r| >= 2
      return T{0};
    }
  }

  // Derivative dN/dr
  // dN/dr = {
  //   -2r + (3/2)*r^2 * sign(r)   for |r| < 1
  //   -sign(r) * (2-|r|)^2 / 2    for 1 <= |r| < 2
  //   0                            for |r| >= 2
  // }
  static T dN(T r) {
    T abs_r = std::abs(r);
    T sign_r = (r < T{0}) ? -T{1} : (r > T{0}) ? T{1} : T{0};

    if (abs_r < T{1}) {
      // Region: |r| < 1
      T r2 = abs_r * abs_r;
      return sign_r * (-T{2} * abs_r + (T{3} / T{2}) * r2);
    } else if (abs_r < T{2}) {
      // Region: 1 <= |r| < 2
      T temp = T{2} - abs_r;
      return -sign_r * temp * temp / T{2};
    } else {
      // Region: |r| >= 2
      return T{0};
    }
  }
};





#endif