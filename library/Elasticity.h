#ifndef HOOKE_ELASTIC_H
#define HOOKE_ELASTIC_H

#include "Matrix.h"
#include <iostream>
#include <string_view>

template <typename T> Matrix<T, 6, 6> tangentOperator3D(T E, T v) {
  // 3D elasticity (6x6) isotropic Hooke's law
  T c0 = E / ((1.0 + v) * (1.0 - 2.0 * v)); // factor for normal components
  T c1 = v;
  T c2 = 1.0 - v;
  T G = E / (2.0 * (1.0 + v)); // shear modulus
  // clang-format off
  return Matrix<T, 6, 6>{c2 * c0, c1 * c0, c1 * c0, 0.0, 0.0, 0.0,
                         c1 * c0, c2 * c0, c1 * c0, 0.0, 0.0, 0.0,
                         c1 * c0, c1 * c0, c2 * c0, 0.0, 0.0, 0.0,
                         0.0,     0.0,     0.0,     G,   0.0, 0.0,
                         0.0,     0.0,     0.0,     0.0, G,   0.0,
                         0.0,     0.0,     0.0,     0.0, 0.0, G};
  // clang-format on
}
// Plane 2D
template <typename T>
Matrix<T, 3, 3> tangentOperator2D(std::string_view type, T E, T v) {
  if (type == "planeStress") {
    return E / (1.0 - v * v) *
           Matrix<T, 3, 3>{1.0, v, 0.0, v, 1.0, 0.0, 0.0, 0.0, (1.0 - v) / 2.0};
  } else if (type == "planeStrain") {
    return E / ((1.0 + v) * (1.0 - 2.0 * v)) *
           Matrix<T, 3, 3>{1.0 - v, v,       0.0,
                           v,       1.0 - v, 0.0,
                           0.0,     0.0,     (1.0 - 2.0 * v) / 2.0};
  } else {
    std::cerr << "Error: Unknown elasticity type \"" << type << "\"\n";
    return Matrix<T, 3, 3>{};
  }
}

#endif // HOOKE_ELASTIC_H
