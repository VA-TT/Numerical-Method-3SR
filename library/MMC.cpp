#ifndef CONTINUUM_MECHANICS_H
#define CONTINUUM_MECHANICS_H

#include "Matrix.h"
#include "Vector.h"
#include "parentElement.h"

using VectorFunction3D = std::function<StaticVector<Dual, 3>(Dual, Dual, Dual)>;

template <typename T, Index n> struct box {
  StaticVector<StaticVector<T, n>, n> equationMotion{};
  Matrix<T, n, n> deformGrad{};
  defomationGradient() {}
  rightStretchTensor() {}
  leftStretchTensor() {}
};

#endif