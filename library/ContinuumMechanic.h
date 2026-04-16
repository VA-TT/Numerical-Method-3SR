#ifndef CONTINUUM_MECHANICS_H
#define CONTINUUM_MECHANICS_H

#include "Coordinates.h"
#include "DualDifferentiation.h"
#include "Matrix.h"
#include "Vector.h"
#include "parentElement.h"
#include <functional>
#include <stdexcept>

// Continuum kinematics with linear motion map: x(X,t) = K(t) * X.
// K(t) is provided as a matrix-valued function of Dual time.
template <typename T, Index n> struct MediumCon {
  using MapMatrix = std::function<Matrix<Dual, n, n>(Dual)>;

  // Member variables
  // Reference configuration: X
  StaticVector<T, n> ref{};
  // Current configuration: x
  StaticVector<T, n> current{};
  // Lagrange's description:
  StaticVector<T, n> velocityL{}, accelerationL;
  // Euler's description of velocity:
  StaticVector<T, n> velocityE{}, accelerationE;

  MapMatrix deformMap{};

  bool hasRef{false};
  bool hasCurrent{false};

private:
  // The deformation gradient's value and its derivatives:
  static Matrix<T, n, n> valuePart(const Matrix<Dual, n, n> &mDual) {
    Matrix<T, n, n> out{};
    for (Index i = 0; i < mDual.length(); ++i) {
      out[i] = static_cast<T>(mDual[i].getVal());
    }
    return out;
  }

  static Matrix<T, n, n> derivativePart(const Matrix<Dual, n, n> &mDual) {
    Matrix<T, n, n> out{};
    for (Index i = 0; i < mDual.length(); ++i) {
      out[i] = static_cast<T>(mDual[i].getDer());
    }
    return out;
  }

  static Matrix<T, n, n> secondDerivativePart(const Matrix<Dual, n, n> &mDual) {
    Matrix<T, n, n> out{};
    for (Index i = 0; i < mDual.length(); ++i) {
      out[i] = static_cast<T>(mDual[i].getDer2());
    }
    return out;
  }

public:
  MediumCon() = default;
  // Material-description constructor
  MediumCon(const MapMatrix &motionEquation,
            const StaticVector<T, n> &refConfig)
      : ref(refConfig), deformMap(motionEquation), hasRef(true) {}

  // Spatial-description constructor
  MediumCon(const MapMatrix &inverseMotionEquation,
            const StaticVector<T, n> &currentConfig, bool isInverseMap)
      : current(currentConfig), hasCurrent(true) {
    if (!isInverseMap)
      throw std::invalid_argument(
          "Set isInverseMap=true when passing inverse motion equation.");

    deformMap = [inverseMotionEquation](Dual t) -> Matrix<Dual, n, n> {
      const Matrix<Dual, n, n> invDual = inverseMotionEquation(t);

      Matrix<T, n, n> invVal{};
      for (Index i = 0; i < invDual.length(); ++i) {
        invVal[i] = static_cast<T>(invDual[i].getVal());
      }

      const Matrix<T, n, n> kVal = invVal.inverse();

      Matrix<Dual, n, n> kDual{};
      for (Index i = 0; i < kDual.length(); ++i) {
        kDual[i] = Dual{static_cast<double>(kVal[i]), 0.0};
      }

      return kDual;
    };
  }

  Matrix<T, n, n> deformMap(T t) const {
    if (!deformMap)
      throw std::logic_error("deformMap is not set.");
    return valuePart(deformMap(Dual{static_cast<double>(t), 0.0}));
  }

  Matrix<T, n, n> inverseMap(T t) const { return deformMap(t).inverse(); }

  // dK/dt from Dual derivative part.
  Matrix<T, n, n> dKdt(T t) const {
    if (!deformMap)
      throw std::logic_error("deformMap is not set.");
    return derivativePart(deformMap(Dual{static_cast<double>(t), 1.0}));
  }

  // d2K/dt2 from Dual second derivative part.
  Matrix<T, n, n> d2Kdt2(T t) const {
    if (!deformMap)
      throw std::logic_error("deformMap is not set.");
    return secondDerivativePart(
        deformMap(Dual{static_cast<double>(t), 1.0, 0.0}));
  }

  // x_P(t) = K(t) * X
  StaticVector<T, n> currentPosition(T t) const {
    if (!hasRef)
      throw std::logic_error("Reference configuration not set.");
    return deformMap(t) * ref;
  }

  // X = K^{-1}(t) * x
  StaticVector<T, n> referencePosition(T t) const {
    if (!hasCurrent)
      throw std::logic_error("Current configuration not set.");
    return inverseMap(t) * current;
  }

  // Material velocity: v(X,t) = dx/dt = dK/dt * X.
  StaticVector<T, n> velocityMaterial(T t) const {
    if (!hasRef)
      throw std::logic_error("Reference configuration not set.");
    return dKdt(t) * ref;
  }

  StaticVector<T, n> displacementMaterial(T t) {};
  StaticVector<T, n> displacementSpatial(T t) {};

  StaticVector<T, n> accelerationMaterial(T t) {};
  // Material acceleration: a(X,t) = d2x/dt2 = d2K/dt2 * X.
  StaticVector<T, n> accelerationMaterial(T t) const {
    if (!hasRef)
      throw std::logic_error("Reference configuration not set.");
    return d2Kdt2(t) * ref;
  }

  // Euler velocity at any spatial point x:
  // v(x,t) = (dK/dt * K^{-1}) * x.
  StaticVector<T, n> velocityEulerAt(T t,
                                     const StaticVector<T, n> &xSpatial) const {
    const Matrix<T, n, n> L = dKdt(t) * inverseMap(t);
    return L * xSpatial;
  }

  // Euler velocity at stored current position.
  StaticVector<T, n> velocityEuler(T t) const {
    if (!hasCurrent)
      throw std::logic_error("Current configuration not set.");
    return velocityEulerAt(t, current);
  }

  // Euler acceleration at any spatial point x:
  // a = (dL/dt + L^2) x, where L = dK/dt * K^{-1}.
  StaticVector<T, n>
  accelerationEulerAt(T t, const StaticVector<T, n> &xSpatial) const {
    const Matrix<T, n, n> L = dKdt(t) * inverseMap(t);
    const Matrix<T, n, n> dLdt = d2Kdt2(t) * inverseMap(t) - L * L;

    return (dLdt + L * L) * xSpatial;
  }

  // Euler acceleration at stored current position.
  StaticVector<T, n> accelerationEuler(T t) const {
    if (!hasCurrent)
      throw std::logic_error("Current configuration not set.");
    return accelerationEulerAt(t, current);
  }

  // Deformation gradient F for this linear map is exactly K(t).
  // Fij = dxi/dXj

  Matrix<T, n, n> deformationGradient(T t) const { return deformMap(t); }

  T Jacobian(T t) { return det(deformationGradient(t)); }

  // An element at point P at time t has direction and magnitude -> find
  // reference configuaration and vice versa: dx = F.dX
  void elementDeformed() {};

  // Nanson's formula: F^T .dA_x = J dA_X and vice versa
  void areaDeformed() {}

  // Material and spatial derivative of material/spatial description of a
  // property field(scalar)

  // Polar decomposition of the deformation Gradient F = RU = VR; U^2 = F^T.F,
  // V^2 = F.F^T

  // left and right Cauchy deformation tensor

  // Green strain tensor

  // infinite tesimal strain tensor

  // traction vector

  // Cauchy stress tensor

  // Right Cauchy-Green tensor C = F^T F.
  Matrix<T, n, n> rightStretchTensor(T t) const {
    const Matrix<T, n, n> F = deformationGradient(t);
    return F.transpose() * F;
  }

  // Left Cauchy-Green tensor B = F F^T.
  Matrix<T, n, n> leftStretchTensor(T t) const {
    const Matrix<T, n, n> F = deformationGradient(t);
    return F * F.transpose();
  }
};

#endif
