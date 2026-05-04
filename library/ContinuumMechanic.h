#ifndef CONTINUUM_MECHANICS_H
#define CONTINUUM_MECHANICS_H

#include "Coordinates.h"
#include "DualDifferentiation.h"
#include "Matrix.h"
#include "ParentElement.h"
#include "Vector.h"
#include <functional>
#include <stdexcept>

//  Homogeneous transformation.
//  Continuum kinematics with linear motion map: x(X,t) = K(t) * X.
//  Transformation gradient therefore is independent with X ->  K(t) is provided
//  as a matrix-valued function of Dual time.
template <typename T, Index n> struct MediumCon {
  using MapMatrix = std::function<Matrix<Dual, n, n>(Dual)>;
  using MotionMap =
      std::function<StaticVector<Dual, n>(const StaticVector<Dual, n> &, Dual)>;

  // Member variables
  // Reference configuration: X
  StaticVector<T, n> ref{};
  // Current configuration: x
  StaticVector<T, n> current{};
  // Lagrange's description:
  StaticVector<T, n> velL{}, accL{};
  // Euler's description of velocity:
  StaticVector<T, n> velE{}, accE{};

  MapMatrix deformMap{};
  MotionMap motionMap{};

  bool hasRef{false};
  bool hasCurrent{false};

private:
  // The deformation gradient's value and its derivatives:
  // K(t) = phi(t)
  static Matrix<T, n, n> valuePart(const Matrix<Dual, n, n> &mDual) {
    Matrix<T, n, n> out{};
    for (Index i = 0; i < mDual.length(); ++i) {
      out[i] = static_cast<T>(mDual[i].getVal());
    }
    return out;
  }

  // dphi(t)/dt
  static Matrix<T, n, n> derivativePart(const Matrix<Dual, n, n> &mDual) {
    Matrix<T, n, n> out{};
    for (Index i = 0; i < mDual.length(); ++i) {
      out[i] = static_cast<T>(mDual[i].getDer());
    }
    return out;
  }

  // d^2(phi) / (dt)^2
  static Matrix<T, n, n> secondDerivativePart(const Matrix<Dual, n, n> &mDual) {
    Matrix<T, n, n> out{};
    for (Index i = 0; i < mDual.length(); ++i) {
      out[i] = static_cast<T>(mDual[i].getDer2());
    }
    return out;
  }

  static StaticVector<T, n> valuePart(const StaticVector<Dual, n> &vDual) {
    StaticVector<T, n> out{};
    for (Index i = 0; i < n; ++i) {
      out[i] = static_cast<T>(vDual[i].getVal());
    }
    return out;
  }

  StaticVector<T, n> evaluateMotionMap(const StaticVector<T, n> &X, T t) const {
    if (!motionMap)
      throw std::logic_error("motionMap is not set.");

    StaticVector<Dual, n> XDual{};
    for (Index i = 0; i < n; ++i) {
      XDual[i] = Dual{static_cast<double>(X[i]), 0.0, 0.0};
    }
    return valuePart(motionMap(XDual, Dual{static_cast<double>(t), 0.0, 0.0}));
  }

  StaticVector<T, n> velocityMaterialAt(const StaticVector<T, n> &X,
                                        T t) const {
    if (deformMap) {
      return dKdt(t) * X;
    }
    if (!motionMap)
      throw std::logic_error("motionMap is not set.");

    StaticVector<Dual, n> XDual{};
    for (Index i = 0; i < n; ++i) {
      XDual[i] = Dual{static_cast<double>(X[i]), 0.0, 0.0};
    }
    const auto xDual = motionMap(XDual, Dual{static_cast<double>(t), 1.0, 0.0});

    StaticVector<T, n> out{};
    for (Index i = 0; i < n; ++i) {
      out[i] = static_cast<T>(xDual[i].getDer());
    }
    return out;
  }

  StaticVector<T, n> accelerationMaterialAt(const StaticVector<T, n> &X,
                                            T t) const {
    if (deformMap) {
      return d2Kdt2(t) * X;
    }
    if (!motionMap)
      throw std::logic_error("motionMap is not set.");

    StaticVector<Dual, n> XDual{};
    for (Index i = 0; i < n; ++i) {
      XDual[i] = Dual{static_cast<double>(X[i]), 0.0, 0.0};
    }
    const auto xDual = motionMap(XDual, Dual{static_cast<double>(t), 1.0, 0.0});

    StaticVector<T, n> out{};
    for (Index i = 0; i < n; ++i) {
      out[i] = static_cast<T>(xDual[i].getDer2());
    }
    return out;
  }

public:
  MediumCon() = default;
  // Material-description constructor
  MediumCon(const MapMatrix &motionEquation,
            const StaticVector<T, n> &refConfig)
      : ref(refConfig), deformMap(motionEquation), hasRef(true) {}

  // General (possibly non-homogeneous) motion map: x = Phi(X,t)
  MediumCon(const MotionMap &motionEquation,
            const StaticVector<T, n> &refConfig)
      : ref(refConfig), motionMap(motionEquation), hasRef(true) {}

  // Spatial-description constructor
  MediumCon(const MapMatrix &inverseMotionEquation,
            const StaticVector<T, n> &currentConfig, bool isInverseMap)
      : current(currentConfig), hasCurrent(true) {
    if (!isInverseMap)
      throw std::invalid_argument(
          "Set isInverseMap=true when passing inverse motion equation.");

    deformMap = [inverseMotionEquation](Dual t) -> Matrix<Dual, n, n> {
      const Matrix<Dual, n, n> invDual = inverseMotionEquation(t);

      Matrix<T, n, n> invVal{valuePart(invDual)};

      const Matrix<T, n, n> kVal = invVal.inverse();

      Matrix<Dual, n, n> kDual{};
      for (Index i = 0; i < kDual.length(); ++i) {
        kDual[i] = Dual{static_cast<double>(kVal[i]), 0.0};
      }

      return kDual;
    };
  }

  Matrix<T, n, n> motionFunction(T t) const {
    if (deformMap)
      return valuePart(deformMap(Dual{static_cast<double>(t), 0.0}));
    if (motionMap && hasRef)
      return deformationGradient(t, ref);
    throw std::logic_error("motionFunction is unavailable: set deformMap or "
                           "set both motionMap and ref.");
  }

  Matrix<T, n, n> inverseMap(T t) const {
    if (!deformMap)
      throw std::logic_error(
          "inverseMap is only available for homogeneous linear map K(t). ");
    return valuePart(deformMap(Dual{static_cast<double>(t), 0.0})).inverse();
  }

  // dK/dt from Dual derivative part.
  Matrix<T, n, n> dKdt(T t) const {
    if (!deformMap)
      throw std::logic_error("dKdt is only available for homogeneous K(t). ");
    return derivativePart(deformMap(Dual{static_cast<double>(t), 1.0}));
  }

  // d2K/dt2 from Dual second derivative part.
  Matrix<T, n, n> d2Kdt2(T t) const {
    if (!deformMap)
      throw std::logic_error("d2Kdt2 is only available for homogeneous K(t). ");
    return secondDerivativePart(
        deformMap(Dual{static_cast<double>(t), 1.0, 0.0}));
  }

  // x_P(t) = K(t) * X
  StaticVector<T, n> currentPosition(T t) const {
    if (!hasRef)
      throw std::logic_error("Reference configuration not set.");
    if (deformMap)
      return motionFunction(t) * ref;
    return evaluateMotionMap(ref, t);
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
    return velocityMaterialAt(ref, t);
  }

  // Material displacement: u(X,t) = x(X,t) - X.
  StaticVector<T, n> displacementMaterial(T t) const {
    if (!hasRef)
      throw std::logic_error("Reference configuration not set.");
    return currentPosition(t) - ref;
  }

  // Spatial displacement at stored current position: u(x,t) = x - X(x,t).
  StaticVector<T, n> displacementSpatial(T t) const {
    if (!hasCurrent)
      throw std::logic_error("Current configuration not set.");
    return current - referencePosition(t);
  }

  // Material acceleration: a(X,t) = d2x/dt2 = d2K/dt2 * X.
  StaticVector<T, n> accelerationMaterial(T t) const {
    if (!hasRef)
      throw std::logic_error("Reference configuration not set.");
    return accelerationMaterialAt(ref, t);
  }

  // Euler velocity at any spatial point x:
  // v(x,t) = (dK/dt * K^{-1}) * x.
  StaticVector<T, n> velEulerAt(T t, const StaticVector<T, n> &xSpatial) const {
    const Matrix<T, n, n> L = dKdt(t) * inverseMap(t);
    return L * xSpatial;
  }

  // Euler velocity at stored current position.
  StaticVector<T, n> velEuler(T t) const {
    if (!hasCurrent)
      throw std::logic_error("Current configuration not set.");
    return velEulerAt(t, current);
  }

  // Euler acceleration at any spatial point x:
  // a = (dL/dt + L^2) x, where L = dK/dt * K^{-1}.
  StaticVector<T, n> accEulerAt(T t, const StaticVector<T, n> &xSpatial) const {
    const Matrix<T, n, n> L = dKdt(t) * inverseMap(t);
    const Matrix<T, n, n> dLdt = d2Kdt2(t) * inverseMap(t) - L * L;

    return (dLdt + L * L) * xSpatial;
  }

  // Euler acceleration at stored current position.
  StaticVector<T, n> accEuler(T t) const {
    if (!hasCurrent)
      throw std::logic_error("Current configuration not set.");
    return accEulerAt(t, current);
  }

  // Deformation gradient F for this linear map is exactly K(t).
  // Fij = dxi/dXj

  Matrix<T, n, n> deformationGradient(T t) const {
    if (deformMap)
      return motionFunction(t);
    if (!hasRef)
      throw std::logic_error(
          "deformationGradient(t) for non-homogeneous map requires ref.");
    return deformationGradient(t, ref);
  }

  // General deformation gradient at arbitrary material point X:
  // F(X,t) = dPhi/dX.
  Matrix<T, n, n> deformationGradient(T t, const StaticVector<T, n> &X) const {
    if (deformMap)
      return motionFunction(t);
    if (!motionMap)
      throw std::logic_error("motionMap is not set.");

    Matrix<T, n, n> F{};
    for (Index j = 0; j < n; ++j) {
      StaticVector<Dual, n> XDual{};
      for (Index k = 0; k < n; ++k) {
        const double seed = (k == j) ? 1.0 : 0.0;
        XDual[k] = Dual{static_cast<double>(X[k]), seed, 0.0};
      }
      const auto xDual =
          motionMap(XDual, Dual{static_cast<double>(t), 0.0, 0.0});
      for (Index i = 0; i < n; ++i) {
        F(i, j) = static_cast<T>(xDual[i].getDer());
      }
    }
    return F;
  }

  T Jacobian(T t) const { return det(deformationGradient(t)); }

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
