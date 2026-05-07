#ifndef CONTINUUM_MECHANICS_H
#define CONTINUUM_MECHANICS_H

#include "Coordinates.h"
#include "DualDifferentiation.h"
#include "Matrix.h"
#include "ParentElement.h"
#include "Vector.h"
#include <functional>
#include <stdexcept>

template <typename T, Index n> struct MediumCon {
  // linear-Homogeneous transformation :  K(t) such that x(X, t) = K(t) * X.
  // Time t takes type Dual as K(t) will need to be derivated by time
  using MapMatrix = std::function<Matrix<Dual, n, n>(Dual)>;
  // General motion equation - non linear: x(X,t) = Phi(X,t)
  using MotionEquation =
      std::function<StaticVector<Dual, n>(const StaticVector<Dual, n> &, Dual)>;

  // Member variables
  // Reference configuration: X
  StaticVector<T, n> ref{};
  // Current configuration: x
  StaticVector<T, n> cur{};
  // Lagrange's description:
  StaticVector<T, n> velL{}, accL{};
  // Euler's description of velocity:
  StaticVector<T, n> velE{}, accE{};

  MapMatrix linearMap{};      // homogeneous part
  MotionEquation deformMap{}; // general part

  bool hasRef{false};
  bool hasCur{false};

private:
  // Note: helper extraction functions for Dual-like scalars are provided as
  // member functions on `Matrix` and `StaticVector` (see Matrix.h / Vector.h).

  StaticVector<T, n> evaluateMotionMap(const StaticVector<T, n> &X, T t) const {
    if (!deformMap)
      throw std::logic_error("deformMap is not set.");

    StaticVector<Dual, n> XDual{};
    for (Index i = 0; i < n; ++i) {
      XDual[i] = Dual{static_cast<double>(X[i]), 0.0, 0.0};
    }
    const auto resDual =
        deformMap(XDual, Dual{static_cast<double>(t), 0.0, 0.0});
    return resDual.template val<T>();
  }

  StaticVector<T, n> velocityMaterialAt(const StaticVector<T, n> &X,
                                        T t) const {
    if (linearMap) {
      return dKdt(t) * X;
    }
    if (!deformMap)
      throw std::logic_error("deformMap is not set.");

    StaticVector<Dual, n> XDual{};
    for (Index i = 0; i < n; ++i) {
      XDual[i] = Dual{static_cast<double>(X[i]), 0.0, 0.0};
    }
    const auto xDual = deformMap(XDual, Dual{static_cast<double>(t), 1.0, 0.0});
    return xDual.template der1<T>();
  }

  StaticVector<T, n> accelerationMaterialAt(const StaticVector<T, n> &X,
                                            T t) const {
    if (linearMap) {
      return d2Kdt2(t) * X;
    }
    if (!deformMap)
      throw std::logic_error("deformMap is not set.");

    StaticVector<Dual, n> XDual{};
    for (Index i = 0; i < n; ++i) {
      XDual[i] = Dual{static_cast<double>(X[i]), 0.0, 0.0};
    }
    const auto xDual = deformMap(XDual, Dual{static_cast<double>(t), 1.0, 0.0});
    return xDual.template der2<T>();
  }

public:
  MediumCon() = default;
  // Material-description constructor
  MediumCon(const MapMatrix &motionEquation,
            const StaticVector<T, n> &refConfig)
      : ref(refConfig), linearMap(motionEquation), hasRef(true) {}

  // General (possibly non-homogeneous) motion map: x = Phi(X,t)
  MediumCon(const MotionEquation &motionEquation,
            const StaticVector<T, n> &refConfig)
      : ref(refConfig), deformMap(motionEquation), hasRef(true) {}

  // Spatial-description constructor
  MediumCon(const MapMatrix &inverseMotionEquation,
            const StaticVector<T, n> &currentConfig, bool isInverseMap)
      : cur(currentConfig), hasCur(true) {
    if (!isInverseMap)
      throw std::invalid_argument(
          "Set isInverseMap=true when passing inverse motion equation.");

    linearMap = [inverseMotionEquation](Dual t) -> Matrix<Dual, n, n> {
      const Matrix<Dual, n, n> invDual = inverseMotionEquation(t);

      Matrix<T, n, n> invVal = invDual.template val<T>();

      const Matrix<T, n, n> kVal = invVal.inverse();

      Matrix<Dual, n, n> kDual{};
      for (Index i = 0; i < kDual.length(); ++i) {
        kDual[i] = Dual{static_cast<double>(kVal[i]), 0.0};
      }

      return kDual;
    };
  }

  Matrix<T, n, n> motionFunction(T t) const {
    if (linearMap)
      return linearMap(Dual{static_cast<double>(t), 0.0}).template val<T>();
    if (deformMap && hasRef)
      return deformationGradient(t, ref);
    throw std::logic_error("motionFunction is unavailable: set linearMap or "
                           "set both deformMap and ref.");
  }

  Matrix<T, n, n> inverseMap(T t) const {
    if (!linearMap)
      throw std::logic_error(
          "inverseMap is only available for homogeneous linear map K(t). ");
    return linearMap(Dual{static_cast<double>(t), 0.0})
        .template val<T>()
        .inverse();
  }

  // dK/dt from Dual derivative part.
  Matrix<T, n, n> dKdt(T t) const {
    if (!linearMap)
      throw std::logic_error("dKdt is only available for homogeneous K(t). ");
    return linearMap(Dual{static_cast<double>(t), 1.0}).template der1<T>();
  }

  // d2K/dt2 from Dual second derivative part.
  Matrix<T, n, n> d2Kdt2(T t) const {
    if (!linearMap)
      throw std::logic_error("d2Kdt2 is only available for homogeneous K(t). ");
    return linearMap(Dual{static_cast<double>(t), 1.0, 0.0}).template der2<T>();
  }

  // x_P(t) = K(t) * X
  StaticVector<T, n> currentPosition(T t) const {
    if (!hasRef)
      throw std::logic_error("Reference configuration not set.");
    if (linearMap)
      return motionFunction(t) * ref;
    return evaluateMotionMap(ref, t);
  }

  // X = K^{-1}(t) * x
  StaticVector<T, n> referencePosition(T t) const {
    if (!hasCur)
      throw std::logic_error("Current configuration not set.");
    return inverseMap(t) * cur;
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

  // Spatial displacement at stored cur position: u(x,t) = x - X(x,t).
  StaticVector<T, n> displacementSpatial(T t) const {
    if (!hasCur)
      throw std::logic_error("Current configuration not set.");
    return cur - referencePosition(t);
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

  // Euler velocity at stored cur position.
  StaticVector<T, n> velEuler(T t) const {
    if (!hasCur)
      throw std::logic_error("Current configuration not set.");
    return velEulerAt(t, cur);
  }

  // Euler acceleration at any spatial point x:
  // a = (dL/dt + L^2) x, where L = dK/dt * K^{-1}.
  StaticVector<T, n> accEulerAt(T t, const StaticVector<T, n> &xSpatial) const {
    const Matrix<T, n, n> L = dKdt(t) * inverseMap(t);
    const Matrix<T, n, n> dLdt = d2Kdt2(t) * inverseMap(t) - L * L;

    return (dLdt + L * L) * xSpatial;
  }

  // Euler acceleration at stored cur position.
  StaticVector<T, n> accEuler(T t) const {
    if (!hasCur)
      throw std::logic_error("Current configuration not set.");
    return accEulerAt(t, cur);
  }

  // Deformation gradient F for this linear map is exactly K(t).
  // Fij = dxi/dXj

  Matrix<T, n, n> deformationGradient(T t) const {
    if (linearMap)
      return motionFunction(t);
    if (!hasRef)
      throw std::logic_error(
          "deformationGradient(t) for non-homogeneous map requires ref.");
    return deformationGradient(t, ref);
  }

  // General deformation gradient at arbitrary material point X:
  // F(X,t) = dPhi/dX.
  Matrix<T, n, n> deformationGradient(T t, const StaticVector<T, n> &X) const {
    if (linearMap)
      return motionFunction(t);
    if (!deformMap)
      throw std::logic_error("deformMap is not set.");

    Matrix<T, n, n> F{};
    for (Index j = 0; j < n; ++j) {
      StaticVector<Dual, n> XDual{};
      for (Index k = 0; k < n; ++k) {
        const double seed = (k == j) ? 1.0 : 0.0;
        XDual[k] = Dual{static_cast<double>(X[k]), seed, 0.0};
      }
      const auto xDual =
          deformMap(XDual, Dual{static_cast<double>(t), 0.0, 0.0});
      const auto deriv = xDual.template der1<T>();
      for (Index i = 0; i < n; ++i) {
        F(i, j) = deriv[i];
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
