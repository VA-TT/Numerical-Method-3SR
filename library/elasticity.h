#ifndef HOOKE_ELASTICITY_H
#define HOOKE_ELASTICITY_H

#include "Matrix.h"
#include "physicConstants.h"
#include <cmath>
#include <iostream>
#include <string_view>
#include <tuple>
#include <utility>

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

template <typename T>
Matrix<T, 3, 3> elasticityMatrix(T E, T nu,
                                 std::string_view type = "planeStrain") {
  return tangentOperator2D(type, E, nu);
}

template <typename T, Index nRows, Index nCols>
StaticVector<T, nRows>
principalStress(const Matrix<T, nRows, nCols> &stress_tensor) {
  assert(stress_tensor.isSymmetric() && "Stress tensor should be symmetric!");
  auto [principalSigmas, V] = stress_tensor.eigen();
  return principalSigmas.sorted(true); // sigma1 >= sigma2 >= sigma3
}

// Drucker-Prager: compute alpha and k from c and phi
template <typename T> std::pair<T, T> druckerPrager(T phi_degree, T c) {
  T phi = constants::pi * phi_degree / T{180.0};
  T alpha =
      T{2.0} * std::sin(phi) / (std::sqrt(T{3.0}) * (T{3.0} - std::sin(phi)));
  T k = (T{6.0} * c * std::cos(phi)) /
        (std::sqrt(T{3.0}) * (T{3.0} - std::sin(phi)));
  return {alpha, k};
}

template <typename T, Index nRows, Index nCols>
std::pair<T, T> invariants(const Matrix<T, nRows, nCols> &stress_tensor,
                           std::string_view type = "planeStrain",
                           T nu = T{0.33}) {
  static_assert(nRows >= 2 && nCols >= 2,
                "invariants: stress tensor must be at least 2x2");

  const T sigma_xx = stress_tensor(0, 0);
  const T tau_xy = stress_tensor(0, 1);
  const T sigma_yy = stress_tensor(1, 1);

  if (type != "planeStrain") {
    const T I1 = sigma_xx + sigma_yy;

    const T s_xx = sigma_xx - I1 / T{2};
    const T s_yy = sigma_yy - I1 / T{2};
    const T s_xy = tau_xy;

    const T J2 = T{0.5} * (s_xx * s_xx + s_yy * s_yy + s_xy * s_xy);
    return {I1, J2};
  }

  const T sigma_zz = nu * (sigma_xx + sigma_yy);
  const T I1 = sigma_xx + sigma_yy + sigma_zz;

  const T s_xx = sigma_xx - I1 / T{3};
  const T s_yy = sigma_yy - I1 / T{3};
  const T s_zz = sigma_zz - I1 / T{3};
  const T s_xy = tau_xy;

  const T J2 = T{0.5} * (s_xx * s_xx + s_yy * s_yy + s_zz * s_zz + s_xy * s_xy);
  return {I1, J2};
}

template <typename T>
std::pair<T, T> invariants(const DynamicVector<T> &stress,
                           std::string_view type = "planeStrain",
                           T nu = T{0.33}) {
  assert(stress.size() >= 3 && "invariants: stress must be [xx, yy, xy]");
  const T sigma_xx = stress[0];
  const T sigma_yy = stress[1];
  const T tau_xy = stress[2];

  if (type != "planeStrain") {
    const T I1 = sigma_xx + sigma_yy;

    const T s_xx = sigma_xx - I1 / T{2};
    const T s_yy = sigma_yy - I1 / T{2};
    const T s_xy = tau_xy;

    const T J2 = T{0.5} * (s_xx * s_xx + s_yy * s_yy + s_xy * s_xy);
    return {I1, J2};
  }

  const T sigma_zz = nu * (sigma_xx + sigma_yy);
  const T I1 = sigma_xx + sigma_yy + sigma_zz;

  const T s_xx = sigma_xx - I1 / T{3};
  const T s_yy = sigma_yy - I1 / T{3};
  const T s_zz = sigma_zz - I1 / T{3};
  const T s_xy = tau_xy;

  const T J2 = T{0.5} * (s_xx * s_xx + s_yy * s_yy + s_zz * s_zz + s_xy * s_xy);
  return {I1, J2};
}

template <typename T>
T druckerPragerYield(const DynamicVector<T> &stress, T alpha, T k,
                     std::string_view stressCondition = "planeStrain",
                     T nu = T{0.33}) {
  const auto [I1, J2] = invariants(stress, stressCondition, nu);
  return alpha * I1 + std::sqrt(J2) - k;
}

template <typename T>
DynamicVector<T>
druckerPragerGradient(const DynamicVector<T> &stress, T alpha,
                      std::string_view stressCondition = "planeStrain",
                      T nu = T{0.33}) {
  assert(stress.size() >= 3 &&
         "druckerPragerGradient: stress must be [xx, yy, xy]");
  const T sigma_x = stress[0];
  const T sigma_y = stress[1];
  const T tau_xy = stress[2];

  const auto [I1, J2] = invariants(stress, stressCondition, nu);
  (void)I1;

  const T denom = std::sqrt(J2);
  if (approximatelyEqualAbsRel(static_cast<double>(denom), 0.0)) {
    return DynamicVector<T>{alpha, alpha, T{0}};
  }

  const T dfdx = alpha + (sigma_x - sigma_y) / denom;
  const T dfdy = alpha - (sigma_x - sigma_y) / denom;
  const T dftau = (T{2} * tau_xy) / denom;
  return DynamicVector<T>{dfdx, dfdy, dftau};
}

template <typename T>
std::tuple<DynamicVector<T>, DynamicVector<T>, T>
updateStressStrainDruckerPrager(
    const DynamicVector<T> &stress_n, const DynamicVector<T> &strain_n,
    const DynamicVector<T> &strain_increment, const Matrix<T, 3, 3> &D, T alpha,
    T k, std::string_view stressCondition = "planeStrain", T nu = T{0.33}) {
  assert(stress_n.size() >= 3 && strain_n.size() >= 3 &&
         strain_increment.size() >= 3 &&
         "updateStressStrainDruckerPrager: stress/strain must be size 3");

  auto applyD = [&](const DynamicVector<T> &v) {
    Matrix<T, 3, 1> v_m{DynamicVector<T>{v[0], v[1], v[2]}};
    Matrix<T, 3, 1> out_m = D * v_m;
    return DynamicVector<T>(out_m);
  };

  const DynamicVector<T> stress_trial =
      DynamicVector<T>{stress_n[0], stress_n[1], stress_n[2]} +
      applyD(strain_increment);

  const T f_trial =
      druckerPragerYield(stress_trial, alpha, k, stressCondition, nu);
  if (f_trial <= T{0}) {
    const DynamicVector<T> strain_updated =
        DynamicVector<T>{strain_n[0], strain_n[1], strain_n[2]} +
        DynamicVector<T>{strain_increment[0], strain_increment[1],
                         strain_increment[2]};
    return {stress_trial, strain_updated, T{0}};
  }

  const DynamicVector<T> df_sigma =
      druckerPragerGradient(stress_trial, alpha, stressCondition, nu);
  const DynamicVector<T> D_df = applyD(df_sigma);
  const T H = dotProduct(df_sigma, D_df);

  if (approximatelyEqualAbsRel(static_cast<double>(H), 0.0)) {
    return {stress_trial,
            DynamicVector<T>{strain_n[0], strain_n[1], strain_n[2]}, T{0}};
  }

  const T delta_lambda = f_trial / H;
  const DynamicVector<T> stress_updated = stress_trial - delta_lambda * D_df;

  const DynamicVector<T> plastic_strain_increment = delta_lambda * df_sigma;
  const DynamicVector<T> strain_updated =
      DynamicVector<T>{strain_n[0], strain_n[1], strain_n[2]} +
      plastic_strain_increment;
  return {stress_updated, strain_updated, delta_lambda};
}
#endif // HOOKE_ELASTIC_H
