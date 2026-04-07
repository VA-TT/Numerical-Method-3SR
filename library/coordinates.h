#ifndef COORDINATE_SYSTEMS_H
#define COORDINATE_SYSTEMS_H

#include "DualDiffrentiation.h"
#include "Matrix.h"
#include "Vector.h"
#include "comparison.h"
#include "physicConstants.h"
#include <cmath>
#include <stdexcept>
#include <string>
#include <string_view>

// linear transformation of a coordinate system != Cartesian
template <typename T, Index n> struct CoordSystem {
  CoordSystem(const StaticVector<StaticVector<T, n>, n> &b_new) : b{b_new} {
    for (Index i{0}; i < n; ++i) {
      for (Index j{i}; j < n; ++j) {
        G(i, j) = dotProduct(b[i], b[j]);
      }
    }
    G.reflect(); // G is symmetric
    for (Index j{0}; j < n; ++j) {
      for (Index i{0}; i < n; ++i) {
        A(i, j) = b_new[j][i];
      }
    }
  }
  CoordSystem(const Matrix<T, n, n> &transformMatrix) : A{transformMatrix} {
    for (Index j{0}; j < n; ++j) {
      b[j] = transformMatrix.getColVector(j);
    }
    for (Index i{0}; i < n; ++i) {
      for (Index j{i}; j < n; ++j) {
        G(i, j) = dotProduct(b[i], b[j]);
      }
    }
    G.reflect(); // G is symmetric
  }

  StaticVector<StaticVector<T, n>, n> b; // new basis
  Matrix<T, n, n> G{};                   // Matrix metric

  Matrix<T, n, n> A{}; // Matrix of the transformation: T(x) = Ax
                       // (standard reprensentation matrix)

  T norm(const StaticVector<T, n> &v) const {
    return std::sqrt(dotProduct(v, G * v));
  }

  static Matrix<T, n, n> changeBasisMatrix(const CoordSystem &B,
                                           const CoordSystem &C) {
    Matrix<T, n, n> P = C.A.inverse() * B.A;
    return P;
  }

  StaticVector<T, n> newBasisCoord(const StaticVector<T, n> &v_B,
                                   const CoordSystem &C) const {
    Matrix<T, n, n> P = changeBasisMatrix(*this, C); // Pc<-b
    return P * v_B;
  }

  StaticVector<T, n> oldBasisCoord(const StaticVector<T, n> &v_C,
                                   const CoordSystem &C) const {
    Matrix<T, n, n> P = changeBasisMatrix(*this, C); // Pc->b
    return P.inverse() * v_C;
  }

  // Change linear-operator matrix from basis B (this) to basis C:
  // [L]_C = P_{C<-B} [L]_B P_{B<-C}
  Matrix<T, n, n> newBasisOperator(const Matrix<T, n, n> &L_B,
                                   const CoordSystem &C) const {
    const Matrix<T, n, n> P = changeBasisMatrix(*this, C);
    return P * L_B * P.inverse();
  }

  // Change linear-operator matrix from basis C back to basis B (this):
  // [L]_B = P_{B<-C} [L]_C P_{C<-B}
  Matrix<T, n, n> oldBasisOperator(const Matrix<T, n, n> &L_C,
                                   const CoordSystem &C) const {
    const Matrix<T, n, n> P = changeBasisMatrix(*this, C);
    return P.inverse() * L_C * P;
  }

  // (R o S) (u) = R(S(u)) = B.A.u
  // We can change from Cartesian -> another -> another
  friend CoordSystem operator*(const CoordSystem &R, const CoordSystem &S) {
    return CoordSystem(R.A * S.A);
  }
};

// 2D polar coordinate
template <typename T> struct PolarCoord {
  // Constructor using (r,t) variables
  PolarCoord(T r, T t) // theta is in degree unit
  {
    this->r = r;
    this->t = t;
    t_rad = t * T{constants::pi} / T{180};
    cost = std::cos(t_rad);
    sint = std::sin(t_rad);
    x = r * cost;
    y = r * sint;
    rUnit = StaticVector<T, 2>{cost, sint};
    tUnit = StaticVector<T, 2>{-sint, cost};
    drUnitdt = tUnit;
    dtUnitdt = -rUnit;
  }

  // Constructor using (x,y, "Cartesian") variables
  PolarCoord(T x, T y, std::string_view coordSystem)
      : PolarCoord(std::sqrt(x * x + y * y),
                   std::atan2(y, x) * T{180} / T{constants::pi}) {
    if (coordSystem != "Cartesian")
      throw std::invalid_argument("PolarCoord: coordSystem must be Cartesian.");

    // Preserve exact Cartesian inputs from caller.
    this->x = x;
    this->y = y;
  }

  // Member variables
  T x{}, y{};
  T r{}, t{};
  T t_rad{};
  T cost{}, sint{};
  StaticVector<T, 2> rUnit, tUnit{};
  StaticVector<T, 2> drUnitdt, dtUnitdt{};
  ScalarFunction2D f; // f = f(r,t)
  VectorFunction2D u; // u = (u1(r,t), u2(r,t))
  StaticVector<T, 2> gradient(ScalarFunction2D f_rt) const {
    if (approximatelyEqualAbsRel(r, T{0}))
      throw std::invalid_argument("PolarCoord::gradient requires r != 0.");

    Dual r_dr{r, 1.0};
    Dual t_dr{t, 0.0};
    auto f_dualR = f_rt(r_dr, t_dr);
    double dfdr = f_dualR.getDer();
    Dual r_dt{r, 0.0};
    Dual t_dt{t, 1.0};
    auto f_dualT = f_rt(r_dt, t_dt);
    double dfdt = f_dualT.getDer();
    return {dfdr, T{1} / r * dfdt}; // polor coordinate
  }

  // Member functions (differential operator)
  T jacobian() const { return r; } // dxdy = |J| drdt

  T div(VectorFunction2D u_rt) const {
    if (approximatelyEqualAbsRel(r, T{0}))
      throw std::invalid_argument("PolarCoord::div requires r != 0.");

    Dual r_dr{r, 1.0};
    Dual t_dr{t, 0.0};
    auto u_dualR = u_rt(r_dr, t_dr);
    double du1dr = u_dualR[0].getDer();
    double u1 = u_dualR[0].getVal();

    Dual r_dt{r, 0.0};
    Dual t_dt{t, 1.0};
    auto u_dualT = u_rt(r_dt, t_dt);
    double du2dt = u_dualT[1].getDer();

    return du1dr + T{1} / r * (u1 + du2dt);
  }

  // 2D curl = scalar z-component of the 3D curl:
  T curl(VectorFunction2D u_rt) const {
    if (approximatelyEqualAbsRel(r, T{0}))
      throw std::invalid_argument("PolarCoord::curl requires r != 0.");

    Dual r_dr{r, 1.0};
    Dual t_dr{t, 0.0};
    auto u_dualR = u_rt(r_dr, t_dr);
    double du2dr = u_dualR[1].getDer();

    Dual r_dt{r, 0.0};
    Dual t_dt{t, 1.0};
    auto u_dualT = u_rt(r_dt, t_dt);
    double du1dt = u_dualT[0].getDer();
    double u2 = u_dualT[1].getVal();

    return du2dr + T{1} / r * (u2 - du1dt); // k component
  }

  // Scalar Laplacian in polar coordinates:
  // nabla^2 f = d2f/dr2 + (1/r) df/dr + (1/r^2) d2f/dt2
  T laplacian(ScalarFunction2D f_rt, T hr = T{1e-5}, T ht = T{1e-5}) const {
    if (approximatelyEqualAbsRel(r, T{0}))
      throw std::invalid_argument("PolarCoord::laplacian requires r != 0.");
    if (hr <= T{0} || ht <= T{0})
      throw std::invalid_argument(
          "PolarCoord::laplacian requires positive steps.");

    auto eval = [&](T rr, T tt) -> T {
      return static_cast<T>(f_rt(Dual{static_cast<double>(rr), 0.0},
                                 Dual{static_cast<double>(tt), 0.0})
                                .getVal());
    };

    const T f0 = eval(r, t);
    const T fr_p = eval(r + hr, t);
    const T fr_m = eval(r - hr, t);
    const T ft_p = eval(r, t + ht);
    const T ft_m = eval(r, t - ht);

    const T d2fdr2 = (fr_p - T{2} * f0 + fr_m) / (hr * hr);
    const T dfdr = (fr_p - fr_m) / (T{2} * hr);
    const T d2fdt2 = (ft_p - T{2} * f0 + ft_m) / (ht * ht);

    return d2fdr2 + T{1} / r * dfdr + T{1} / (r * r) * d2fdt2;
  }

  // Vector Laplacian in polar coordinates (component-wise).
  StaticVector<T, 2> laplacian(VectorFunction2D u_rt, T hr = T{1e-5},
                               T ht = T{1e-5}) const {
    auto u1 = [u_rt](Dual rr, Dual tt) -> Dual { return u_rt(rr, tt)[0]; };
    auto u2 = [u_rt](Dual rr, Dual tt) -> Dual { return u_rt(rr, tt)[1]; };
    return {laplacian(u1, hr, ht), laplacian(u2, hr, ht)};
  }

  // Function f = f(x,y) df/dr and df/dt value
  T dfdr(ScalarFunction2D f_xy) {
    auto [dfdx, dfdy] = grad2D(f_xy, x, y);
    return cost * dfdx + sint * dfdy;
  }
  T dfdt(ScalarFunction2D f_xy) {
    auto [dfdx, dfdy] = grad2D(f_xy, x, y);
    return -r * sint * dfdx + r * cost * dfdy;
  }
};

// 3D Cylindrical coordinates
template <typename T> struct CylinCoord {
  // Constructor using cylindrical variables (r, p, z), p in degrees.
  CylinCoord(T r, T p, T z) {
    this->r = r;
    this->p = p;
    this->z = z;

    p_rad = p * T{constants::pi} / T{180};
    cosp = std::cos(p_rad);
    sinp = std::sin(p_rad);

    x = r * cosp;
    y = r * sinp;

    rUnit = StaticVector<T, 3>{cosp, sinp, T{0}};
    pUnit = StaticVector<T, 3>{-sinp, cosp, T{0}};
    zUnit = StaticVector<T, 3>{T{0}, T{0}, T{1}};
    position = r * rUnit + z * zUnit;
  }

  // Constructor using Cartesian variables (x, y, z).
  CylinCoord(T x, T y, T z, std::string_view coordSystem) {
    if (coordSystem != "Cartesian")
      throw std::invalid_argument("CylinCoord: coordSystem must be Cartesian.");

    this->x = x;
    this->y = y;
    this->z = z;

    r = std::sqrt(x * x + y * y);
    p = std::atan2(y, x) * T{180} / T{constants::pi};

    p_rad = p * T{constants::pi} / T{180};
    cosp = std::cos(p_rad);
    sinp = std::sin(p_rad);

    rUnit = StaticVector<T, 3>{cosp, sinp, T{0}};
    pUnit = StaticVector<T, 3>{-sinp, cosp, T{0}};
    zUnit = StaticVector<T, 3>{T{0}, T{0}, T{1}};
    position = r * rUnit + z * zUnit;
  }

  // Member variables
  T x{}, y{}, z{};
  T r{}, p{};
  T p_rad{};
  T cosp{}, sinp{};
  StaticVector<T, 3> position{};
  StaticVector<T, 3> rUnit, pUnit{}, zUnit{};
  ScalarFunction3D f; // f = f(r,p,z)
  VectorFunction3D u; // u = (u1(r,p,z), u2(r,p,z), u3(r,p,z))

  // Differential operators
  T jacobian() const { return r; } // dxdydz = |J| drdpdz

  StaticVector<T, 3> gradient(ScalarFunction3D f_rpz) const {
    if (approximatelyEqualAbsRel(r, T{0}))
      throw std::invalid_argument("CylinCoord::gradient requires r != 0.");

    const T dfdr = static_cast<T>(
        f_rpz(Dual{r, 1.0}, Dual{p, 0.0}, Dual{z, 0.0}).getDer());
    const T dfdp = static_cast<T>(
        f_rpz(Dual{r, 0.0}, Dual{p, 1.0}, Dual{z, 0.0}).getDer());
    const T dfdz = static_cast<T>(
        f_rpz(Dual{r, 0.0}, Dual{p, 0.0}, Dual{z, 1.0}).getDer());

    return {dfdr, T{1} / r * dfdp, dfdz};
  }

  T div(VectorFunction3D u_rpz) const {
    if (approximatelyEqualAbsRel(r, T{0}))
      throw std::invalid_argument("CylinCoord::div requires r != 0.");

    const auto u_dr = u_rpz(Dual{r, 1.0}, Dual{p, 0.0}, Dual{z, 0.0});
    const T du1dr = static_cast<T>(u_dr[0].getDer());
    const T u1 = static_cast<T>(u_dr[0].getVal());

    const auto u_dp = u_rpz(Dual{r, 0.0}, Dual{p, 1.0}, Dual{z, 0.0});
    const T du2dp = static_cast<T>(u_dp[1].getDer());

    const auto u_dz = u_rpz(Dual{r, 0.0}, Dual{p, 0.0}, Dual{z, 1.0});
    const T du3dz = static_cast<T>(u_dz[2].getDer());

    return du1dr + T{1} / r * (u1 + du2dp) + du3dz;
  }

  // Curl in cylindrical components (r, p, z)
  StaticVector<T, 3> curl(VectorFunction3D u_rpz) const {
    if (approximatelyEqualAbsRel(r, T{0}))
      throw std::invalid_argument("CylinCoord::curl requires r != 0.");

    const auto u_dr = u_rpz(Dual{r, 1.0}, Dual{p, 0.0}, Dual{z, 0.0});
    const T du2dr = static_cast<T>(u_dr[1].getDer());
    const T du3dr = static_cast<T>(u_dr[2].getDer());

    const auto u_dp = u_rpz(Dual{r, 0.0}, Dual{p, 1.0}, Dual{z, 0.0});
    const T du1dp = static_cast<T>(u_dp[0].getDer());
    const T du3dp = static_cast<T>(u_dp[2].getDer());
    const T u2 = static_cast<T>(u_dp[1].getVal());

    const auto u_dz = u_rpz(Dual{r, 0.0}, Dual{p, 0.0}, Dual{z, 1.0});
    const T du1dz = static_cast<T>(u_dz[0].getDer());
    const T du2dz = static_cast<T>(u_dz[1].getDer());

    const T curl_r = T{1} / r * du3dp - du2dz;
    const T curl_p = du1dz - du3dr;
    const T curl_z = du2dr + T{1} / r * (u2 - du1dp);

    return {curl_r, curl_p, curl_z};
  }

  // Scalar Laplacian:
  // nabla^2 f = d2f/dr2 + (1/r)df/dr + (1/r^2)d2f/dp2 + d2f/dz2
  T laplacian(ScalarFunction3D f_rpz, T hr = T{1e-5}, T hp = T{1e-5},
              T hz = T{1e-5}) const {
    if (approximatelyEqualAbsRel(r, T{0}))
      throw std::invalid_argument("CylinCoord::laplacian requires r != 0.");
    if (hr <= T{0} || hp <= T{0} || hz <= T{0})
      throw std::invalid_argument(
          "CylinCoord::laplacian requires positive steps.");

    auto eval = [&](T rr, T pp, T zz) -> T {
      return static_cast<T>(f_rpz(Dual{static_cast<double>(rr), 0.0},
                                  Dual{static_cast<double>(pp), 0.0},
                                  Dual{static_cast<double>(zz), 0.0})
                                .getVal());
    };

    const T f0 = eval(r, p, z);
    const T fr_p = eval(r + hr, p, z);
    const T fr_m = eval(r - hr, p, z);
    const T fp_p = eval(r, p + hp, z);
    const T fp_m = eval(r, p - hp, z);
    const T fz_p = eval(r, p, z + hz);
    const T fz_m = eval(r, p, z - hz);

    const T d2fdr2 = (fr_p - T{2} * f0 + fr_m) / (hr * hr);
    const T dfdr = (fr_p - fr_m) / (T{2} * hr);
    const T d2fdp2 = (fp_p - T{2} * f0 + fp_m) / (hp * hp);
    const T d2fdz2 = (fz_p - T{2} * f0 + fz_m) / (hz * hz);

    return d2fdr2 + T{1} / r * dfdr + T{1} / (r * r) * d2fdp2 + d2fdz2;
  }

  // Vector Laplacian (component-wise)
  StaticVector<T, 3> laplacian(VectorFunction3D u_rpz, T hr = T{1e-5},
                               T hp = T{1e-5}, T hz = T{1e-5}) const {
    auto u1 = [u_rpz](Dual rr, Dual pp, Dual zz) -> Dual {
      return u_rpz(rr, pp, zz)[0];
    };
    auto u2 = [u_rpz](Dual rr, Dual pp, Dual zz) -> Dual {
      return u_rpz(rr, pp, zz)[1];
    };
    auto u3 = [u_rpz](Dual rr, Dual pp, Dual zz) -> Dual {
      return u_rpz(rr, pp, zz)[2];
    };

    return {laplacian(u1, hr, hp, hz), laplacian(u2, hr, hp, hz),
            laplacian(u3, hr, hp, hz)};
  }
};

// 3D Spherical coordinates
template <typename T> struct SphereCoord {
public:
  // Constructor using (r,t,p) where:
  // t in [0,180] deg is polar angle (theta),
  // p in [0,360] deg is azimuth angle (phi).
  SphereCoord(T r, T t, T p) {
    if (t < T{0} || t > T{180})
      throw std::invalid_argument(
          "SphereCoord: theta must be in [0, 180] degrees.");
    if (p < T{0} || p > T{360})
      throw std::invalid_argument(
          "SphereCoord: phi must be in [0, 360] degrees.");

    this->r = r;
    this->t = t;
    this->p = p;
    t_rad = t * T{constants::pi} / T{180};
    p_rad = p * T{constants::pi} / T{180};

    cost = std::cos(t_rad);
    sint = std::sin(t_rad);
    cosp = std::cos(p_rad);
    sinp = std::sin(p_rad);

    x = r * sint * cosp;
    y = r * sint * sinp;
    z = r * cost;

    rUnit = StaticVector<T, 3>{sint * cosp, sint * sinp, cost};
    tUnit = StaticVector<T, 3>{cost * cosp, cost * sinp, -sint};
    pUnit = StaticVector<T, 3>{-sinp, cosp, T{0}};
    position = r * rUnit;

    // Derivatives of spherical unit vectors:
    drUnitdt = tUnit;
    dtUnitdt = -rUnit;
    dpUnitdt = StaticVector<T, 3>{T{0}, T{0}, T{0}};
    drUnitdp = sint * pUnit;
    dtUnitdp = cost * pUnit;
    dpUnitdp = -sint * rUnit - cost * tUnit;
  }

  // Constructor from Cartesian coordinates.
  SphereCoord(T x, T y, T z, std::string_view coordSystem) {
    if (coordSystem != "Cartesian")
      throw std::invalid_argument(
          "SphereCoord: coordSystem must be Cartesian.");

    this->x = x;
    this->y = y;
    this->z = z;

    r = std::sqrt(x * x + y * y + z * z);
    if (approximatelyEqualAbsRel(r, T{0})) {
      t = T{0};
      p = T{0};
    } else {
      t = std::acos(z / r) * T{180} / T{constants::pi};
      p = std::atan2(y, x) * T{180} / T{constants::pi};
      if (p < T{0})
        p += T{360};
    }

    if (t < T{0} || t > T{180})
      throw std::invalid_argument(
          "SphereCoord: theta must be in [0, 180] degrees.");
    if (p < T{0} || p > T{360})
      throw std::invalid_argument(
          "SphereCoord: phi must be in [0, 360] degrees.");

    t_rad = t * T{constants::pi} / T{180};
    p_rad = p * T{constants::pi} / T{180};

    cost = std::cos(t_rad);
    sint = std::sin(t_rad);
    cosp = std::cos(p_rad);
    sinp = std::sin(p_rad);

    rUnit = StaticVector<T, 3>{sint * cosp, sint * sinp, cost};
    tUnit = StaticVector<T, 3>{cost * cosp, cost * sinp, -sint};
    pUnit = StaticVector<T, 3>{-sinp, cosp, T{0}};
    position = r * rUnit;

    drUnitdt = tUnit;
    dtUnitdt = -rUnit;
    dpUnitdt = StaticVector<T, 3>{T{0}, T{0}, T{0}};
    drUnitdp = sint * pUnit;
    dtUnitdp = cost * pUnit;
    dpUnitdp = -sint * rUnit - cost * tUnit;
  }

  T x{}, y{}, z{};
  T r{}, t{}, p{};
  T t_rad{}, p_rad{};
  T cost{}, sint{}, cosp{}, sinp{};
  StaticVector<T, 3> position{};
  StaticVector<T, 3> rUnit, tUnit{}, pUnit{};
  StaticVector<T, 3> drUnitdt, dtUnitdt{}, dpUnitdt{};
  StaticVector<T, 3> drUnitdp, dtUnitdp{}, dpUnitdp{};
  ScalarFunction3D f;
  VectorFunction3D u;

  // dxdydz = |J| drdtdp
  T jacobian() const { return r * r * sint; }

  StaticVector<T, 3> gradient(ScalarFunction3D f_rtp) const {
    if (approximatelyEqualAbsRel(r, T{0}))
      throw std::invalid_argument("SphereCoord::gradient requires r != 0.");
    if (approximatelyEqualAbsRel(sint, T{0}))
      throw std::invalid_argument(
          "SphereCoord::gradient undefined for sin(theta)=0.");

    Dual r_dr{r, 1.0};
    Dual t_dr{t, 0.0};
    Dual p_dr{p, 0.0};
    double dfdr = f_rtp(r_dr, t_dr, p_dr).getDer();

    Dual r_dt{r, 0.0};
    Dual t_dt{t, 1.0};
    Dual p_dt{p, 0.0};
    double dfdt = f_rtp(r_dt, t_dt, p_dt).getDer();

    Dual r_dp{r, 0.0};
    Dual t_dp{t, 0.0};
    Dual p_dp{p, 1.0};
    double dfdp = f_rtp(r_dp, t_dp, p_dp).getDer();

    return {dfdr, T{1} / r * dfdt, T{1} / (r * sint) * dfdp};
  }

  T div(VectorFunction3D u_rtp) const {
    if (approximatelyEqualAbsRel(r, T{0}))
      throw std::invalid_argument("SphereCoord::div requires r != 0.");
    if (approximatelyEqualAbsRel(sint, T{0}))
      throw std::invalid_argument(
          "SphereCoord::div undefined for sin(theta)=0.");

    auto r2ur = [u_rtp](Dual r_var, Dual t_var, Dual p_var) {
      auto val = u_rtp(r_var, t_var, p_var);
      return r_var * r_var * val[0];
    };
    double d_r2ur_dr = r2ur(Dual{r, 1.0}, Dual{t, 0.0}, Dual{p, 0.0}).getDer();

    auto sint_ut = [u_rtp](Dual r_var, Dual t_var, Dual p_var) {
      auto val = u_rtp(r_var, t_var, p_var);
      return sin(t_var * Dual{constants::pi / 180.0}) * val[1];
    };
    double d_sintut_dt =
        sint_ut(Dual{r, 0.0}, Dual{t, 1.0}, Dual{p, 0.0}).getDer();

    auto up = u_rtp(Dual{r, 0.0}, Dual{t, 0.0}, Dual{p, 1.0})[2];
    double dup_dp = up.getDer();

    return (T{1} / (r * r)) * d_r2ur_dr + (T{1} / (r * sint)) * d_sintut_dt +
           (T{1} / (r * sint)) * dup_dp;
  }

  StaticVector<T, 3> curl(VectorFunction3D u_rtp) const {
    if (approximatelyEqualAbsRel(r, T{0}))
      throw std::invalid_argument("SphereCoord::curl requires r != 0.");
    if (approximatelyEqualAbsRel(sint, T{0}))
      throw std::invalid_argument(
          "SphereCoord::curl undefined for sin(theta)=0.");

    auto sint_up = [u_rtp](Dual r_var, Dual t_var, Dual p_var) {
      auto val = u_rtp(r_var, t_var, p_var);
      return sin(t_var * Dual{constants::pi / 180.0}) * val[2];
    };
    double d_sintup_dt =
        sint_up(Dual{r, 0.0}, Dual{t, 1.0}, Dual{p, 0.0}).getDer();
    double dut_dp = u_rtp(Dual{r, 0.0}, Dual{t, 0.0}, Dual{p, 1.0})[1].getDer();
    T curl_r = (T{1} / (r * sint)) * (d_sintup_dt - dut_dp);

    auto r_up = [u_rtp](Dual r_var, Dual t_var, Dual p_var) {
      auto val = u_rtp(r_var, t_var, p_var);
      return r_var * val[2];
    };
    double d_rup_dr = r_up(Dual{r, 1.0}, Dual{t, 0.0}, Dual{p, 0.0}).getDer();
    double dur_dp = u_rtp(Dual{r, 0.0}, Dual{t, 0.0}, Dual{p, 1.0})[0].getDer();
    T curl_t = (T{1} / r) * ((T{1} / sint) * dur_dp - d_rup_dr);

    auto r_ut = [u_rtp](Dual r_var, Dual t_var, Dual p_var) {
      auto val = u_rtp(r_var, t_var, p_var);
      return r_var * val[1];
    };
    double d_rut_dr = r_ut(Dual{r, 1.0}, Dual{t, 0.0}, Dual{p, 0.0}).getDer();
    double dur_dt = u_rtp(Dual{r, 0.0}, Dual{t, 1.0}, Dual{p, 0.0})[0].getDer();
    T curl_p = (T{1} / r) * (d_rut_dr - dur_dt);

    return {curl_r, curl_t, curl_p};
  }

  T laplacian(ScalarFunction3D f_rtp, T hr = T{1e-5}, T ht = T{1e-5},
              T hp = T{1e-5}) const {
    if (approximatelyEqualAbsRel(r, T{0}))
      throw std::invalid_argument("SphereCoord::laplacian requires r != 0.");
    if (approximatelyEqualAbsRel(sint, T{0}))
      throw std::invalid_argument(
          "SphereCoord::laplacian undefined for sin(theta)=0.");
    if (hr <= T{0} || ht <= T{0} || hp <= T{0})
      throw std::invalid_argument(
          "SphereCoord::laplacian requires positive steps.");

    auto eval = [&](T r_eval, T t_eval, T p_eval) -> T {
      return static_cast<T>(f_rtp(Dual{static_cast<double>(r_eval), 0.0},
                                  Dual{static_cast<double>(t_eval), 0.0},
                                  Dual{static_cast<double>(p_eval), 0.0})
                                .getVal());
    };

    const T f0 = eval(r, t, p);
    const T fr_p = eval(r + hr, t, p);
    const T fr_m = eval(r - hr, t, p);
    const T ft_p = eval(r, t + ht, p);
    const T ft_m = eval(r, t - ht, p);
    const T fp_p = eval(r, t, p + hp);
    const T fp_m = eval(r, t, p - hp);

    const T d2fdr2 = (fr_p - T{2} * f0 + fr_m) / (hr * hr);
    const T dfdr = (fr_p - fr_m) / (T{2} * hr);
    const T d2fdt2 = (ft_p - T{2} * f0 + ft_m) / (ht * ht);
    const T dfdt = (ft_p - ft_m) / (T{2} * ht);
    const T d2fdp2 = (fp_p - T{2} * f0 + fp_m) / (hp * hp);

    const T cot_t = cost / sint;
    return d2fdr2 + (T{2} / r) * dfdr +
           (T{1} / (r * r)) *
               (d2fdt2 + cot_t * dfdt + (T{1} / (sint * sint)) * d2fdp2);
  }

  StaticVector<T, 3> laplacian(VectorFunction3D u_rtp, T hr = T{1e-5},
                               T ht = T{1e-5}, T hp = T{1e-5}) const {
    auto u1 = [u_rtp](Dual r_var, Dual t_var, Dual p_var) -> Dual {
      return u_rtp(r_var, t_var, p_var)[0];
    };
    auto u2 = [u_rtp](Dual r_var, Dual t_var, Dual p_var) -> Dual {
      return u_rtp(r_var, t_var, p_var)[1];
    };
    auto u3 = [u_rtp](Dual r_var, Dual t_var, Dual p_var) -> Dual {
      return u_rtp(r_var, t_var, p_var)[2];
    };
    return {laplacian(u1, hr, ht, hp), laplacian(u2, hr, ht, hp),
            laplacian(u3, hr, ht, hp)};
  }
};

#endif