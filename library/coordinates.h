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

void changeCoord() {}

// 2D polar coordinate
template <typename T> struct polarCoor {
  // Constructor using (r,t) variables
  polarCoor(T r, T t) // theta is in degree unit
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
  polarCoor(T x, T y, std::string_view coordSystem)
      : polarCoor(std::sqrt(x * x + y * y),
                  std::atan2(y, x) * T{180} / T{constants::pi}) {
    if (coordSystem != "Cartesian")
      throw std::invalid_argument("polarCoor: coordSystem must be Cartesian.");

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
      throw std::invalid_argument("polarCoor::gradient requires r != 0.");

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
      throw std::invalid_argument("polarCoor::div requires r != 0.");

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
      throw std::invalid_argument("polarCoor::curl requires r != 0.");

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
      throw std::invalid_argument("polarCoor::laplacian requires r != 0.");
    if (hr <= T{0} || ht <= T{0})
      throw std::invalid_argument(
          "polarCoor::laplacian requires positive steps.");

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
template <typename T> struct cylinCoor {
  // Cylindrical system constructors
  cylinCoor(T r, T p, T z) // phi is in degree unit
  {
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

  // Constructor using Cartesian variables.
  cylinCoor(T x, T y, T z, std::string_view coordSystem)
      : cylinCoor(std::sqrt(x * x + y * y),
                  std::atan2(y, x) * T{180} / T{constants::pi}, z) {
    if (coordSystem != "Cartesian")
      throw std::invalid_argument("cylinCoor: coordSystem must be Cartesian.");

    this->x = x;
    this->y = y;
    this->z = z;
  }

  // Member variables
  T x{}, y{};
  T z{};      // Shared component between Cartesian and cylindrical
  T r{}, p{}; // Rho and Phi - Equivalent to r and theta in 2D
  T p_rad{};
  T cosp{}, sinp{};
  StaticVector<T, 3> position{};
  StaticVector<T, 3> rUnit, pUnit{}, zUnit{};
  ScalarFunction3D f; // f = f(r,p,z)
  VectorFunction3D u; // u = (u1(r,p,z), u2(r,p,z), u3(r,p,z))

  // Member functions (differential operators)
  T jacobian() const { return r; } // dx.dy.dz = |J| dr.dp.dz

  StaticVector<T, 3> gradient(ScalarFunction3D f_rpz) const {
    if (approximatelyEqualAbsRel(r, T{0}))
      throw std::invalid_argument("cylinCoor::gradient requires rho != 0.");

    Dual r_dr{r, 1.0};
    Dual p_dr{p, 0.0};
    Dual z_dr{z, 0.0};
    auto f_dualR = f_rpz(r_dr, p_dr, z_dr);
    double dfdr = f_dualR.getDer();

    Dual r_dp{r, 0.0};
    Dual p_dp{p, 1.0};
    Dual z_dp{z, 0.0};
    auto f_dualP = f_rpz(r_dp, p_dp, z_dp);
    double dfdp = f_dualP.getDer();

    Dual r_dz{r, 0.0};
    Dual p_dz{p, 0.0};
    Dual z_dz{z, 1.0};
    auto f_dualZ = f_rpz(r_dz, p_dz, z_dz);
    double dfdz = f_dualZ.getDer();

    return {dfdr, T{1} / r * dfdp, dfdz}; // cylindrical coordinate
  }

  // Member functions (differential operators)
  T div(VectorFunction3D u_rpz) const {
    if (approximatelyEqualAbsRel(r, T{0}))
      throw std::invalid_argument("cylinCoor::div requires r != 0.");

    Dual r_dr{r, 1.0};
    Dual p_dr{p, 0.0};
    Dual z_dr{z, 0.0};
    auto u_dualR = u_rpz(r_dr, p_dr, z_dr);
    double du1dr = u_dualR[0].getDer();
    double u1 = u_dualR[0].getVal();

    Dual r_dp{r, 0.0};
    Dual p_dp{p, 1.0};
    Dual z_dp{z, 0.0};
    auto u_dualP = u_rpz(r_dp, p_dp, z_dp);
    double du2dp = u_dualP[1].getDer();

    Dual r_dz{r, 0.0};
    Dual p_dz{p, 0.0};
    Dual z_dz{z, 1.0};
    auto u_dualZ = u_rpz(r_dz, p_dz, z_dz);
    double du3dz = u_dualZ[2].getDer();

    return du1dr + T{1} / r * (u1 + du2dp) + du3dz;
  }

  // Cylindrical curl (r, p, z components).
  StaticVector<T, 3> curl(VectorFunction3D u_rpz) const {
    if (approximatelyEqualAbsRel(r, T{0}))
      throw std::invalid_argument("cylinCoor::curl requires r != 0.");

    Dual r_dr{r, 1.0};
    Dual p_dr{p, 0.0};
    Dual z_dr{z, 0.0};
    auto u_dualR = u_rpz(r_dr, p_dr, z_dr);
    double du2dr = u_dualR[1].getDer();
    double du3dr = u_dualR[2].getDer();

    Dual r_dp{r, 0.0};
    Dual p_dp{p, 1.0};
    Dual z_dp{z, 0.0};
    auto u_dualP = u_rpz(r_dp, p_dp, z_dp);
    double du1dp = u_dualP[0].getDer();
    double du3dp = u_dualP[2].getDer();
    double u2 = u_dualP[1].getVal();

    Dual r_dz{r, 0.0};
    Dual p_dz{p, 0.0};
    Dual z_dz{z, 1.0};
    auto u_dualZ = u_rpz(r_dz, p_dz, z_dz);
    double du1dz = u_dualZ[0].getDer();
    double du2dz = u_dualZ[1].getDer();

    return {T{1} / r * du3dp - du2dz, du1dz - du3dr,
            du2dr + T{1} / r * (u2 - du1dp)};
  }

  // Scalar Laplacian in cylindrical coordinates:
  // nabla^2 f = d2f/dr2 + (1/r)df/dr + (1/r^2)d2f/dp2 + d2f/dz2
  T laplacian(ScalarFunction3D f_rpz, T hr = T{1e-5}, T hp = T{1e-5},
              T hz = T{1e-5}) const {
    if (approximatelyEqualAbsRel(r, T{0}))
      throw std::invalid_argument("cylinCoor::laplacian requires r != 0.");
    if (hr <= T{0} || hp <= T{0} || hz <= T{0})
      throw std::invalid_argument(
          "cylinCoor::laplacian requires positive steps.");

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

  // Vector Laplacian in cylindrical coordinates (component-wise).
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

#endif