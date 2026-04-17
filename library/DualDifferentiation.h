#ifndef DUAL_DIFFERENTIATION_H
#define DUAL_DIFFERENTIATION_H

#include "Matrix.h"
#include "Vector.h"
#include <cmath>
#include <functional>
#include <iostream>

// Define an epsilon that epsilon^3 = 0
// Dual expansion: f(x+epsilon) = f(x) + epsilon.f'(x) + epsilon^2.f''(x)
//                              = a    + epsilon.b     + epsilon.c
//  => a = f(x); b = f'(x); c = f''(x)
// Note: Must seed Dual x{x0,1.0} / x{x0,1.0,0.0}
class Dual {
private:
  double m_val{1.0};
  double m_der{1.0};  // first-order derivative
  double m_der2{0.0}; // second-order derivative

public:
  Dual(double val, double der, double der2)
      : m_val{val}, m_der{der}, m_der2{der2} {}
  Dual(double val, double der) : m_val{val}, m_der{der}, m_der2{0.0} {}
  Dual(double val) : m_val{val}, m_der{0.0}, m_der2{0.0} {}
  Dual() = default;
  Dual(const Dual &) = default;
  Dual(Dual &&) = default;
  Dual &operator=(const Dual &) = default;
  Dual &operator=(Dual &&) = default;
  ~Dual() = default;

  constexpr double getVal() const { return m_val; }
  constexpr double getDer() const { return m_der; }
  constexpr double getDer2() const { return m_der2; }
  void setVal(double val) { m_val = val; }
  void setDer(double der) { m_der = der; }
  void setDer2(double der2) { m_der2 = der2; }
  friend Dual operator*(double k, const Dual &d) {
    Dual result{d};
    result.m_val *= k;
    result.m_der *= k;
    result.m_der2 *= k;
    return result;
  }
  friend Dual operator*(const Dual &d, double k) { return k * d; }
  friend Dual operator+(const Dual &d, double a) {
    Dual result{d};
    result.m_val += a;
    return result;
  }
  friend Dual operator+(double a, const Dual &d) { return d + a; }
  friend Dual operator+(const Dual &d1, const Dual &d2) {
    Dual result{d1};
    result.m_val += d2.m_val;
    result.m_der += d2.m_der;
    result.m_der2 += d2.m_der2;
    return result;
  }
  friend Dual operator-(const Dual &d1, const Dual &d2) { return d1 + (-d2); }
  friend Dual operator-(const Dual &d) {
    return Dual{-d.m_val, -d.m_der, -d.m_der2};
  }
  friend Dual operator*(const Dual &d1, const Dual &d2) {
    Dual result{};
    result.m_val = d1.m_val * d2.m_val;
    result.m_der = d1.m_val * d2.m_der + d1.m_der * d2.m_val;
    result.m_der2 =
        d1.m_der2 * d2.m_val + 2.0 * d1.m_der * d2.m_der + d1.m_val * d2.m_der2;
    return result;
  }

  friend Dual inverse(const Dual &d) {
    const double inv = 1.0 / d.m_val;
    const double inv2 = inv * inv;
    const double inv3 = inv2 * inv;
    return Dual{inv, -d.m_der * inv2,
                (2.0 * d.m_der * d.m_der - d.m_val * d.m_der2) * inv3};
  }

  friend Dual operator/(const Dual &d1, const Dual &d2) {
    return d1 * inverse(d2);
  }

  // Comparison operators - compare value parts only
  friend bool operator<(const Dual &d, double b) { return d.m_val < b; }
  friend bool operator<(double a, const Dual &d) { return a < d.m_val; }
  friend bool operator>(const Dual &d, double b) { return d.m_val > b; }
  friend bool operator>(double a, const Dual &d) { return a > d.m_val; }
  friend bool operator<=(const Dual &d, double b) { return d.m_val <= b; }
  friend bool operator<=(double a, const Dual &d) { return a <= d.m_val; }
  friend bool operator>=(const Dual &d, double b) { return d.m_val >= b; }
  friend bool operator>=(double a, const Dual &d) { return a >= d.m_val; }
  friend bool operator==(const Dual &d, double b) { return d.m_val == b; }
  friend bool operator==(double a, const Dual &d) { return a == d.m_val; }
  friend bool operator!=(const Dual &d, double b) { return d.m_val != b; }
  friend bool operator!=(double a, const Dual &d) { return a != d.m_val; }

  // Dual vs Dual comparisons
  friend bool operator<(const Dual &d1, const Dual &d2) {
    return d1.m_val < d2.m_val;
  }
  friend bool operator>(const Dual &d1, const Dual &d2) {
    return d1.m_val > d2.m_val;
  }
  friend bool operator<=(const Dual &d1, const Dual &d2) {
    return d1.m_val <= d2.m_val;
  }
  friend bool operator>=(const Dual &d1, const Dual &d2) {
    return d1.m_val >= d2.m_val;
  }
  friend bool operator==(const Dual &d1, const Dual &d2) {
    return d1.m_val == d2.m_val;
  }
  friend bool operator!=(const Dual &d1, const Dual &d2) {
    return d1.m_val != d2.m_val;
  }

  // Compound assignment operators
  Dual &operator+=(const Dual &d) {
    m_val += d.m_val;
    m_der += d.m_der;
    m_der2 += d.m_der2;
    return *this;
  }
  Dual &operator+=(double a) {
    m_val += a;
    return *this;
  }
  Dual &operator-=(const Dual &d) {
    m_val -= d.m_val;
    m_der -= d.m_der;
    m_der2 -= d.m_der2;
    return *this;
  }
  Dual &operator-=(double a) {
    m_val -= a;
    return *this;
  }
  Dual &operator*=(const Dual &d) {
    // (a + ae*E + ae2*E^2) * (b + be*E + be2*E^2)
    m_der2 = m_der2 * d.m_val + 2.0 * m_der * d.m_der + m_val * d.m_der2;
    m_der = m_val * d.m_der + m_der * d.m_val;
    m_val = m_val * d.m_val;
    return *this;
  }
  Dual &operator*=(double k) {
    m_val *= k;
    m_der *= k;
    m_der2 *= k;
    return *this;
  }
  Dual &operator/=(const Dual &d) {
    *this = *this / d;
    return *this;
  }
  Dual &operator/=(double k) {
    m_val /= k;
    m_der /= k;
    m_der2 /= k;
    return *this;
  }

  friend std::ostream &operator<<(std::ostream &out, const Dual &d) {
    out << d.m_val << " + " << d.m_der << "E + " << d.m_der2 << "E2\n";
    return out;
  }
};

Dual sin(const Dual &d) {
  const double s = std::sin(d.getVal());
  const double c = std::cos(d.getVal());
  return Dual{s, c * d.getDer(),
              -s * d.getDer() * d.getDer() + c * d.getDer2()};
}
Dual cos(const Dual &d) {
  const double s = std::sin(d.getVal());
  const double c = std::cos(d.getVal());
  return Dual{c, -s * d.getDer(),
              -c * d.getDer() * d.getDer() - s * d.getDer2()};
}

inline Dual log(const Dual &d) {
  // d/dx log(x) = 1/x
  const double inv = 1.0 / d.getVal();
  return Dual{std::log(d.getVal()), d.getDer() * inv,
              d.getDer2() * inv - d.getDer() * d.getDer() * inv * inv};
}

inline Dual pow(const Dual &d, double n) {
  // d/dx x^n = n*x^(n-1)
  const double v_n = std::pow(d.getVal(), n);
  const double v_n1 = std::pow(d.getVal(), n - 1.0);
  const double v_n2 = std::pow(d.getVal(), n - 2.0);
  return Dual{v_n, n * v_n1 * d.getDer(),
              n * (n - 1.0) * v_n2 * d.getDer() * d.getDer() +
                  n * v_n1 * d.getDer2()};
}

inline Dual pow(const Dual &d, int n) {
  // d/dx x^n = n*x^(n-1)
  return pow(d, static_cast<double>(n));
}

inline Dual operator-(const Dual &d, double b) {
  return Dual{d.getVal() - b, d.getDer(), d.getDer2()};
}
inline Dual operator-(double a, const Dual &d) {
  return Dual{a - d.getVal(), -d.getDer(), -d.getDer2()};
}
inline Dual operator/(double a, const Dual &d) { return Dual{a} / d; }

// Generic automatic differentiation
template <typename Func, typename T = double>
auto automaticDiff(Func func, T x0) -> double {
  Dual d{static_cast<double>(x0), 1.0};
  auto result = func(d);
  return static_cast<double>(result.getDer());
}

template <typename Func, typename T = double>
auto automaticDiff2(Func func, T x0) -> double {
  Dual d{static_cast<double>(x0), 1.0, 0.0};
  auto result = func(d);
  return static_cast<double>(result.getDer2());
}

// Compute gradient of a scalar function f(Dual x,Dual y) ->vector
using ScalarFunction2D = std::function<Dual(Dual, Dual)>;
StaticVector<double, 2> grad2D(const ScalarFunction2D &f,
                               const StaticVector<double, 2> &initValue) {
  const double x0 = initValue[0];
  const double y0 = initValue[1];

  Dual r_x = f({x0, 1.0}, {y0, 0.0});
  double dfdx = r_x.getDer();

  Dual r_y = f({x0, 0.0}, {y0, 1.0});
  double dfdy = r_y.getDer();

  return {dfdx, dfdy};
}

// Compute gradient of vector function -> tensor aka Jacobian matrix
using VectorFunction2D = std::function<StaticVector<Dual, 2>(Dual, Dual)>;
Matrix<double, 2, 2> grad2D(const VectorFunction2D &u,
                            const StaticVector<double, 2> &initValue) {
  const double x0 = initValue[0];
  const double y0 = initValue[1];

  auto u_dualX = u({x0, 1.0}, {y0, 0.0});
  double du1dx = u_dualX.x().getDer();
  double du2dx = u_dualX.y().getDer();

  auto u_dualY = u({x0, 0.0}, {y0, 1.0});
  double du1dy = u_dualY.x().getDer();
  double du2dy = u_dualY.y().getDer();

  return {du1dx, du1dy, du2dx, du2dy};
}

Matrix<double, 2, 2> Jacobian(const VectorFunction2D &u,
                              const StaticVector<double, 2> &initValue) {
  return grad2D(u, initValue);
}

double div2D(const VectorFunction2D &u,
             const StaticVector<double, 2> &initValue) {
  const double x0 = initValue[0];
  const double y0 = initValue[1];

  auto u_dualX = u({x0, 1.0}, {y0, 0.0});

  auto u_dualY = u({x0, 0.0}, {y0, 1.0});

  return u_dualX.x().getDer() + u_dualY.y().getDer();
}

// 2D curl = scalar z-component of the 3D curl:
double curl2D(const VectorFunction2D &u,
              const StaticVector<double, 2> &initValue) {
  const double x0 = initValue[0];
  const double y0 = initValue[1];

  auto u_dualX = u({x0, 1.0}, {y0, 0.0});

  auto u_dualY = u({x0, 0.0}, {y0, 1.0});

  return u_dualX.y().getDer() - u_dualY.x().getDer();
}

using ScalarFunction3D = std::function<Dual(Dual, Dual, Dual)>;
StaticVector<double, 3> grad3D(const ScalarFunction3D &f,
                               const StaticVector<double, 3> &initValue) {
  const double x0 = initValue[0];
  const double y0 = initValue[1];
  const double z0 = initValue[2];

  // df/dx
  Dual r_x = f({x0, 1.0}, {y0, 0.0}, {z0, 0.0});
  double dfdx = r_x.getDer();

  // df/dy
  Dual r_y = f({x0, 0.0}, {y0, 1.0}, {z0, 0.0});
  double dfdy = r_y.getDer();

  // df/dz
  Dual r_z = f({x0, 0.0}, {y0, 0.0}, {z0, 1.0});
  double dfdz = r_z.getDer();

  return {dfdx, dfdy, dfdz};
}

using VectorFunction3D = std::function<StaticVector<Dual, 3>(Dual, Dual, Dual)>;

Matrix<double, 3, 3> grad3D(const VectorFunction3D &u,
                            const StaticVector<double, 3> &initValue) {
  const double x0 = initValue[0];
  const double y0 = initValue[1];
  const double z0 = initValue[2];

  auto u_dualX = u({x0, 1.0}, {y0, 0.0}, {z0, 0.0});
  double du1dx = u_dualX.x().getDer();
  double du2dx = u_dualX.y().getDer();
  double du3dx = u_dualX.z().getDer();

  auto u_dualY = u({x0, 0.0}, {y0, 1.0}, {z0, 0.0});
  double du1dy = u_dualY.x().getDer();
  double du2dy = u_dualY.y().getDer();
  double du3dy = u_dualY.z().getDer();

  auto u_dualZ = u({x0, 0.0}, {y0, 0.0}, {z0, 1.0});
  double du1dz = u_dualZ.x().getDer();
  double du2dz = u_dualZ.y().getDer();
  double du3dz = u_dualZ.z().getDer();

  return {du1dx, du1dy, du1dz, du2dx, du2dy, du2dz, du3dx, du3dy, du3dz};
}

Matrix<double, 3, 3> Jacobian(const VectorFunction3D &u,
                              const StaticVector<double, 3> &initValue) {
  return grad3D(u, initValue);
}

double div3D(const VectorFunction3D &u,
             const StaticVector<double, 3> &initValue) {
  const double x0 = initValue[0];
  const double y0 = initValue[1];
  const double z0 = initValue[2];

  // u_x,x
  auto u_dualX = u({x0, 1.0}, {y0, 0.0}, {z0, 0.0});
  double du1dx = u_dualX.x().getDer();

  // u_y,y
  auto u_dualY = u({x0, 0.0}, {y0, 1.0}, {z0, 0.0});
  double du2dy = u_dualY.y().getDer();

  // u_z,z
  auto u_dualZ = u({x0, 0.0}, {y0, 0.0}, {z0, 1.0});
  double du3dz = u_dualZ.z().getDer();

  return du1dx + du2dy + du3dz;
}

StaticVector<double, 3> curl3D(const VectorFunction3D &u,
                               const StaticVector<double, 3> &initValue) {
  const double x0 = initValue[0];
  const double y0 = initValue[1];
  const double z0 = initValue[2];

  // du_i/dx
  auto u_dualX = u({x0, 1.0}, {y0, 0.0}, {z0, 0.0});
  double du2dx = u_dualX.y().getDer();
  double du3dx = u_dualX.z().getDer();

  // du_i/dy
  auto u_dualY = u({x0, 0.0}, {y0, 1.0}, {z0, 0.0});
  double du1dy = u_dualY.x().getDer();
  double du3dy = u_dualY.z().getDer();

  // du_i/dz
  auto u_dualZ = u({x0, 0.0}, {y0, 0.0}, {z0, 1.0});
  double du1dz = u_dualZ.x().getDer();
  double du2dz = u_dualZ.y().getDer();

  return {du3dy - du2dz, du1dz - du3dx, du2dx - du1dy};
}

// Laplacian scalar
double laplacian3D(const ScalarFunction3D &f,
                   const StaticVector<double, 3> &initValue) {
  const double x0 = initValue[0];
  const double y0 = initValue[1];
  const double z0 = initValue[2];

  const double d2fdx2 = f({x0, 1.0}, {y0, 0.0}, {z0, 0.0}).getDer2();
  const double d2fdy2 = f({x0, 0.0}, {y0, 1.0}, {z0, 0.0}).getDer2();
  const double d2fdz2 = f({x0, 0.0}, {y0, 0.0}, {z0, 1.0}).getDer2();

  return d2fdx2 + d2fdy2 + d2fdz2;
}
// Laplacian vector
StaticVector<double, 3> laplacian3D(const VectorFunction3D &u,
                                    const StaticVector<double, 3> &initValue) {
  const double x0 = initValue[0];
  const double y0 = initValue[1];
  const double z0 = initValue[2];

  StaticVector<double, 3> uLap;
  for (int i = 0; i < 3; i++) {
    auto ui = [u, i](Dual x, Dual y, Dual z) -> Dual { return u(x, y, z)[i]; };
    uLap[i] = laplacian3D(ui, {x0, y0, z0}); // laplacian scalar
  }
  return uLap;
}

#endif

#if 0
int main() {
  double x0{std::numbers::pi * 0.5};
  std::cout << "f(x0) = " << myFunc(x0) << '\n';
  std::cout << "f'(x0) = " << automaticDiff(x0) << '\n';
  auto F = [](Dual x, Dual y, Dual z) {
    return StaticVector<Dual, 3>{
        x * x, // Fx
        y * y, // Fy
        z * z  // Fz
    };
  };

  double div = div3D(F, {1.0, 2.0, 3.0});

  auto u = [](Dual x, Dual y, Dual z) {
    return StaticVector<Dual, 3>{y, -x, z * z};
  };

  auto c = curl3D(u, {1.0, 2.0, 3.0});
  return 0;
}
#endif
