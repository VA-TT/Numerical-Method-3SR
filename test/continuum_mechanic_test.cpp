#include "../library/ContinuumMechanic.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

bool near(double a, double b, double tol = 1e-10) {
  return std::abs(a - b) <= tol;
}

void test_homogeneous_2d() {
  using MC = MediumCon<double, 2>;

  MC::MapMatrix K = [](Dual t) {
    Matrix<Dual, 2, 2> out{};
    out(0, 0) = Dual{1.0, 0.0} + t;
    out(0, 1) = Dual{0.0, 0.0};
    out(1, 0) = Dual{0.0, 0.0};
    out(1, 1) = Dual{1.0, 0.0} - t;
    return out;
  };

  StaticVector<double, 2> X{2.0, 3.0};
  MC mc(K, X);

  const double t = 0.2;

  const auto x = mc.currentPosition(t);
  assert(near(x.x(), 2.4));
  assert(near(x.y(), 2.4));

  const auto v = mc.velocityMaterial(t);
  assert(near(v.x(), 2.0));
  assert(near(v.y(), -3.0));

  const auto a = mc.accelerationMaterial(t);
  assert(near(a.x(), 0.0));
  assert(near(a.y(), 0.0));

  const auto F = mc.deformationGradient(t);
  assert(near(F(0, 0), 1.2));
  assert(near(F(0, 1), 0.0));
  assert(near(F(1, 0), 0.0));
  assert(near(F(1, 1), 0.8));

  assert(near(mc.Jacobian(t), 0.96));

  const auto u = mc.displacementMaterial(t);
  assert(near(u.x(), 0.4));
  assert(near(u.y(), -0.6));
}

void test_non_homogeneous_2d() {
  using MC = MediumCon<double, 2>;

  MC::MotionMap phi = [](const StaticVector<Dual, 2> &X, Dual t) {
    StaticVector<Dual, 2> out{};
    out.x() = X.x() + t * X.x() * X.y();
    out.y() = X.y() + t * t * X.x();
    return out;
  };

  StaticVector<double, 2> X{1.5, 2.0};
  MC mc(phi, X);

  const double t = 0.3;

  const auto x = mc.currentPosition(t);
  assert(near(x.x(), 2.4));
  assert(near(x.y(), 2.135));

  const auto v = mc.velocityMaterial(t);
  assert(near(v.x(), 3.0));
  assert(near(v.y(), 0.9));

  const auto a = mc.accelerationMaterial(t);
  assert(near(a.x(), 0.0));
  assert(near(a.y(), 3.0));

  const auto F = mc.deformationGradient(t, X);
  assert(near(F(0, 0), 1.6));
  assert(near(F(0, 1), 0.45));
  assert(near(F(1, 0), 0.09));
  assert(near(F(1, 1), 1.0));

  bool threw_inverse = false;
  try {
    (void)mc.inverseMap(t);
  } catch (const std::logic_error &) {
    threw_inverse = true;
  }
  assert(threw_inverse);
}

} // namespace

int main() {
  test_homogeneous_2d();
  test_non_homogeneous_2d();
  std::cout << "ContinuumMechanic tests passed\n";
  return 0;
}
