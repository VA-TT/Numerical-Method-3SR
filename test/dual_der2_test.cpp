#include "../library/DualDifferentiation.h"
#include <cassert>
#include <cmath>
#include <iostream>

namespace {

bool almostEqual(double a, double b, double tol = 1e-9) {
  return std::abs(a - b) <= tol;
}

} // namespace

int main() {
  {
    // f(x) = x^3 => f''(x) = 6x => f''(2) = 12
    auto f = [](Dual x) -> Dual { return x * x * x; };
    const double d2 = automaticDiff2(f, 2.0);
    assert(almostEqual(d2, 12.0));
  }

  {
    // f(x) = sin(x) => f''(x) = -sin(x)
    const double x0 = std::numbers::pi / 6.0;
    auto f = [](Dual x) -> Dual { return sin(x); };
    const double d2 = automaticDiff2(f, x0);
    assert(almostEqual(d2, -std::sin(x0), 1e-10));
  }

  {
    // Seed directly with Dual(val, der1, der2)
    // g(x) = x^4 => g''(x)=12x^2 => g''(2)=48
    Dual x{2.0, 1.0, 0.0};
    Dual g = x * x * x * x;
    assert(almostEqual(g.getDer2(), 48.0));
  }

  {
    // Scalar Laplacian of x^2+y^2+z^2 is 6.
    auto f = [](Dual x, Dual y, Dual z) -> Dual {
      return x * x + y * y + z * z;
    };
    const double lap = laplacian3D(f, 1.25, -2.0, 0.3);
    assert(almostEqual(lap, 6.0));
  }

  {
    // Vector Laplacian of [x^2, y^2, z^2] is [2,2,2].
    auto u = [](Dual x, Dual y, Dual z) -> StaticVector<Dual, 3> {
      return {x * x, y * y, z * z};
    };
    const auto lapU = laplacian3D(u, 0.5, -1.0, 3.0);
    assert(almostEqual(lapU[0], 2.0));
    assert(almostEqual(lapU[1], 2.0));
    assert(almostEqual(lapU[2], 2.0));
  }

  std::cout << "dual_der2_test passed\n";
  return 0;
}
