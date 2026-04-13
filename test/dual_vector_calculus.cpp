#include "../library/DualDifferentiation.h"
#include <cassert>
#include <cmath>
#include <iostream>

int main() {
  const double tol = 1e-7;

  {
    auto f2 = [](Dual x, Dual y) -> Dual { return x * x + x * y + y * y; };
    auto g = grad2D(f2, 1.0, 2.0);
    assert(std::abs(g.x() - 4.0) < tol);
    assert(std::abs(g.y() - 5.0) < tol);
  }

  {
    VectorFunction2D v2 = [](Dual x, Dual y) -> StaticVector<Dual, 2> {
      return {x * x, x * y};
    };
    const double div = div2D(v2, 1.0, 2.0);
    const double curl = curl2D(v2, 1.0, 2.0);
    assert(std::abs(div - 3.0) < tol);  // 2x + x = 3 at x=1
    assert(std::abs(curl - 2.0) < tol); // d(xy)/dx - d(x^2)/dy = y - 0
  }

  {
    ScalarFunction3D f3 = [](Dual x, Dual y, Dual z) -> Dual {
      return x * x + y * y + z * z;
    };
    auto g3 = grad3D(f3, 1.0, 2.0, 3.0);
    assert(std::abs(g3.x() - 2.0) < tol);
    assert(std::abs(g3.y() - 4.0) < tol);
    assert(std::abs(g3.z() - 6.0) < tol);

    double lap = laplacian3D(f3, 1.0, 2.0, 3.0);
    assert(std::abs(lap - 6.0) < 5e-4);
  }

  {
    VectorFunction3D u = [](Dual x, Dual y, Dual z) -> StaticVector<Dual, 3> {
      return {x + y, y + z, z + x};
    };
    const auto J = grad3D(u, 1.0, 2.0, 3.0);
    assert(std::abs(J(0, 0) - 1.0) < tol);
    assert(std::abs(J(0, 1) - 1.0) < tol);
    assert(std::abs(J(0, 2) - 0.0) < tol);
    assert(std::abs(J(1, 0) - 0.0) < tol);
    assert(std::abs(J(1, 1) - 1.0) < tol);
    assert(std::abs(J(1, 2) - 1.0) < tol);
    assert(std::abs(J(2, 0) - 1.0) < tol);
    assert(std::abs(J(2, 1) - 0.0) < tol);
    assert(std::abs(J(2, 2) - 1.0) < tol);
  }

  {
    VectorFunction3D u = [](Dual x, Dual y, Dual z) -> StaticVector<Dual, 3> {
      return {x, y, z};
    };
    const double div = div3D(u, 0.5, -1.2, 4.0);
    assert(std::abs(div - 3.0) < tol);
  }

  {
    VectorFunction3D v = [](Dual x, Dual y, Dual z) -> StaticVector<Dual, 3> {
      (void)z;
      return {y, -x, x * y};
    };
    auto c = curl3D(v, 1.5, -2.0, 0.0);
    assert(std::abs(c.x() - 1.5) < tol);
    assert(std::abs(c.y() - 2.0) < tol);
    assert(std::abs(c.z() + 2.0) < tol);
  }

  {
    VectorFunction3D v = [](Dual x, Dual y, Dual z) -> StaticVector<Dual, 3> {
      return {x * x, y * y, z * z};
    };
    auto l = laplacian3D(v, 2.0, -3.0, 1.5);
    assert(std::abs(l.x() - 2.0) < 5e-4);
    assert(std::abs(l.y() - 2.0) < 5e-4);
    assert(std::abs(l.z() - 2.0) < 5e-4);
  }

  std::cout << "Dual vector calculus tests passed.\n";
  return 0;
}
