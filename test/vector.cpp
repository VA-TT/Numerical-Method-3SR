#include "../library/Vector.h"
#include <cassert>
#include <iostream>
#include <stdexcept>

int main() {
  StaticVector<double, 3> v1{1.0, 2.0, 3.0};
  StaticVector<double, 3> v2{4.0, 5.0, 6.0};

  // Test truy cập phần tử
  assert(v1[0] == 1.0 && v1[1] == 2.0 && v1[2] == 3.0);

  // Test toán tử +
  auto v3 = v1 + v2;
  assert(v3[0] == 5.0 && v3[1] == 7.0 && v3[2] == 9.0);

  // Test toán tử *
  auto v4 = 2.0 * v1;
  assert(v4[0] == 2.0 && v4[1] == 4.0 && v4[2] == 6.0);

  // Test dot product
  double dot = dotProduct(v1, v2);
  assert(dot == 32.0);

  // Test cross product
  auto v5 = crossProduct(v1, v2);
  assert(v5[0] == -3.0 && v5[1] == 6.0 && v5[2] == -3.0);

  // Test normalize
  auto v6 = normalize(v1);
  double mag = magnitude(v6);
  assert(std::abs(mag - 1.0) < 1e-12);

  // ===== New feature tests: Dynamic <-> Static conversion =====
  DynamicVector<double> dv1{7.0, 8.0, 9.0};

  // Dynamic -> Static (explicit constructor)
  StaticVector<double, 3> sv_from_dv{dv1};
  assert(sv_from_dv[0] == 7.0 && sv_from_dv[1] == 8.0 && sv_from_dv[2] == 9.0);

  // Static -> Dynamic (explicit constructor)
  DynamicVector<double> dv_from_sv{v1};
  assert(dv_from_sv.size() == 3);
  assert(dv_from_sv[0] == 1.0 && dv_from_sv[1] == 2.0 && dv_from_sv[2] == 3.0);

  // Static -> Dynamic (assignment)
  DynamicVector<double> dv_assign;
  dv_assign = v2;
  assert(dv_assign.size() == 3);
  assert(dv_assign[0] == 4.0 && dv_assign[1] == 5.0 && dv_assign[2] == 6.0);

  // Dynamic -> Static (assignment)
  StaticVector<double, 3> sv_assign{};
  sv_assign = dv1;
  assert(sv_assign[0] == 7.0 && sv_assign[1] == 8.0 && sv_assign[2] == 9.0);

  // Static -> Dynamic (explicit cast operator)
  DynamicVector<double> dv_cast = static_cast<DynamicVector<double>>(v1);
  assert(dv_cast.size() == 3);
  assert(dv_cast[0] == 1.0 && dv_cast[1] == 2.0 && dv_cast[2] == 3.0);

  // Runtime size validation (Dynamic -> Static mismatch must throw)
  bool mismatchCaught = false;
  try {
    DynamicVector<double> bad{1.0, 2.0};
    StaticVector<double, 3> shouldThrow{bad};
    (void)shouldThrow;
  } catch (const std::invalid_argument &) {
    mismatchCaught = true;
  }
  assert(mismatchCaught && "Expected size-mismatch conversion to throw");

  std::cout << "All vector conversion tests passed!\n";
  return 0;
}