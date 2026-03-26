#include "../library/Vector.h"
#include <cassert>
#include <iostream>

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

  std::cout << "All static vectors tests passed!\n";
  return 0;
}