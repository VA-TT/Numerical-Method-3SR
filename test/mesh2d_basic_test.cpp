// Test đơn giản cho Mesh2D
#include "../versionSoA/Mesh.h"
#include <cassert>
#include <iostream>

int main() {
  using T = double;
  T length = 4.0, height = 3.0;
  int nx = 5, ny = 4;
  int nMPperEle = 2;
  Mesh2D<T> mesh(length, height, nx, ny, nMPperEle);

  // Kiểm tra số node, element
  assert(mesh.getNumNodes() == nx * ny);
  assert(mesh.getNumElements() == (nx - 1) * (ny - 1));

  // Kiểm tra số MP
  int expectedMP = (nx - 1) * (ny - 1) * nMPperEle * nMPperEle;
  assert(mesh.getNumMPs() == expectedMP);

  // Kiểm tra vị trí MP nằm trong domain
  for (int i = 0; i < mesh.getNumMPs(); ++i) {
    auto pos = mesh.getMPpos(i);
    assert(pos.first >= 0.0 && pos.first <= length);
    assert(pos.second >= 0.0 && pos.second <= height);
    int elem = mesh.getMPelementID(i);
    assert(elem >= 0 && elem < mesh.getNumElements());
  }
  std::cout << "Mesh2D basic test PASSED!\n";
  return 0;
}
