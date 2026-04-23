// Test đơn giản cho Mesh1D
#include "../versionSoA/Mesh.h"
#include <cassert>
#include <iostream>

int main() {
  using T = double;
  T length = 10.0;
  int nNodes = 6;
  int nMPperEle = 2;
  Mesh1D<T> mesh(length, nNodes, nMPperEle);

  // Kiểm tra số node, element, MP
  assert(mesh.getNumNodes() == nNodes);
  assert(mesh.getNumElements() == nNodes - 1);
  assert(mesh.getNumMPs() == (nNodes - 1) * nMPperEle);

  // Kiểm tra vị trí node đầu/cuối
  auto nodes = mesh.nodeCoords();
  assert(nodes[0] == 0.0);
  assert(nodes[nNodes - 1] == length);

  // Kiểm tra vị trí MP nằm trong domain
  for (int i = 0; i < mesh.getNumMPs(); ++i) {
    auto x = mesh.getMPpos(i);
    assert(x >= 0.0 && x <= length);
    int elem = mesh.getMPelementID(i);
    assert(elem >= 0 && elem < mesh.getNumElements());
  }
  std::cout << "Mesh1D basic test PASSED!\n";
  return 0;
}
