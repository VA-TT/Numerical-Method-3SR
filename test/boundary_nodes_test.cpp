#include "../library/Mesh.h"
#include <cassert>
#include <iostream>

static bool contains(const DynamicVector<Index> &ids, Index v) {
  for (Index i = 0; i < ids.size(); ++i) {
    if (ids[i] == v)
      return true;
  }
  return false;
}

int main() {
  // 3x3 nodes on 1x1 domain:
  // y ^
  //   |  6--7--8
  //   |  3--4--5
  //   |  0--1--2
  //   +---------> x
  Mesh2D<double> mesh(1.0, 1.0, 3, 3);

  const auto bottom = mesh.bottomNodes();
  const auto top = mesh.topNodes();
  const auto left = mesh.leftNodes();
  const auto right = mesh.rightNodes();

  // Bottom: 0,1,2
  assert(contains(bottom, 0));
  assert(contains(bottom, 1));
  assert(contains(bottom, 2));

  // Top: 6,7,8
  assert(contains(top, 6));
  assert(contains(top, 7));
  assert(contains(top, 8));

  // Left: 0,3,6
  assert(contains(left, 0));
  assert(contains(left, 3));
  assert(contains(left, 6));

  // Right: 2,5,8
  assert(contains(right, 2));
  assert(contains(right, 5));
  assert(contains(right, 8));

  std::cout << "OK\n";
  return 0;
}
