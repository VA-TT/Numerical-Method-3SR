#include "../library/Mesh.h"
#include <iostream>

int main() {
  // Domain [0,1]x[0,1], simple grid
  // MPs: generate in [0,1]x[0,1] with spacing 0.5 -> centers at 0.25 and 0.75
  Mesh2D<double> mesh(1.0, 1.0, 2, 2, 0.0, 0.0, 1.0, 1.0, 0.5);

  // Default convention: R = MP_size/2
  const auto mask = mesh.getMPBoundaryContactMask();

  std::cout << "MP count: " << mesh.getNumMPs() << "\n";
  for (Index p = 0; p < mesh.getNumMPs(); ++p) {
    const auto [x, y] = mesh.getMPpos(p);
    std::cout << "MP " << p << " at (" << x << "," << y
              << ") mask=" << static_cast<int>(mask[p]) << "\n";
  }

  // Expected: all MPs at 0.25 or 0.75 with R=0.25 touch at least one boundary.
  // (0.25-0.25==0) => left/bottom; (0.75+0.25==1) => right/top.
  bool allTouch = true;
  for (Index p = 0; p < mesh.getNumMPs(); ++p) {
    if (!mesh.isMPContactBound(p))
      allTouch = false;
  }

  if (!allTouch) {
    std::cout << "FAILED: some MPs do not touch boundary but expected to.\n";
    return 1;
  }

  std::cout << "OK\n";
  return 0;
}
