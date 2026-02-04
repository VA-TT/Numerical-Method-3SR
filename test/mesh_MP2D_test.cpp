#include "../library/Mesh.h"
#include <iomanip>
#include <iostream>

int main() {
  std::cout << std::fixed << std::setprecision(4);

  // Test 1: 2x2 grid (1 element), 2x2 MPs per element
  std::cout << "=== Test 1: 1.0x1.0 domain, 2x2 nodes (1 element), 2x2 "
               "MPs/element ==="
            << '\n';
  Mesh2D<double> mesh1(1.0, 1.0, 2, 2, 2);
  mesh1.print();

  std::cout << "\nMaterial Points coordinates:\n";
  std::cout << "Total MPs: " << mesh1.getNumMPs() << '\n';
  const auto &MPcoords = mesh1.getMPCoords();
  for (Index p = 0; p < MPcoords.size(); ++p) {
    Index elemID = mesh1.findCageID(MPcoords[p].first, MPcoords[p].second);
    std::cout << "  MP " << p << ": (" << MPcoords[p].first << ", "
              << MPcoords[p].second << ") in element " << elemID << '\n';
  }

  // Test 2: 3x3 grid (4 elements), 3x3 MPs per element
  std::cout << "\n=== Test 2: 2.0x1.5 domain, 3x3 nodes (4 elements), 3x3 "
               "MPs/element ==="
            << '\n';
  Mesh2D<double> mesh2(2.0, 1.5, 3, 3, 3);

  std::cout << "Material Points coordinates (first 9 only):\n";
  std::cout << "Total MPs: " << mesh2.getNumMPs() << '\n';
  const auto &MPcoords2 = mesh2.getMPCoords();
  for (Index p = 0; p < std::min(Index(9), Index(MPcoords2.size())); ++p) {
    Index elemID = mesh2.findCageID(MPcoords2[p].first, MPcoords2[p].second);
    std::cout << "  MP " << p << ": (" << MPcoords2[p].first << ", "
              << MPcoords2[p].second << ") in element " << elemID << '\n';
  }

  // Verify all MPs are inside elements
  std::cout << "\n=== Verification ===" << '\n';
  bool allInside = true;
  for (Index p = 0; p < MPcoords2.size(); ++p) {
    Index elemID = mesh2.findCageID(MPcoords2[p].first, MPcoords2[p].second);
    if (elemID == -1) {
      std::cout << "ERROR: MP " << p << " at (" << MPcoords2[p].first << ", "
                << MPcoords2[p].second << ") is outside domain!" << '\n';
      allInside = false;
    }
  }
  if (allInside) {
    std::cout << "SUCCESS: All " << mesh2.getNumMPs()
              << " MPs are inside their elements!" << '\n';
  }

  return 0;
}
