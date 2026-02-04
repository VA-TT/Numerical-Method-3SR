#include "../library/Mesh.h"
#include <iomanip>
#include <iostream>

int main() {
  std::cout << std::fixed << std::setprecision(4);

  // Test 1: 3 elements, 2 MPs per element
  std::cout << "=== Test 1: L=1.0, 4 nodes (3 elements), 2 MPs/element ==="
            << '\n';
  Mesh1D<double> mesh1(1.0, 4, 2);
  mesh1.print();

  std::cout << "\nMaterial Points coordinates:\n";
  std::cout << "Total MPs: " << mesh1.getNumMPs() << '\n';
  const auto &MPcoords = mesh1.getMPCoords();
  for (Index p = 0; p < MPcoords.size(); ++p) {
    Index elemID = mesh1.findCageID(MPcoords[p]);
    std::cout << "  MP " << p << ": x = " << MPcoords[p] << " (in element "
              << elemID << ")" << '\n';
  }

  // Test 2: 5 elements, 3 MPs per element
  std::cout << "\n=== Test 2: L=2.0, 6 nodes (5 elements), 3 MPs/element ==="
            << '\n';
  Mesh1D<double> mesh2(2.0, 6, 3);

  std::cout << "Material Points coordinates:\n";
  std::cout << "Total MPs: " << mesh2.getNumMPs() << '\n';
  const auto &MPcoords2 = mesh2.getMPCoords();
  for (Index p = 0; p < MPcoords2.size(); ++p) {
    Index elemID = mesh2.findCageID(MPcoords2[p]);
    std::cout << "  MP " << p << ": x = " << MPcoords2[p] << " (in element "
              << elemID << ")" << '\n';
  }

  // Verify MPs are inside elements
  std::cout << "\n=== Verification ===" << '\n';
  bool allInside = true;
  for (Index p = 0; p < MPcoords2.size(); ++p) {
    Index elemID = mesh2.findCageID(MPcoords2[p]);
    if (elemID == -1) {
      std::cout << "ERROR: MP " << p << " is outside domain!" << '\n';
      allInside = false;
    }
  }
  if (allInside) {
    std::cout << "SUCCESS: All MPs are inside their elements!" << '\n';
  }

  return 0;
}
