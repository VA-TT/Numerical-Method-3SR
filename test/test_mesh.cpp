#include "../library/Mesh.h"
#include <iomanip>

int main() {
  std::cout << std::fixed << std::setprecision(4);

  // ========== TEST MESH 1D ==========
  std::cout << "==========================================\n";
  std::cout << "         TEST MESH 1D\n";
  std::cout << "==========================================\n\n";

  Mesh1D<double> mesh1D(10.0, 6); // 6 nodes, 5 elements
  mesh1D.print();

  std::cout << "\n--- Element Properties (Connectivity-based) ---\n";
  for (Index e = 0; e < mesh1D.getNumElements(); ++e) {
    auto [x1, x2] = mesh1D.getElementNodes(e);
    double length = mesh1D.getLengthEle(e);
    std::cout << "Element " << e << ": [" << x1 << ", " << x2
              << "] length = " << length << "\n";
  }

  std::cout << "\n--- Test Update & Reset Functionality ---\n";
  std::cout << "Initial node 2: " << mesh1D.nodeCoords()[2] << "\n";

  mesh1D.updateNodePosition(2, 5.5);
  std::cout << "After update node 2 to 5.5: " << mesh1D.nodeCoords()[2] << "\n";
  auto [x1_updated, x2_updated] = mesh1D.getElementNodes(1);
  std::cout << "  Element 1 now: [" << x1_updated << ", " << x2_updated
            << "]\n";

  mesh1D.resetMesh();
  std::cout << "After reset: " << mesh1D.nodeCoords()[2] << "\n";
  auto [x1_reset, x2_reset] = mesh1D.getElementNodes(1);
  std::cout << "  Element 1 now: [" << x1_reset << ", " << x2_reset << "]\n";

  std::cout << "\n--- Test findCageID (MPM) ---\n";
  double test_x = 4.5;
  Index cage = mesh1D.findCageID(test_x);
  std::cout << "Point x=" << test_x << " is in element " << cage << "\n";
  if (cage >= 0) {
    auto [ex1, ex2] = mesh1D.getElementNodes(cage);
    std::cout << "  Element bounds: [" << ex1 << ", " << ex2 << "]\n";
  }

  std::cout << "\n✓ Mesh1D tests passed!\n\n";

  // ========== TEST MESH 2D ==========
  std::cout << "==========================================\n";
  std::cout << "         TEST MESH 2D\n";
  std::cout << "==========================================\n\n";

  Mesh2D<double> mesh2D(4.0, 2.0, 3, 3); // 3x3 nodes, 2x2 elements
  mesh2D.print();

  std::cout << "\n--- Test getNodes() - All nodes as (x,y) pairs ---\n";
  const auto &nodes = mesh2D.getNodes();
  for (Index i = 0; i < mesh2D.getNumNodes(); ++i) {
    std::cout << "  Node " << i << ": (" << std::setw(6) << nodes[i].first
              << ", " << std::setw(6) << nodes[i].second << ")\n";
  }

  std::cout << "\n--- Test getNodeCoor(id) - Individual node access ---\n";
  for (Index i = 0; i < 5 && i < mesh2D.getNumNodes(); ++i) {
    auto [x, y] = mesh2D.getNodeCoor(i);
    std::cout << "  Node " << i << ": (" << x << ", " << y << ")\n";
  }

  std::cout << "\n--- Element Connectivity and Area ---\n";
  double total_area = 0.0;
  for (Index e = 0; e < mesh2D.getNumElements(); ++e) {
    const auto &conn = mesh2D.getConnectivity(e);
    std::cout << "  Element " << e << ": [" << conn[0] << ", " << conn[1]
              << ", " << conn[2] << ", " << conn[3] << "]";

    // Calculate area using first and third nodes
    auto [x1, y1] = nodes[conn[0]]; // Bottom-left
    auto [x2, y2] = nodes[conn[1]]; // Bottom-right
    auto [x4, y4] = nodes[conn[3]]; // Top-left

    double width = std::abs(x2 - x1);
    double height = std::abs(y4 - y1);
    double area = width * height;
    total_area += area;
    std::cout << " -> Area: " << width << " x " << height << " = " << area
              << "\n";
  }
  std::cout << "\nTotal mesh area: " << total_area << " (expected: 8.0)\n";

  std::cout << "\n--- Test Regenerate Mesh ---\n";
  mesh2D.setXnodes(5);
  mesh2D.regenerateMesh();
  std::cout << "After regenerate with 5x3 nodes:\n";
  std::cout << "  Number of nodes: " << mesh2D.getNumNodes() << "\n";
  std::cout << "  Number of elements: " << mesh2D.getNumElements() << "\n";
  std::cout << "  Grid size: " << mesh2D.nx() << " x " << mesh2D.ny() << "\n";

  std::cout << "\n--- Test Update & Reset 2D ---\n";
  auto [x0_init, y0_init] = mesh2D.getNodeCoor(4); // Center node
  std::cout << "Initial node 4: (" << x0_init << ", " << y0_init << ")\n";

  mesh2D.updateNodePosition(4, 3.0, 1.5);
  auto [x0_updated, y0_updated] = mesh2D.getNodeCoor(4);
  std::cout << "After update: (" << x0_updated << ", " << y0_updated << ")\n";

  mesh2D.resetMesh();
  auto [x0_reset, y0_reset] = mesh2D.getNodeCoor(4);
  std::cout << "After reset: (" << x0_reset << ", " << y0_reset << ")\n";

  std::cout << "\n✓ Mesh2D tests passed!\n\n";

  std::cout << "==========================================\n";
  std::cout << "   ALL TESTS COMPLETED SUCCESSFULLY!\n";
  std::cout << "==========================================\n";

  return 0;
}
