#include "../library/mesh2.h"
#include <cassert>
#include <cmath>
#include <iostream>

int main() {
  std::cout << "Testing Particle-based Mesh refactor...\n";

  // ===== Test Mesh1D with Particle1D =====
  std::cout << "\n[Test 1] Mesh1D with Particle1D<double>\n";
  {
    Mesh1D<double> mesh(10.0, 5, 2); // length=10, nodes=5, MPs=2 per element
    assert(mesh.getNumMPs() == 6 &&
           "Mesh1D should have 6 particles (2 per element x 3 elements)");

    // Access particles
    const auto &particles = mesh.getMPs();
    assert(particles.size() == 6 && "getMPs should return vector of size 6");

    // Get individual particle
    auto p0 = mesh.getMP(0);
    assert(p0.pos > 0 && "First particle should have pos > 0");

    // Get just position (backward compatible)
    double coord = mesh.getMPpos(0);
    assert(coord == p0.pos && "getMPpos should match particle.pos");

    std::cout << "  ✓ Particle1D access works correctly\n";
  }

  // ===== Test Mesh2D with Particle2D =====
  std::cout << "\n[Test 2] Mesh2D with Particle2D<double> - square grid\n";
  {
    Mesh2D<double> mesh(10.0, 10.0, 3, 3); // 3x3 grid
    assert(mesh.getNumMPs() == 0 && "Initial mesh should have 0 MPs");

    // Generate square grid of MPs
    mesh.generateSquareGridMP(1.0, 1.0, 9.0, 9.0, 2.0);
    Index nMPs = mesh.getNumMPs();
    assert(nMPs > 0 && "Should have generated MPs");
    std::cout << "  Generated " << nMPs << " particles\n";

    // Access particles via new API
    const auto &particles = mesh.getParticles();
    assert(particles.size() == nMPs &&
           "getParticles() should return all particles");

    // Check particle attributes
    auto p0 = mesh.getMP(0);
    assert(p0.pos.x() >= 1.0 && p0.pos.x() <= 9.0 &&
           "MP x should be in domain");
    assert(p0.pos.y() >= 1.0 && p0.pos.y() <= 9.0 &&
           "MP y should be in domain");
    std::cout << "  Particle 0: pos = (" << p0.pos.x() << ", " << p0.pos.y()
              << "), volume = " << p0.volume << "\n";

    // Test backward-compatible getMPCoords() returning pairs
    const auto &mp_pairs = mesh.getMPCoords();
    assert(mp_pairs.size() == nMPs &&
           "Backward-compat getMPCoords should work");

    std::cout
        << "  ✓ Particle2D access and square grid generation work correctly\n";
  }

  // ===== Test Mesh2D circular MP generation =====
  std::cout << "\n[Test 3] Mesh2D with Particle2D<double> - circular region\n";
  {
    Mesh2D<double> mesh(10.0, 10.0, 3, 3);

    // Generate circular region of MPs
    mesh.generateCircleMPgrid({5.0, 5.0}, 3.0, 0.5);
    Index nMPs = mesh.getNumMPs();
    assert(nMPs > 0 && "Should have generated circular MPs");
    std::cout << "  Generated " << nMPs << " particles in circle\n";

    // Verify all particles are within circle
    const auto &particles = mesh.getParticles();
    for (Index p = 0; p < nMPs; ++p) {
      double dx = particles[p].pos.x() - 5.0;
      double dy = particles[p].pos.y() - 5.0;
      double dist_sq = dx * dx + dy * dy;
      assert(dist_sq <= 3.0 * 3.0 + 1e-6 && "Particle should be in circle");
    }

    std::cout << "  ✓ Circular MP generation maintains geometric constraint\n";
  }

  // ===== Test MP modification =====
  std::cout << "\n[Test 4] Modifying particle positions and properties\n";
  {
    Mesh2D<double> mesh(10.0, 10.0, 3, 3);
    mesh.generateSquareGridMP(2.0, 2.0, 8.0, 8.0, 2.0);
    Index nMPs = mesh.getNumMPs();

    // Modify particle position
    mesh.setMPCoord(0, 1.5, 1.5);
    auto p0 = mesh.getMP(0);
    assert(std::abs(p0.pos.x() - 1.5) < 1e-10 && "setMPCoord should update x");
    assert(std::abs(p0.pos.y() - 1.5) < 1e-10 && "setMPCoord should update y");

    // Modify all particles
    DynamicVector<double> new_x(nMPs), new_y(nMPs);
    for (Index p = 0; p < nMPs; ++p) {
      new_x[p] = 5.0 + p * 0.1;
      new_y[p] = 5.0 - p * 0.1;
    }
    mesh.setMPCoords(new_x, new_y);

    auto p_check = mesh.getMP(0);
    assert(std::abs(p_check.pos.x() - new_x[0]) < 1e-10 &&
           "setMPCoords should update all");

    std::cout << "  ✓ Particle position modification works\n";
  }

  // ===== Test MP reset functionality =====
  std::cout << "\n[Test 5] Reset MPs to initial configuration\n";
  {
    Mesh2D<double> mesh(10.0, 10.0, 3, 3);
    mesh.generateSquareGridMP(2.0, 2.0, 8.0, 8.0, 1.0);
    Index nMPs = mesh.getNumMPs();

    // Store initial state
    auto p0_initial = mesh.getMP(0);
    double x_initial = p0_initial.pos.x();
    double y_initial = p0_initial.pos.y();

    // Modify
    mesh.setMPCoord(0, 9.0, 9.0);
    assert(std::abs(mesh.getMP(0).pos.x() - 9.0) < 1e-10 &&
           "Should be modified");

    // Reset
    mesh.resetMPs();
    auto p0_reset = mesh.getMP(0);
    assert(std::abs(p0_reset.pos.x() - x_initial) < 1e-10 &&
           "resetMPs should restore initial pos");
    assert(std::abs(p0_reset.pos.y() - y_initial) < 1e-10 &&
           "resetMPs should restore initial pos");

    std::cout << "  ✓ MP reset functionality works\n";
  }

  // ===== Test activate nodes/elements =====
  std::cout << "\n[Test 6] Activate nodes and elements containing MPs\n";
  {
    Mesh2D<double> mesh(10.0, 10.0, 5, 5);
    mesh.generateSquareGridMP(3.0, 3.0, 7.0, 7.0, 1.0);

    // Activate nodes
    mesh.activateNodes();
    const auto &active_nodes = mesh.getActiveNodes();
    assert(active_nodes.size() == mesh.getNumNodes() &&
           "Should have node activity data");

    // Count active nodes
    Index count_active = 0;
    for (Index i = 0; i < mesh.getNumNodes(); ++i) {
      if (mesh.isActiveNode(i))
        count_active++;
    }
    assert(count_active > 0 && "Should have some active nodes");
    std::cout << "  Active nodes: " << count_active << " out of "
              << mesh.getNumNodes() << "\n";

    // Activate elements
    mesh.activateElements();
    Index count_active_ele = 0;
    for (Index e = 0; e < mesh.getNumElements(); ++e) {
      if (mesh.isActiveElement(e))
        count_active_ele++;
    }
    assert(count_active_ele > 0 && "Should have some active elements");
    std::cout << "  Active elements: " << count_active_ele << " out of "
              << mesh.getNumElements() << "\n";

    std::cout << "  ✓ Activation of nodes and elements works\n";
  }

  std::cout << "\n✅ All Particle refactor tests passed!\n";
  return 0;
}
