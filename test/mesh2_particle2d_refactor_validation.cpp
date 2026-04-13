#include "../library/mesh2.h"
#include <cassert>
#include <cmath>
#include <iostream>

int main() {
  using T = double;

  // 1) Particle2D creation and pos initialization
  {
    Mesh2D<T> mesh(1.0, 1.0, 3, 3, 0.0, 0.0, 1.0, 1.0, 0.25);
    assert(mesh.getNumMPs() > 0);

    const auto &particles = mesh.getParticles();
    assert(particles.size() > 0);
    const auto p0 = mesh.getMP(0);
    assert(std::isfinite(static_cast<double>(p0.pos.x())));
    assert(std::isfinite(static_cast<double>(p0.pos.y())));
  }

  // 2) MP generation methods
  {
    Mesh2D<T> meshPerElement(2.0, 2.0, 5, 5, 2);
    assert(meshPerElement.getNumMPs() == meshPerElement.getNumElements() * 4);

    Mesh2D<T> meshSquare(2.0, 2.0, 5, 5, 0.2, 0.2, 1.8, 1.8, 0.2);
    assert(meshSquare.getNumMPs() > 0);

    Mesh2D<T> meshCircle(2.0, 2.0, 5, 5, std::pair<T, T>{1.0, 1.0}, 0.5, 0.1);
    assert(meshCircle.getNumMPs() > 0);
    const auto &parts = meshCircle.getParticles();
    for (Index p = 0; p < meshCircle.getNumMPs(); ++p) {
      const T dx = parts[p].pos.x() - 1.0;
      const T dy = parts[p].pos.y() - 1.0;
      assert(dx * dx + dy * dy <= 0.25 + 1e-12);
    }
  }

  // 3) resetMPs copies from m_MPs_initial
  {
    Mesh2D<T> mesh(2.0, 2.0, 5, 5, 0.2, 0.2, 1.8, 1.8, 0.2);
    const auto before = mesh.getMP(0);

    mesh.updateMPCoord(0, 1.7, 1.7);
    const auto moved = mesh.getMP(0);
    assert(std::abs(static_cast<double>(moved.pos.x() - before.pos.x())) >
           1e-10);

    mesh.resetMPs();
    const auto after = mesh.getMP(0);
    assert(std::abs(static_cast<double>(after.pos.x() - before.pos.x())) <
           1e-12);
    assert(std::abs(static_cast<double>(after.pos.y() - before.pos.y())) <
           1e-12);
  }

  // 4) Particle properties persist when setters update only pos
  {
    Mesh2D<T> mesh(2.0, 2.0, 5, 5, 0.2, 0.2, 1.8, 1.8, 0.2);

    auto &particles = mesh.getParticles();
    particles[0].mass = 42.0;
    particles[0].volume = 2.5;
    particles[0].vel.x() = 3.0;
    particles[0].vel.y() = -1.0;

    mesh.setMPCoord(0, 1.0, 1.0);
    const auto p = mesh.getMP(0);
    assert(std::abs(static_cast<double>(p.pos.x() - 1.0)) < 1e-12);
    assert(std::abs(static_cast<double>(p.pos.y() - 1.0)) < 1e-12);
    // Existing non-position attributes should remain untouched.
    assert(std::abs(static_cast<double>(p.mass - 42.0)) < 1e-12);
    assert(std::abs(static_cast<double>(p.volume - 2.5)) < 1e-12);
    assert(std::abs(static_cast<double>(p.vel.x() - 3.0)) < 1e-12);
    assert(std::abs(static_cast<double>(p.vel.y() + 1.0)) < 1e-12);
  }

  // Backward compatibility helper
  {
    Mesh2D<T> mesh(2.0, 2.0, 5, 5, 0.2, 0.2, 1.8, 1.8, 0.2);
    const auto &coords = mesh.getMPCoords();
    assert(static_cast<Index>(coords.size()) == mesh.getNumMPs());
  }

  std::cout << "mesh2 Particle2D refactor validation passed\n";
  return 0;
}
