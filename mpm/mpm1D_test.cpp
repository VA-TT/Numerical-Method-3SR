#include "../library/MPM.h"
#include "../library/clock.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

int main() {
  Timer t;

  // Parameters matching Python code
  // Match Sample.py
  const double L = 1.0;
  const double v0 = 0.1;
  const double E = 4.0 * constants::pi * constants::pi;
  const double rho = 1.0;
  const double dt = 0.01;
  const double duration = 10.0;
  const double xloc = 0.5; // same as Sample.py (only used for optional compare)

  using Beam = MPM1D<double, 2, 1>; // 2 nodes, 1 MP/element
  Beam beam(E, rho, L, v0, dt, duration, xloc);
  // (dt, duration, rho, length, xloc, v0, a0) // v0 = a0 = 0 by default
  beam.setG(0.0); // G = constants::gravity if gravity is considered
  beam.setE(E);   // Module Young
  beam.setComportmentLaw({});

  const double omega = (1.0 / L) * std::sqrt(E / rho);
  auto analytic_v = [&](double time) { return v0 * std::cos(omega * time); };
  auto analytic_x = [&](double time) {
    return xloc * std::exp((v0 / (L * omega)) * std::sin(omega * time));
  };

  double time = 0.0;
  for (Index step{0}; step < beam.getNumSteps(); ++step) {
    beam.setupMP();
    beam.p2n();
    beam.nodalEquilibrium();
    // Python-like BC: clamp left node momentum and total force
    beam.applyNodalMomentumConstraint(0, 0.0);
    beam.applyNodalForceConstraint(0, 0.0);
    beam.applyNodalVeloConstraint(0, 0.0);
    beam.n2p();
    beam.resetMesh();
    // beam.exportResult();

    // State is advanced by dt after n2p()
    const double t_state = time + dt;
    if (step < 5 || step == beam.getNumSteps() - 1) {
      std::cout << std::fixed << std::setprecision(6);
      std::cout << "Step " << step << " t=" << t_state
                << " | x_num=" << beam.getMPposition(0)
                << " v_num=" << beam.getMPvelocity(0)
                << " | x_ana=" << analytic_x(t_state)
                << " v_ana=" << analytic_v(t_state) << '\n';
    }

    time += dt;
  }

  std::cout << "\nComputation time: " << t.elapsed() << " seconds\n";
  return 0;
}
