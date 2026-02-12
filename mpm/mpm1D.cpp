#include "../library/MPM.h"
#include "../library/clock.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

int main() {
  Timer t;

  // Input
  const double L = 1.0;
  const double v0 = 0.1;
  const double E = 4.0 * constants::pi * constants::pi;
  const double rho = 1.0;
  const double dt = 0.01;
  const double duration = 10.0;
  const double xloc = 0.5;

  // Analytic result
  const double omega = (1.0 / L) * std::sqrt(E / rho);
  auto analytic_v = [&](double time) { return v0 * std::cos(omega * time); };
  auto analytic_x = [&](double time) {
    return xloc * std::exp((v0 / (L * omega)) * std::sin(omega * time));
  };

  // Set up MPM1D grid with 2 nodes and 1 MP
  using Beam = MPM1D<double, 2, 1>;
  Beam beam(E, rho, L, v0, dt, duration, xloc);
  // (dt, duration, rho, length, xloc, v0, a0) // v0 = a0 = 0 by default
  beam.setG(0.0);             // G = constants::gravity if gravity is considered
  beam.setE(E);               // Module Young
  beam.setComportmentLaw({}); // Tangential Operator

  // Output files
  std::ofstream hist("mpm1D_history.txt");
  hist << "# time\tx_num\tv_num\tx_ana\tv_ana\n";

  std::ofstream stress_strain("mpm1D_stress_strain.txt");
  stress_strain << "# time\tstress\tstrain\n";

  // Set boundary conditions ONCE
  beam.setNodalVeloConstraint(0, 0.0);
  beam.setNodalAccConstraint(0, 0.0);
  beam.setNodalMomentumConstraint(0, 0.0);
  beam.setNodalForceConstraint(0, 0.0);

  // Find MP closest to xloc for tracking
  Index tracked_mp = 0;
  double min_dist = std::abs(beam.getMPposition(0) - xloc);
  for (Index p = 1; p < beam.getNumMps(); ++p) {
    double dist = std::abs(beam.getMPposition(p) - xloc);
    if (dist < min_dist) {
      min_dist = dist;
      tracked_mp = p;
    }
  }

  std::cout << "Tracking MP[" << tracked_mp
            << "] at x=" << beam.getMPposition(tracked_mp) << " (xloc=" << xloc
            << ")\n\n";

  // Initial state (t=0) to the output file
  {
    const double x_num = beam.getMPposition(tracked_mp);
    const double v_num = beam.getMPvelocity(tracked_mp);
    const double x_ana = analytic_x(0.0);
    const double v_ana = analytic_v(0.0);
    hist << std::fixed << std::setprecision(10) << 0.0 << '\t' << x_num << '\t'
         << v_num << '\t' << x_ana << '\t' << v_ana << '\n';

    stress_strain << std::fixed << std::setprecision(10) << 0.0 << '\t'
                  << beam.getMPstress(tracked_mp) << '\t'
                  << beam.getMPstrain(tracked_mp) << '\n';
  }

  // Time integration loop
  for (Index step = 0; step < beam.getNumSteps(); ++step) {
    double time = (step + 1) * beam.getTimeStep();

    // MPM algorithm
    beam.setupMP(); // Update MP elements and activate nodes
    beam.p2n();
    beam.nodalEquilibrium();
    beam.n2p();
    beam.resetMesh();

    // Record results
    const double x_num = beam.getMPposition(tracked_mp);
    const double v_num = beam.getMPvelocity(tracked_mp);
    const double x_ana = analytic_x(time);
    const double v_ana = analytic_v(time);

    hist << std::fixed << std::setprecision(10) << time << '\t' << x_num << '\t'
         << v_num << '\t' << x_ana << '\t' << v_ana << '\n';

    stress_strain << std::fixed << std::setprecision(10) << time << '\t'
                  << beam.getMPstress(tracked_mp) << '\t'
                  << beam.getMPstrain(tracked_mp) << '\n';

    if (step < 5 || step == beam.getNumSteps() - 1) {
      std::cout << std::fixed << std::setprecision(6);
      std::cout << "t=" << time << " | Num: x=" << x_num << " v=" << v_num
                << " | Ana: x=" << x_ana << " v=" << v_ana
                << " | Err: Δx=" << std::abs(x_num - x_ana)
                << " Δv=" << std::abs(v_num - v_ana) << '\n';
    }
  }

  hist.close();
  stress_strain.close();
  std::cout << "\n✓ Results saved to: mpm1D_history.txt\n";
  std::cout << "✓ Stress-strain data saved to: mpm1D_stress_strain.txt\n";
  std::cout << "✓ Computation time: " << t.elapsed() << " s\n";
  return 0;
}
