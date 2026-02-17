#include "../library/MPM.h"
#include "../library/clock.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

int main() {
  Timer t;
  // Material properties
  const double E = 100000.0;
  const double v = 0.3;
  const double rho = 3600.0;
  const double mu = 0.385;
  const double phi = 30;
  const double c = 1;
  const double K0 = 0.5;

  // Analysis setting
  const double dt = 0.0001;
  const double duration = 1.0;
  const double interval = 200;
  const double t = 0;
  const double g = 4.8;

  const double L = 1.0;
  const double nx = 20;
  const double ny = 20;

  // Computational parameters
  const double c = std::sqrt(E / rho);
  const Index nelements = 13;
  const double dx = L / nelements;
  const double dt_crit = dx / c;
  const double dt = 0.1 * dt_crit;

  // Particles per cell
  const Index ppc = 2;
  const Index nparticles = nelements * ppc; // 26 particles total

  // Position to track (particle closest to center, left of center)
  const Index pmid = 12;       // MP[12] at x=12.1795 (closest to center)
  const double xloc = 0.5 * L; // x = 12.5 (center of domain)

  // Wave parameters
  const double beta1 = constants::pi / (2.0 * L);
  const double omega1 = beta1 * c;

  // Analytic result for standing wave
  auto analytic_v = [&](double x, double time) {
    return v0 * std::sin(beta1 * x) * std::cos(omega1 * time);
  };
  auto analytic_x = [&](double x0, double time) {
    return x0 + (v0 / omega1) * std::sin(beta1 * x0) * std::sin(omega1 * time);
  };

  // Set up MPM1D grid: 14 nodes, 13 elements, 2 MPs per element
  using Beam = MPM1D<double, 14, 2>;
  Beam beam(E, rho, L, v0, dt, duration, xloc);
  beam.setG(0.0);
  beam.setE(E);
  beam.setComportmentLaw({});

  // Output files
  std::ofstream hist("mpm1Dc_history.txt");
  hist << "# time\tx_num\tv_num\tx_ana\tv_ana\n";

  std::ofstream stress_strain("mpm1Dc_stress_strain.txt");
  stress_strain << "# time\tstress\tstrain\n";

  // Set non-uniform initial velocity: v(x) = v0 * sin(beta1 * x)
  for (Index p = 0; p < beam.getNumMps(); ++p) {
    double x_p = beam.getMPposition(p);
    beam.setMPvelocity(p, v0 * std::sin(beta1 * x_p));
  }

  // Set boundary conditions
  beam.setNodalVeloConstraint(0, 0.0);
  beam.setNodalMomentumConstraint(0, 0.0);
  beam.setNodalForceConstraint(0, 0.0);

  // Track middle particle
  const Index tracked_mp = pmid;
  const double x0_tracked = beam.getMPposition(tracked_mp);

  std::cout << "=== C++ MPM Standing Wave (mpm1Dc - " << ppc
            << " particles/element) ===\n";
  std::cout << "L=" << L << ", E=" << E << ", rho=" << rho << ", v0=" << v0
            << "\n";
  std::cout << "c=" << c << ", beta1=" << beta1 << ", omega1=" << omega1
            << "\n";
  std::cout << "dt=" << dt << ", nsteps=" << beam.getNumSteps() << "\n";
  std::cout << "Total particles: " << beam.getNumMps() << " (" << ppc
            << " per element)\n";
  std::cout << "Tracking MP[" << tracked_mp << "] at x=" << x0_tracked
            << " (domain center=" << xloc << ")\n";
  std::cout << "Initial velocity: v=" << beam.getMPvelocity(tracked_mp)
            << " (expected: " << v0 * std::sin(beta1 * x0_tracked) << ")\n\n";

  // Initial state (t=0)
  {
    const double x_num = beam.getMPposition(tracked_mp);
    const double v_num = beam.getMPvelocity(tracked_mp);
    const double x_ana = x0_tracked;                        // Initial position
    const double v_ana = v0 * std::sin(beta1 * x0_tracked); // Initial velocity
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
    beam.setupMP();
    beam.p2n();
    beam.nodalEquilibrium();
    beam.n2p();
    beam.resetMesh();

    // Record results
    const double x_num = beam.getMPposition(tracked_mp);
    const double v_num = beam.getMPvelocity(tracked_mp);
    const double x_ana = analytic_x(x0_tracked, time);
    const double v_ana = analytic_v(x0_tracked, time);

    hist << std::fixed << std::setprecision(10) << time << '\t' << x_num << '\t'
         << v_num << '\t' << x_ana << '\t' << v_ana << '\n';

    stress_strain << std::fixed << std::setprecision(10) << time << '\t'
                  << beam.getMPstress(tracked_mp) << '\t'
                  << beam.getMPstrain(tracked_mp) << '\n';

    if (step < 5 || step % 1000 == 0 || step == beam.getNumSteps() - 1) {
      std::cout << std::fixed << std::setprecision(6);
      std::cout << "t=" << time << " | Num: x=" << x_num << " v=" << v_num
                << " | Ana: x=" << x_ana << " v=" << v_ana
                << " | Err: Δx=" << std::abs(x_num - x_ana)
                << " Δv=" << std::abs(v_num - v_ana) << '\n';
    }
  }

  hist.close();
  stress_strain.close();
  std::cout << "\n✓ Results saved to: mpm1Dc_history.txt\n";
  std::cout << "✓ Stress-strain data saved to: mpm1Dc_stress_strain.txt\n";
  std::cout << "✓ Computation time: " << t.elapsed() << " s\n";
  return 0;
}
