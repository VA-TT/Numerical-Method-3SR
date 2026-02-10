#include "../library/MPM.h"
#include "../library/clock.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

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

  std::ofstream hist("mpm1D_history.txt");
  if (!hist.is_open()) {
    std::cerr << "Error: Cannot open mpm1D_history.txt for writing\n";
    return 1;
  }
  hist << "# time\t x_num\t v_num\t x_ana\t v_ana\n";

  // Debug export for step-by-step comparison with the Python reference
  std::ofstream dbg("mpm1D_debug.txt");
  if (!dbg.is_open()) {
    std::cerr << "Error: Cannot open mpm1D_debug.txt for writing\n";
    return 1;
  }
  dbg << "# step\t time\t x_p_n\t v_p_n\t N0\t N1\t mass_n0\t mass_n1\t "
         "mv_n0\t mv_n1\t f_total0\t f_total1\t mv_n0_after\t mv_n1_after\t "
         "a_grid\t v_pic\t x_p_np1\t v_p_np1\n";

  // Initial state at t = 0 (before any update)
  {
    const double t_state = 0.0;
    const double x_num = beam.getMPposition(0);
    const double v_num = beam.getMPvelocity(0);
    const double x_ana = analytic_x(t_state);
    const double v_ana = analytic_v(t_state);
    hist << std::fixed << std::setprecision(10) << t_state << '\t' << x_num
         << '\t' << v_num << '\t' << x_ana << '\t' << v_ana << '\n';
  }

  double time = 0.0;
  for (Index step{0}; step < beam.getNumSteps(); ++step) {
    beam.setupMP();
    beam.p2n();
    beam.nodalEquilibrium();
    // Python-like BC: clamp left node momentum and total force
    beam.applyNodalMomentumConstraint(0, 0.0);
    beam.applyNodalForceConstraint(0, 0.0);
    beam.applyNodalVeloConstraint(0, 0.0);
    // Capture debug quantities at the same point in the step as the Python
    // code:
    // - N, mass_n, mv_n, f_total computed from (x_p^n, v_p^n)
    // - then mv_n is updated by f_total*dt before updating (x_p, v_p)
    const double t_n = time;
    const double x_p_n = beam.getMPposition(0);
    const double v_p_n = beam.getMPvelocity(0);

    const double N0 = 1.0 - std::abs(x_p_n - 0.0) / L;
    const double N1 = 1.0 - std::abs(x_p_n - L) / L;

    const double mass_n0 = beam.getNodalMass(0);
    const double mass_n1 = beam.getNodalMass(1);
    const double mv_n0 = beam.getNodalMomentum(0);
    const double mv_n1 = beam.getNodalMomentum(1);
    const double f_total0 = beam.getNodalTotalForce(0);
    const double f_total1 = beam.getNodalTotalForce(1);

    // Apply the same nodal momentum update as in MPM (and as in the Python
    // snippet) Node 0 is momentum-constrained in this test.
    const double mv_n0_after = mv_n0;
    const double mv_n1_after = mv_n1 + f_total1 * dt;

    const double a_grid =
        (std::abs(mass_n0) < 1e-14 ? 0.0 : (N0 * f_total0 / mass_n0)) +
        (std::abs(mass_n1) < 1e-14 ? 0.0 : (N1 * f_total1 / mass_n1));

    const double v_pic =
        (std::abs(mass_n0) < 1e-14 ? 0.0 : (N0 * mv_n0_after / mass_n0)) +
        (std::abs(mass_n1) < 1e-14 ? 0.0 : (N1 * mv_n1_after / mass_n1));

    beam.n2p();

    const double x_p_np1 = beam.getMPposition(0);
    const double v_p_np1 = beam.getMPvelocity(0);

    dbg << std::fixed << std::setprecision(10) << step << '\t' << t_n << '\t'
        << x_p_n << '\t' << v_p_n << '\t' << N0 << '\t' << N1 << '\t' << mass_n0
        << '\t' << mass_n1 << '\t' << mv_n0 << '\t' << mv_n1 << '\t' << f_total0
        << '\t' << f_total1 << '\t' << mv_n0_after << '\t' << mv_n1_after
        << '\t' << a_grid << '\t' << v_pic << '\t' << x_p_np1 << '\t' << v_p_np1
        << '\n';

    beam.resetMesh();
    // beam.exportResult();

    // State is advanced by dt after n2p()
    const double t_state = time + dt;
    const double x_num = beam.getMPposition(0);
    const double v_num = beam.getMPvelocity(0);
    const double x_ana = analytic_x(t_state);
    const double v_ana = analytic_v(t_state);

    hist << std::fixed << std::setprecision(10) << t_state << '\t' << x_num
         << '\t' << v_num << '\t' << x_ana << '\t' << v_ana << '\n';
    if (step < 5 || step == beam.getNumSteps() - 1) {
      std::cout << std::fixed << std::setprecision(6);
      std::cout << "Step " << step << " t=" << t_state << " | x_num=" << x_num
                << " v_num=" << v_num << " | x_ana=" << x_ana
                << " v_ana=" << v_ana << '\n';
    }

    time += dt;
  }

  hist.close();
  dbg.close();
  std::cout << "\nWrote time history to: mpm1D_history.txt\n";
  std::cout << "Wrote debug history to: mpm1D_debug.txt\n";

  std::cout << "\nComputation time: " << t.elapsed() << " seconds\n";
  return 0;
}
