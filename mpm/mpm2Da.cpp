#include "../library/MPM.h"
#include "../library/clock.h"
#include "../library/ioDirectory.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

int main() {
  Timer timer;
  // Material properties
  const double E = 100000.0;
  const double nu = 0.3;
  const double rho = 3600.0;
  const double mu = 0.385; // wall-obstacle
  const double phi = 30;   // internal friction andgle (in degree)
  const double cohesion = 1.0;
  const double K0 = 0.5;

  // Analysis setting
  const double duration = 1.0;

  // Grid definition (must be constexpr to be used in MPM2D template args)
  static constexpr double L = 1.0;
  static constexpr double H = 1.0;
  static constexpr Index nx = 20;
  static constexpr Index ny = 20;
  static constexpr double dx = L / static_cast<double>(nx);
  static constexpr double dy = H / static_cast<double>(ny);

  // Particles per cell
  static constexpr Index ppc = 2;
  static constexpr double MP_size = dx / static_cast<double>(ppc);

  // MP domain
  const StaticVector<double, 2> minCorner{0.0, 0.0};
  const StaticVector<double, 2> maxCorner{0.3, 0.3};

  // Computational parameters (CFL-ish)
  const double c_wave = std::sqrt(E / rho);
  const double dt_crit = std::min(dx, dy) / c_wave;
  const double dt = 0.1 * dt_crit;

  const double v0 = 0.0;

  // Set up MPM2D grid: 14 nodes, 13 elements, 2 MPs per element
  using collapse2D = MPM2D<double, L, H, nx, ny, MP_size>;
  collapse2D collumn(rho, E, nu, phi, cohesion, mu, minCorner, maxCorner, dt,
                     duration, v0);
  collumn.setG(-9.81); // G=-9.81
  collumn.setE(E);
  collumn.setComportmentLaw(std::function<double(double)>{});
  collumn.setShape(shapePolicy::cubicBSpline);

  // Output files
  std::ofstream hist("mpm2Da_history.txt");
  hist << "# time\tx\ty\tstress_xx\tstrain_xx\n";

  // Simple boundary condition: clamp bottom active nodes
  for (const Index nodeID : collumn.getMesh().bottomActiveNodes()) {
    collumn.setNodalVeloConstraint(nodeID, 0.0);
  }

  // Track a representative particle
  const Index tracked_mp =
      collumn.getNumMPs() > 0 ? collumn.getNumMPs() / 2 : 0;

  std::cout << "=== C++ MPM2D (mpm2Da) ===\n";
  std::cout << "L=" << L << ", H=" << H << ", E=" << E << ", rho=" << rho
            << ", mu=" << mu << "\n";
  std::cout << "dt=" << dt << ", nsteps=" << collumn.getNumSteps() << "\n";
  std::cout << "Total particles: " << collumn.getNumMPs() << "\n";
  std::cout << "Tracking MP[" << tracked_mp << "]\n\n";

  // VTK output (ParaView): mpm/data2D/particles_000000.vtk and mesh_000000.vtk
  // (relative to where the executable is).
  const Index vtkInterval = 10;
  const std::filesystem::path vtkDir = ioFile::vtkOutputDir("data2D");

  // Initial state (t=0)
  {
    const auto mp = collumn.getMesh().getMP(tracked_mp).pos;
    hist << std::fixed << std::setprecision(10) << 0.0 << '\t' << mp.x() << '\t'
         << mp.y() << '\t' << collumn.getMPstress(tracked_mp) << '\t'
         << collumn.getMPstrain(tracked_mp) << '\n';

    collumn.exportVTKFrame(vtkDir, 0);
  }

  // Time integration loop
  for (Index step = 0; step < collumn.getNumSteps(); ++step) {
    double time = (step + 1) * collumn.getTimeStep();

    // MPM algorithm
    collumn.setupMP();
    collumn.p2n();
    collumn.nodalEquilibrium();
    collumn.n2p();
    collumn.resetMesh();

    if (vtkInterval > 0 && ((step + 1) % vtkInterval == 0)) {
      collumn.exportVTKFrame(vtkDir, step + 1);
    }

    // Record results
    const auto mp = collumn.getMesh().getMP(tracked_mp).pos;
    hist << std::fixed << std::setprecision(10) << time << '\t' << mp.x()
         << '\t' << mp.y() << '\t' << collumn.getMPstress(tracked_mp) << '\t'
         << collumn.getMPstrain(tracked_mp) << '\n';

    if (step < 5 || step % 1000 == 0 || step == collumn.getNumSteps() - 1) {
      std::cout << std::fixed << std::setprecision(6);
      std::cout << "t=" << time << " | MP[" << tracked_mp << "]: x=" << mp.x()
                << " y=" << mp.y()
                << " | stress_xx=" << collumn.getMPstress(tracked_mp) << '\n';
    }
  }

  hist.close();
  std::cout << "\n✓ Results saved to: mpm2Da_history.txt\n";
  std::cout << "✓ Computation time: " << timer.elapsed() << " s\n";
  return 0;
}
