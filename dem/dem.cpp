#include "../library/DEM.h"
#include <iostream>

int main() {
  constexpr Index ngl = 5;
  constexpr Index ngh = 5;
  double density = 2700.0; // density of sand particles (2700 kg/m^3)
  double rmin = 0.5e-3;    // smallest diameter is 1 mm
  double rmax = 1e-3;      // largest diameter is 2 mm
  double kn = 1e4;
  double kt = kn;
  double muy = 0.5;
  double cohesion = 0.0; // density of sand particles (2700 kg/m^3)
  double viscoRate = 0.95;

  double dt = 1e-5;
  double duration = 0.2;

  DEM2D<double, ngl, ngh> rigidWall(density, rmin, rmax, kn, kt, muy, cohesion,
                                    viscoRate, dt, duration);
  rigidWall.setG();
  rigidWall.setAlpha(0.01);
  rigidWall.setDeta();
  rigidWall.runFreeFall("data2D", 20);

  std::cout << "Done. Open dem/data2D/particles_*.vtk in ParaView.\n";
  return 0;
}