#include "../library/MPM.h"
#include "../library/clock.h"
#include <iostream>

int main() {
  Timer t;
  double E = 4 * constants::pi * constants::pi;
  MPM1D<double, 2, 1> beam1D(E, 1.0, 1.0, 0.1, 0.001, 1.0, 0.5);
  // MPM1D<type, nNodes, nMPperEle>(E,rho,length,v0,dt,duration,xloc)

  double L = beam1D.getLength();
  double v0 = beam1D.getIniVelo();
  double omega = (1.0 / L) * std::sqrt(beam1D.getE() / beam1D.getRho());
  double xloc = beam1D.getXloc();

  // Analytical solution: x(t) = xloc * exp(v0/(L*omega) * sin(omega*t))
  beam1D.setAnalyticSolution([L, v0, omega, xloc](double t) {
    return xloc * std::exp((v0 / (L * omega)) * std::sin(omega * t));
  });

  std::cout << "Starting MPM simulation...\n";
  std::cout << "E = " << E << ", rho = 1.0, L = " << L << '\n';
  std::cout << "omega = " << omega << " rad/s\n";
  std::cout << "v0 = " << v0 << " m/s, dt = " << beam1D.getTimeStep()
            << " s\n\n";

  beam1D.timeIntegration();
  beam1D.compareAnalytic();
  beam1D.exportResult();

  std::cout << "\nTime elapsed: " << t.elapsed() << " seconds\n";
  return 0;
}
