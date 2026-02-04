#ifndef MATERIAL_POINT_METHOD_H
#define MATERIAL_POINT_METHOD_H

#include "Matrix.h"
#include "Mesh.h"
#include "Vector.h"
#include "gaussQuadrature.h"
#include "parentElement.h"
#include <cassert>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>

// 1D Finite Element Method
//  EA U_xx + f(x) = 0, 0 < x < L
//  Boundary Conditions:
//    - Dirichlet: u(x) = prescribed value (essential BC)
//    - Neumann: EA * u'(x) = prescribed force (natural BC)

template <typename T, Index nNodes, Index nMPperEle> class MPM1D {
private:
  // Physical properties
  T m_E{10.0};     // Module Young
  T m_length{1.0}; // Domain length
  T m_volume{1.0}; // Volume
  T m_rho{1000};   // Density

  // Simulation properties
  T m_currentTime{0.0}; // Current time
  T m_dt{0.0};          // Time step
  T m_duration{0.0};    // Duration of simulation
  Inden m_nSteps{0.0};  // Steps of simulation
  T m_xloc{0.0};        // Position of surveying
  T m_v0{0.0};          // Position of surveying

  //  Mesh
  Mesh1D<T> m_mesh{};

  // Nodes n
  Vector<T> mass_n(nNodes);
  Vector<T> position_n(nNodes);
  Vector<T> velocity_n(nNodes);
  Vector<T> momentum_n(nNodes);
  // Nodal external forces
  Vector<double> f_ext_n(nNodes), f_int_n(nNodes), f_total_n(nNodes);

  // Material Points p
  Vector<T> volume_p((nNodes - 1) * nMPperEle);   // Volume
  Vector<T> mass_p((nNodes - 1) * nMPperEle);     // Mass
  Vector<T> position_p((nNodes - 1) * nMPperEle); // Position
  Vector<T> velocity_p((nNodes - 1) * nMPperEle); // Velocity
  Vector<T> momentum_p((nNodes - 1) * nMPperEle); // Momentum
  Vector<T> N_p((nNodes - 1) * nMPperEle),
      B_p((nNodes - 1) * nMPperEle); // Shape function
  Vector<T> stress_p((nNodes - 1) * nMPperEle),
      strain_p((nNodes - 1) * nMPperEle),
      strain_rate_p((nNodes - 1) * nMPperEle),
      dStrain_p((nNodes - 1) * nMPperEle); // Stress + strain

  // Analytical solution (if available)
  std::function<T(T)> m_analyticSolution;

public:
  // Constructor
  MPM1D(T E, T rho, T length, T v0, T dt, T duration, T xloc)
      : m_E{E}, m_rho{rho}, m_length{length}, m_v0{v0}, m_dt{dt},
        m_duration{duration}, m_xloc{xloc} {
    // Check critical time
    double c{std::sqrt{E / rho}};
    double dt_crit{m_length / c};
    assert((dt_crit / 10.0) >= dt &&
           "Time step isn't satisfied CFL condition (too big)");
    m_nSteps = duration / dt;

    // System set up
    m_volume = length * 1.0 * 1.0 * rho; // 1D: B = H = 1.0

    // Set up mesh
    m_mesh = Mesh1D<T>{length, nNodes, nMPperEle};
    m_mesh.print();

    // Set up MPs
    Index nMPs = m_mesh.getNumMPs;
    volume_p = Vector<T>(nMPs, m_volume / nMPs);
    mass_p = rho * volume_p;
    velocity_p = Vector<T>(nMPs, v0);
    momentum_p = mass_p * v0;
  };

  // Other defaults
  MPM1D() = default;
  MPM1D(const MPM1D &) = default;
  MPM1D(MPM1D &&) = default;
  MPM1D &operator=(const MPM1D &) = default;
  MPM1D &operator=(MPM1D &&) = default;
  ~MPM1D() = default;

  // Getters
  T getE() const { return m_E; }
  T getRho() const { return m_rho; }
  T getLength() const { return m_length; }
  T getVolume() const { return m_volume; }
  T getCurrentTime() const { return m_currentTime; }
  T getTimeStep() const { return m_dt; }
  T getDuration() const { return m_duration; }
  T getNumSteps() const { return m_nSteps; }
  T getSurLoc() const { return m_xloc; }
  T getIniVelo() const { return m_v0; }

  const Mesh1D<T> &getMesh() const { return m_mesh; }
  Index getNumNodes() const { return m_mesh.getNumNodes(); }
  Index getNumElements() const { return m_mesh.getNumElements(); }
  Index getNumMps() const { return m_mesh.getNumMPs(); }

  getNodalMass(Index i) { return mass_n[i]; }
  getNodalVelocity(Index i) { return velocity_n[i]; }
  getNodalMomentum(Index i) { return momentum_n[i]; }
  getNodalExtForce(Index i) { return f_ext_n[i]; }
  getNodalIntForce(Index i) { return f_int_n[i]; }
  getNodalTotalForce(Index i) { return f_total_n[i]; }

  getMPvolume(Index p) { return volume_p[p]; }
  getMPmass(Index p) { return mass_p[p]; }
  getMPvelocity(Index p) { return velocity_p[p]; }
  getMPposition(Index p) { return position_p[p]; }
  getMPmomentum(Index p) { return momentum_p[p]; }
  getMPstrain(Index p) { return strain_p[p]; }
  getMPstrainRate(Index p) { return strain_rate_p[p]; }
  getMPdStrain(Index p) { return dStrain_p[p]; }
  getMPstress(Index p) { return stress[p]; }
  getMPtotalForce(Index p) { return f_total_n[p]; }

  // Setters
  void setE(T E) { m_E = E; }
  void setAnalyticSolution(std::function<T(T)> sol) {
    m_analyticSolution = sol;
  }

  void p2n() {}

  void computeNodalForce() {}

  void applyBC() {}

  void n2p() {}

  void timeIntegration() {
    for (Index step{0}; step < nSteps; step++) {
      m_currentTime = step * m_dt;
      times.push_back(current_time);
      p2n();
      computeNodalForce();
      applyBC();
      n2p();
    }
  }

  void compareAnalytic() {}

  void exportResult(const std::string &filename = "MPM1D_results.txt") {}
};

#endif // MPM_H

#include "../library/clock.h"
#include <iostream>
int main() {
  Timer t;
  std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";
}
