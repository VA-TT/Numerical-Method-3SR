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
  Index m_nSteps{0.0};  // Steps of simulation
  T m_xloc{0.0};        // Position of surveying
  T m_v0{0.0};          // Initial velocity

  // Analytical solution (if available)
  std::function<T(T)> m_analyticSolution;

  //  Mesh
  Mesh1D<T> m_mesh{};
  Index m_nMPs{};

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
    m_nMPs = m_mesh.getNumMPs();
    m_mesh.print();
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
  T getXloc() const { return m_xloc; }
  T getIniVelo() const { return m_v0; }

  const Mesh1D<T> &getMesh() const { return m_mesh; }
  Index getNumNodes() const { return m_mesh.getNumNodes(); }
  Index getNumElements() const { return m_mesh.getNumElements(); }
  Index getNumMps() const { return m_mesh.getNumMPs(); }

  Vector<T> getNodalMass(Index i) { return mass_n[i]; }
  Vector<T> getNodalVelocity(Index i) { return velocity_n[i]; }
  Vector<T> getNodalMomentum(Index i) { return momentum_n[i]; }
  Vector<T> getNodalExtForce(Index i) { return f_ext_n[i]; }
  Vector<T> getNodalIntForce(Index i) { return f_int_n[i]; }
  Vector<T> getNodalTotalForce(Index i) { return f_total_n[i]; }

  Vector<T> getMPvolume(Index p) { return volume_p[p]; }
  Vector<T> getMPmass(Index p) { return mass_p[p]; }
  Vector<T> getMPvelocity(Index p) { return velocity_p[p]; }
  Vector<T> getMPposition(Index p) { return position_p[p]; }
  Vector<T> getMPmomentum(Index p) { return momentum_p[p]; }
  Vector<T> getMPstrain(Index p) { return strain_p[p]; }
  Vector<T> getMPstrainRate(Index p) { return strain_rate_p[p]; }
  Vector<T> getMPdStrain(Index p) { return dStrain_p[p]; }
  Vector<T> getMPstress(Index p) { return stress[p]; }
  Vector<T> getMPtotalForce(Index p) { return f_total_n[p]; }

  // Setters
  void setE(T E) { m_E = E; }
  void setAnalyticSolution(std::function<T(T)> sol) {
    m_analyticSolution = sol;
  }

  void setupMP() {                               // Set up MPs
    volume_p = Vector<T>(nMPs, m_volume / nMPs); // m_volume might change
    velocity_p = Vector<T>(nMPs, m_v0);
    momentum_p = mass_p * m_v0;
  }

  void p2n() {
    for (Index p{0}; p < m_nMPs; ++p) {
      Index e = findCageID(p);
      if (e != -1) {
        Index n1 = m_mesh.m_connectivity[e][0];
        Index n2 = m_connectivity[e][1];
      }
    }
  }

  void computeNodalForce() { // f_i,ext = b_i + t_i (body force + nodal force)
    f_ext_i.resetZero();
    // f_i,int = V_p*B_i(x_p)*stress_p (body force + nodal force)
    f_int_i.resetZero();
    for (Index p{0}; p < nMPs; p++) {
      for (Index i{0}; i < nNodes; i++) {
        f_int_i[i] += -B_p[i] * V_p *
                      modelParameters::stress_p; // m_i = Σ_p N_i(x_p) * mass_p
      }
    }
    f_total_i = f_int_i + f_ext_i; // Unbalanced at this point}
  }

  void applyBC() {
    mv_i[0] = 0;      // v[0] = 0
    f_total_i[0] = 0; // a[0] = 0}
                      // Update nodal momentum
    mv_i += f_total_i * dt;

    void n2p() { // Update particle velocity: v_p^{n+1} = v_p^n + Δt * Σ_i
                 // N_i(x_p) * a_i
      // where a_i = f_i / m_i
      // Update particle position: x_p^{n+1} = x_p^n + Δt * Σ_i N_i(x_p) *
      // v_i^{n+1} where v_i^{n+1} = mv_i^{n+1} / m_i
      for (Index i{0}; i < nNodes; i++) {
        v_p += N_p[i] * (f_total_i[i] / m_i[i]) * dt;
        x_p += N_p[i] * (mv_i[i] / m_i[i]) * dt;
      }

      // Check for NaN or out of bounds
      if (std::isnan(x_p) || std::isnan(v_p)) {
        std::cerr << "ERROR at step " << step << ": NaN detected!\n";
        std::cerr << "  x_p = " << x_p << ", v_p = " << v_p << "\n";
        std::cerr << "  m_i = " << m_i << "\n";
        std::cerr << "  N_p = " << N_p << "\n";
        break;
      }

      if (x_p < a || x_p > b) {
        std::cerr << "WARNING at step " << step
                  << ": Particle out of domain!\n";
        std::cerr << "  x_p = " << x_p << " (domain: [" << a << ", " << b
                  << "])\n";
      }

      // Update particle momentum
      mv_p = mass_p * v_p;
      // Update nodal velocity
      for (Index i{0}; i < nNodes; i++) {
        if (m_i[i] > 1e-12) {
          v_i[i] = mass_p * v_p * N_p[i] / m_i[i];
          v_i[i] = mv_i[i] / m_i[i]; // Causing losing in energy
        } else {
          v_i[i] = 0.0;
        }
      }
      // Apply BC
      v_i[0] = 0;

      // Compute strain rate and stress
      strain_rate_p = 0.0; // Reset strain rate
      for (Index i{0}; i < nNodes; i++) {
        strain_rate_p += B_p[i] * v_i[i]; // Mapping strain rate from the node
      }
      dStrain_p = strain_rate_p * dt;
      stress_p += E * dStrain_p; // Elastic
    }
    std::cout << "Mapped mass to nodes (last step): " << m_i << '\n';
  }

  void resetMesh() {
    // Reset nodal masses and momentum to zero
    for (Index i{0}; i < nNodes; i++) {
      m_i[i] = mv_i[i] = 0.0;
      m_mesh.resetMesh();
    }
  }

  void timeIntegration() {
    for (Index step{0}; step < nSteps; step++) {
      m_currentTime = step * m_dt;
      times.push_back(current_time);
      p2n();
      computeNodalForce();
      applyBC();
      n2p();
      resetMesh();
    }
  }

  void compareAnalytic() {
    std::cout << "\nComparison at t=0:\n";
    std::cout << "  Position: analytical=" << positions[0]
              << ", MPM=" << positions_mpm[0]
              << ", error=" << constexpr_fabs(positions[0] - positions_mpm[0])
              << "\n";
    std::cout << "  Velocity: analytical=" << velocities[0]
              << ", MPM=" << velocities_mpm[0]
              << ", error=" << constexpr_fabs(velocities[0] - velocities_mpm[0])
              << "\n";
    std::cout << "\nComparison at final time:\n";
    std::cout << "  Position: analytical=" << positions.back()
              << ", MPM=" << positions_mpm.back() << ", error="
              << constexpr_fabs(positions.back() - positions_mpm.back())
              << "\n";
    std::cout << "  Velocity: analytical=" << velocities.back()
              << ", MPM=" << velocities_mpm.back() << ", error="
              << constexpr_fabs(velocities.back() - velocities_mpm.back())
              << "\n";
  }

  void exportResult(
      std::ofstream txtFile("mpm1D.txt");
      txtFile << std::fixed << std::setprecision(6);
      txtFile << "#Time \t x_analytical \t v_analytical \t x_MPM \t v_MPM \t "
                 "error_x \t error_v\n";
      for (Index i = 0; i < nSteps; ++i) {
        double error_x = constexpr_fabs(positions[i] - positions_mpm[i]);
        double error_v = constexpr_fabs(velocities[i] - velocities_mpm[i]);
        txtFile << times[i] << " \t " << positions[i] << " \t " << velocities[i]
                << " \t " << positions_mpm[i] << " \t " << velocities_mpm[i]
                << " \t " << error_x << " \t " << error_v << "\n";
      } txtFile.close();) {}
};

#endif // MPM_H

#include "../library/clock.h"
#include <iostream>
int main() {
  Timer t;
  double E = 4 * constants::pi * constants::pi;
  MPM1D<double, 2, 1> beam1D(
      E, 1.0, 1.0, 0.1, 0.01, 10.0,
      0.5); // MPM1D<type, nPoints, nMPperEle>(E,rho,length,v0,dt,duration,
            // xloc)
  double omega =
      1 / (beam1D.getLength()) * std::sqrt(beam1D.getE() / beam1D.getRho());
  beam1D.setAnalyticSolution([]() {
    return beam1D.getXloc() / std::exp(beam1D.getIniVelo() / (beam1D.getLength()*w) * std::sin(w * beam1D.getCurrentTime())};
std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";
}
