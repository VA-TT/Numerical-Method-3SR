#ifndef MATERIAL_POINT_METHOD_H
#define MATERIAL_POINT_METHOD_H

#include "Matrix.h"
#include "Mesh.h"
#include "Vector.h"
#include "gaussQuadrature.h"
#include "parentElement.h"
#include "physicConstants.h"
#include <cassert>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>

template <typename T, Index nNodes, Index nMPperEle> class MPM1D {
private:
  // Physical properties
  T m_E{1.0};      // Module Young
  T m_length{1.0}; // Domain length
  T m_volume{1.0}; // Volume
  T m_rho{1000};   // Density
  T m_mass{1000};  // mass
  T m_G{0.0};      // Gravity Acceleration

  // Simulation properties
  T m_currentTime{0.0}; // Current time
  T m_dt{0.0};          // Time step
  T m_duration{10.0};   // Duration of simulation
  Index m_nSteps{0};    // Number of steps
  T m_xloc{0.0};        // Position of surveying
  T m_v0{0.0};          // Initial velocity

  // Analytical solution (if available)
  std::function<T(T)> m_analyticSolution;
  // Behavior law
  std::function<T(T)> m_law;
  //  Mesh
  Mesh1D<T> m_mesh{};

  // Nodes n
  Vector<T> mass_n{};
  // Vector<T> position_n{}; // Don't get used
  Vector<T> velocity_n{};
  Vector<T> acceleration_n{};
  Vector<T> momentum_n{};
  Vector<char> velocityConstrained_n{};
  Vector<char> accelerationConstrained_n{};
  Vector<char> momentumConstrained_n{};
  Vector<char> forceConstrained_n{};

  // Constraint values (so constraints can be set anytime and enforced later)
  Vector<T> velocityConstraintValue_n{};
  Vector<T> accelerationConstraintValue_n{};
  Vector<T> momentumConstraintValue_n{};
  Vector<T> forceConstraintValue_n{};
  Vector<T> bodyForce_n{}, tractionForce_n{};
  // Nodal external forces
  Vector<double> forceExternal_n{}, forceInternal_n{}, totalForce_n{};

  // Material Points p
  Vector<T> volume_p{};   // Volume
  Vector<T> mass_p{};     // Mass
  Vector<T> position_p{}; // Position
  Vector<T> velocity_p{}; // Velocity
  Vector<T> momentum_p{}; // Momentum
  Vector<T> stress_p{}, strain_p{}, strain_rate_p{},
      dStrain_p{}; // Stress + strain

public:
  // Constructor
  MPM1D(T E, T rho, T length, T v0, T dt, T duration, T xloc)
      : m_E{E}, m_rho{rho}, m_length{length}, m_v0{v0}, m_dt{dt},
        m_duration{duration}, m_xloc{xloc}, m_mass{rho * length},
        m_volume{length * 1.0}, m_mesh{Mesh1D<T>{length, nNodes, nMPperEle}} {
    // Check critical time
    T c = std::sqrt(m_E / rho);
    T dt_crit = m_length / c;
    assert((dt_crit / 10.0) >= dt &&
           "Time step isn't satisfied CFL condition (too big)");
    m_nSteps = static_cast<Index>(duration / dt);

    // Initialize nodal vectors
    mass_n.resize(nNodes, T{});
    // position_n.resize(nNodes, T{});
    velocity_n.resize(nNodes, T{});
    acceleration_n.resize(nNodes, T{});
    momentum_n.resize(nNodes, T{});

    velocityConstrained_n.resize(nNodes, 0);
    accelerationConstrained_n.resize(nNodes, 0);
    momentumConstrained_n.resize(nNodes, 0);
    forceConstrained_n.resize(nNodes, 0);
    velocityConstraintValue_n.resize(nNodes, T{});
    accelerationConstraintValue_n.resize(nNodes, T{});
    momentumConstraintValue_n.resize(nNodes, T{});
    forceConstraintValue_n.resize(nNodes, T{});

    bodyForce_n.resize(nNodes, T{});
    tractionForce_n.resize(nNodes, T{});
    forceExternal_n.resize(nNodes, T{});
    forceInternal_n.resize(nNodes, T{});
    totalForce_n.resize(nNodes, T{});

    // Initialize MP vectors
    Index nMPs = m_mesh.getNumMPs();
    volume_p.resize(nMPs, m_volume / nMPs);
    mass_p.resize(nMPs, m_mass / nMPs); // MP's mass is constant
    position_p.resize(nMPs, T{});
    velocity_p.resize(nMPs, m_v0); // Initial velocity
    momentum_p.resize(nMPs, T{});
    stress_p.resize(nMPs, T{});
    strain_p.resize(nMPs, T{});
    strain_rate_p.resize(nMPs, T{});
    dStrain_p.resize(nMPs, T{});

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
  T getG() const { return m_G; }
  T getRho() const { return m_rho; }
  T getMass() const { return m_mass; }
  T getLength() const { return m_length; }
  T getVolume() const { return m_volume; }
  T getCurrentTime() const { return m_currentTime; }
  T getTimeStep() const { return m_dt; }
  T getDuration() const { return m_duration; }
  T getNumSteps() const { return m_nSteps; }
  T getSurveyLoc() const { return m_xloc; }
  T getIniVelo() const { return m_v0; }

  const Mesh1D<T> &getMesh() const { return m_mesh; }
  Index getNumNodes() const { return m_mesh.getNumNodes(); }
  Index getNumElements() const { return m_mesh.getNumElements(); }
  Index getNumMps() const { return m_mesh.getNumMPs(); }

  T getNodalMass(Index i) const { return mass_n[i]; }
  T getNodalVelocity(Index i) const { return velocity_n[i]; }
  T getNodalMomentum(Index i) const { return momentum_n[i]; }
  T getNodalExtForce(Index i) const { return forceExternal_n[i]; }
  T getNodalIntForce(Index i) const { return forceInternal_n[i]; }
  T getNodalTotalForce(Index i) const { return totalForce_n[i]; }

  T getMPvolume(Index p) const { return volume_p[p]; }
  T getMPmass(Index p) const { return mass_p[p]; }
  T getMPvelocity(Index p) const { return velocity_p[p]; }
  T getMPposition(Index p) const { return position_p[p]; }
  T getMPmomentum(Index p) const { return momentum_p[p]; }
  T getMPstrain(Index p) const { return strain_p[p]; }
  T getMPstrainRate(Index p) const { return strain_rate_p[p]; }
  T getMPdStrain(Index p) const { return dStrain_p[p]; }
  T getMPstress(Index p) const { return stress_p[p]; }

  // Setters
  void setE(T E) { m_E = E; }
  void setG(T G) { m_G = G; } // Set G if considering gravity
  void setAnalyticSolution(std::function<T(T)> sol) {
    m_analyticSolution = sol;
  }
  void setComportmentLaw(std::function<T(T)> law) { m_law = law; }

  // Setters: store constraint values + mark constrained
  void setNodalVeloConstraint(Index i, T value) {
    velocityConstraintValue_n[i] = value;
    velocityConstrained_n[i] = 1;
  }
  void setNodalAccConstraint(Index i, T value) {
    accelerationConstraintValue_n[i] = value;
    accelerationConstrained_n[i] = 1;
  }
  void setNodalMomentumConstraint(Index i, T value) {
    momentumConstraintValue_n[i] = value;
    momentumConstrained_n[i] = 1;
  }
  void setNodalForceConstraint(Index i, T value) {
    forceConstraintValue_n[i] = value;
    forceConstrained_n[i] = 1;
  }

  // Apply stored constraints to current nodal state
  void applyNodalVeloConstraint() {
    for (Index i{0}; i < getNumNodes(); ++i) {
      if (velocityConstrained_n[i] != 0) {
        velocity_n[i] = velocityConstraintValue_n[i];
      }
    }
  }
  void applyNodalAccConstraint() {
    for (Index i{0}; i < getNumNodes(); ++i) {
      if (accelerationConstrained_n[i] != 0) {
        acceleration_n[i] = accelerationConstraintValue_n[i];
      }
    }
  }
  void applyNodalMomentumConstraint() {
    for (Index i{0}; i < getNumNodes(); ++i) {
      if (momentumConstrained_n[i] != 0) {
        momentum_n[i] = momentumConstraintValue_n[i];
      }
    }
  }
  void applyNodalForceConstraint() {
    for (Index i{0}; i < getNumNodes(); ++i) {
      if (forceConstrained_n[i] != 0) {
        totalForce_n[i] = forceConstraintValue_n[i];
      }
    }
  }

  void setupMP() {
    position_p = m_mesh.getMPCoords(); // Initialize position of MPs
    m_mesh.updateMPElementIds();
    m_mesh.activateNodes();
  }

  void p2n() {
    Index nMPs = m_mesh.getNumMPs();
    for (Index p = 0; p < nMPs; ++p) {
      volume_p[p] = m_volume / nMPs;
      // mass_p[p] = m_mass / nMPs;
      momentum_p[p] = mass_p[p] * velocity_p[p];
    }
    // Map nodal mass + momentum
    for (Index p{0}; p < nMPs; ++p) {
      Index e = m_mesh.getMPElementId(p);
      if (e != -1) {
        T x_p = position_p[p];
        Index n1 = m_mesh.getEleConnectivity(e)[0];
        Index n2 = m_mesh.getEleConnectivity(e)[1];
        auto [x1, x2] = m_mesh.getElementNodes(e);
        T xi = parentCoor(x_p, x1, x2);
        mass_n[n1] += N1_ref(xi) * mass_p[p];
        mass_n[n2] += N2_ref(xi) * mass_p[p];
        momentum_n[n1] += N1_ref(xi) * momentum_p[p];
        momentum_n[n2] += N2_ref(xi) * momentum_p[p];
      }
    }

    applyNodalMomentumConstraint();
  }

  void nodalEquilibrium() {
    // f^ext = b + t
    for (Index p{0}; p < m_mesh.getNumMPs(); ++p) {
      Index e = m_mesh.getMPElementId(p);
      if (e != -1) {
        T x_p = position_p[p];
        Index n1 = m_mesh.getEleConnectivity(e)[0];
        Index n2 = m_mesh.getEleConnectivity(e)[1];
        auto [x1, x2] = m_mesh.getElementNodes(e);
        T xi = parentCoor(x_p, x1, x2);
        bodyForce_n[n1] += m_G * N1_ref(xi) * mass_p[p];
        bodyForce_n[n2] += m_G * N2_ref(xi) * mass_p[p];
        // Traction force t_i (to be implemented)
        forceInternal_n[n1] -= volume_p[p] * dN1_dx(x1, x2) * stress_p[p];
        forceInternal_n[n2] -= volume_p[p] * dN2_dx(x1, x2) * stress_p[p];
      }
    }
    forceExternal_n = bodyForce_n + tractionForce_n;
    totalForce_n = forceExternal_n + forceInternal_n;

    // Enforce any stored nodal force constraints after assembly
    applyNodalForceConstraint();
  }

  void n2p() {
    // Update momentum at nodes
    for (Index i{0}; i < getNumNodes(); ++i) {
      if (m_mesh.isActiveNode(i)) {
        acceleration_n[i] = totalForce_n[i] / mass_n[i];
        momentum_n[i] += totalForce_n[i] * m_dt;
        // velocity_n[i] += acceleration_n[i] * m_dt;
        // applyNodalVeloConstraint();
        // position_n[i] += velocity_n[i] * m_dt; //Grid nodes don't change
        // position?
      }
    }

    // Update particle position and velocity

    // Map back to MPs
    for (Index p{0}; p < m_mesh.getNumMPs(); ++p) {
      Index e = m_mesh.getMPElementId(p);
      if (e != -1) {
        T x_p = position_p[p];
        Index n1 = m_mesh.getEleConnectivity(e)[0];
        Index n2 = m_mesh.getEleConnectivity(e)[1];
        auto [x1, x2] = m_mesh.getElementNodes(e);
        T xi = parentCoor(x_p, x1, x2);
        velocity_p[p] += m_dt * N1_ref(xi) * totalForce_n[n1] / mass_n[n1] +
                         m_dt * N2_ref(xi) * totalForce_n[n2] / mass_n[n2];
        // position_p[p] += velocity_p[p] * m_dt;
        position_p[p] += (N1_ref(xi) * momentum_n[n1] / mass_n[n1] +
                          N2_ref(xi) * momentum_n[n2] / mass_n[n2]) *
                         m_dt;
        momentum_p[p] = mass_p[p] * velocity_p[p];
      }
    }

    // Map updated nodal velocity to calculate strain
    // for (Index i{0}; i < getNumNodes(); ++i) {
    //   if (m_mesh.isActiveNode(i)) {
    //     velocity_n[i] = momentum_p[p] * velocity_p[p] * m_dt;
    //   }
    // }

    // Update Nodal velocity to calculate strain
    // for (Index i{0}; i < getNumNodes(); ++i) {
    //   if (m_mesh.isActiveNode(i)) {
    //     velocity_n[i] = momentum_n[i] / mass_n[i];
    //   }
    // }
    for (Index p{0}; p < m_mesh.getNumMPs(); ++p) {
      Index e = m_mesh.getMPElementId(p);
      if (e != -1) {
        T x_p = position_p[p];
        Index n1 = m_mesh.getEleConnectivity(e)[0];
        Index n2 = m_mesh.getEleConnectivity(e)[1];
        auto [x1, x2] = m_mesh.getElementNodes(e);
        T xi = parentCoor(x_p, x1, x2);
        velocity_n[n1] += mass_p[p] * velocity_p[p] * N1_ref(xi) / mass_n[n1];
        velocity_n[n2] += mass_p[p] * velocity_p[p] * N2_ref(xi) / mass_n[n2];
      }
    }
    applyNodalVeloConstraint();

    for (Index p{0}; p < m_mesh.getNumMPs(); ++p) {
      Index e = m_mesh.getMPElementId(p);
      if (e != -1) {
        T x_p = position_p[p];
        Index n1 = m_mesh.getEleConnectivity(e)[0];
        Index n2 = m_mesh.getEleConnectivity(e)[1];
        auto [x1, x2] = m_mesh.getElementNodes(e);
        T xi = parentCoor(x_p, x1, x2);
        // Attention: x1-x2 belong to nodes (their positions don't get
        // updated), while the velocity is measured at MPs (updated at t+dt)
        strain_rate_p[p] =
            dN1_dx(x1, x2) * velocity_n[n1] + dN2_dx(x1, x2) * velocity_n[n2];
        dStrain_p[p] = strain_rate_p[p] * m_dt;
        strain_p[p] += dStrain_p[p];

        // Constitutive law:
        if (m_law) {
          stress_p[p] += m_law(dStrain_p[p]);
        } else {
          stress_p[p] += m_E * dStrain_p[p]; // Default linear elastic
        }
      }
    }

    m_mesh.setMPCoords(position_p); // Saving the updated position to Mesh

    // Update volume
    // volume_p[p] *= (1.0 + dStrain_p[p]);
  }

  void resetMesh() {
    mass_n.resetZero();
    momentum_n.resetZero();
    velocity_n.resetZero();
    acceleration_n.resetZero();

    bodyForce_n.resetZero();
    tractionForce_n.resetZero();
    forceExternal_n.resetZero();
    forceInternal_n.resetZero();
    totalForce_n.resetZero();
    m_mesh.nodalReset();
  }

  void applyBC() {
    // Apply boundary conditions
  }

  void timeIntegration() {}

  void compareAnalytic() {
    if (!m_analyticSolution) {
      std::cout << "\nNo analytical solution provided.\n";
      return;
    }

    std::cout << "\n=== Comparison with Analytical Solution ===\n";
    std::cout << std::setw(10) << "Time" << std::setw(15) << "x_numerical"
              << std::setw(15) << "x_exact" << std::setw(15) << "Error" << '\n';
    std::cout << std::string(55, '-') << '\n';

    // Find MP closest to m_xloc
    T min_dist = std::abs(position_p[0] - m_xloc);
    Index closest_p = 0;
    for (Index p = 1; p < m_mesh.getNumMPs(); ++p) {
      T dist = std::abs(position_p[p] - m_xloc);
      if (dist < min_dist) {
        min_dist = dist;
        closest_p = p;
      }
    }
    T x_numerical = position_p[closest_p];

    T x_exact = m_analyticSolution(m_currentTime);
    T error = std::abs(x_numerical - x_exact);

    std::cout << std::setw(10) << std::fixed << std::setprecision(4)
              << m_currentTime << std::setw(15) << std::setprecision(6)
              << x_numerical << std::setw(15) << x_exact << std::setw(15)
              << std::scientific << error << '\n';
  }

  void exportResult(const std::string &filename = "mpm1D_results.txt") {
    std::cout << "\n=== Exporting Results ===\n";
    std::ofstream file(filename);
    if (!file.is_open()) {
      std::cerr << "Error: Cannot open file " << filename << '\n';
      return;
    }

    file << "# MPM 1D Results\n";
    file << "# E = " << m_E << ", rho = " << m_rho << ", L = " << m_length
         << '\n';
    file << "# MPs = " << m_mesh.getNumMPs() << ", Nodes = " << nNodes << '\n';
    file << "# Columns: MP, x, v, strain, stress\n";
    file << std::setw(8) << "MP" << std::setw(15) << "x" << std::setw(15) << "v"
         << std::setw(15) << "strain" << std::setw(15) << "stress" << '\n';

    for (Index p = 0; p < m_mesh.getNumMPs(); ++p) {
      file << std::setw(8) << p << std::setw(15) << std::fixed
           << std::setprecision(6) << position_p[p] << std::setw(15)
           << velocity_p[p] << std::setw(15) << strain_p[p] << std::setw(15)
           << stress_p[p] << '\n';
    }

    file.close();
    std::cout << "Results exported to: " << filename << '\n';
  }
};

#endif // MATERIAL_POINT_METHOD_H
