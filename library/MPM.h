#ifndef MATERIAL_POINT_METHOD_H
#define MATERIAL_POINT_METHOD_H

#include "Elasticity.h"
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
#include <vector>

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
  Vector<T> forceExternal_n{}, forceInternal_n{}, forceTotal{};

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
    forceTotal.resize(nNodes, T{});

    // Initialize MP vectors
    Index nMPs = m_mesh.getNumMPs();
    volume_p.resize(nMPs, m_volume / nMPs);
    mass_p.resize(nMPs, m_mass / nMPs); // MP's mass is constant
    // position_p.resize(nMPs);
    position_p = m_mesh.getMPCoords(); // Initialize position of MPs
    velocity_p.resize(nMPs, m_v0);     // Initial velocity
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
  T getNodalTotalForce(Index i) const { return forceTotal[i]; }

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
  void setCurrentTime(T value) { m_currentTime = value; }
  void setE(T E) { m_E = E; }
  void setG(T G) { m_G = G; } // Set G if considering gravity
  void setAnalyticSolution(std::function<T(T)> sol) {
    m_analyticSolution = sol;
  }
  void setComportmentLaw(std::function<T(T)> law) { m_law = law; }
  void setMPvelocity(Index p, T value) { velocity_p[p] = value; }

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
        forceTotal[i] = forceConstraintValue_n[i];
      }
    }
  }

  void setupMP() {
    m_mesh.updateMPElementIds();
    m_mesh.activateNodes();
  }

  void p2n() {
    Index nMPs = m_mesh.getNumMPs();
    // Map nodal mass + momentum
    for (Index p{0}; p < nMPs; ++p) {
      momentum_p[p] = mass_p[p] * velocity_p[p];
      Index e = m_mesh.getMPelementID(p);
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
      Index e = m_mesh.getMPelementID(p);
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
    forceTotal = forceExternal_n + forceInternal_n;

    // Enforce any stored nodal force constraints after assembly
    applyNodalForceConstraint();
  }

  void n2p() {
    // Update momentum at nodes
    for (Index i{0}; i < getNumNodes(); ++i) {
      if (m_mesh.isActiveNode(i)) {
        momentum_n[i] += forceTotal[i] * m_dt;
        acceleration_n[i] = forceTotal[i] / mass_n[i];

        // velocity_n[i] += acceleration_n[i] * m_dt;
        // applyNodalVeloConstraint();
        // position_n[i] += velocity_n[i] * m_dt; //Grid nodes don't change
        // position?
      }
    }
    applyNodalAccConstraint(); // No need to apply AccConstraint here as force
                               // constraint is already applied

    // Update particle position and velocity

    // Map back to MPs
    for (Index p{0}; p < m_mesh.getNumMPs(); ++p) {
      Index e = m_mesh.getMPelementID(p);
      if (e != -1) {
        T x_p = position_p[p];
        Index n1 = m_mesh.getEleConnectivity(e)[0];
        Index n2 = m_mesh.getEleConnectivity(e)[1];
        auto [x1, x2] = m_mesh.getElementNodes(e);
        T xi = parentCoor(x_p, x1, x2);
        // Update velocity and position using FLIP style
        velocity_p[p] += (N1_ref(xi) * acceleration_n[n1] +
                          N2_ref(xi) * acceleration_n[n2]) *
                         m_dt;
        // Hybrid
        // T v_pic = N1 * velocity_n[n1] + N2 * velocity_n[n2];
        // T v_flip = velocity_p[p] + (N1 * a_n1 + N2 * a_n2) * dt;
        // velocity_p[p] = alpha * v_pic + '(1-alpha)' * v_flip;
        // position_p[p] += velocity_p[p] * m_dt;
        position_p[p] += (N1_ref(xi) * momentum_n[n1] / mass_n[n1] +
                          N2_ref(xi) * momentum_n[n2] / mass_n[n2]) *
                         m_dt;
        momentum_p[p] = mass_p[p] * velocity_p[p];
      }
    }

    for (Index p{0}; p < m_mesh.getNumMPs(); ++p) {
      Index e = m_mesh.getMPelementID(p);
      if (e != -1) {
        T x_p = position_p[p];
        Index n1 = m_mesh.getEleConnectivity(e)[0];
        Index n2 = m_mesh.getEleConnectivity(e)[1];
        auto [x1, x2] = m_mesh.getElementNodes(e);
        T xi = parentCoor(x_p, x1, x2);
        velocity_n[n1] += momentum_p[p] * N1_ref(xi) / mass_n[n1];
        velocity_n[n2] += momentum_p[p] * N2_ref(xi) / mass_n[n2];
      }
    }
    applyNodalVeloConstraint();

    for (Index p{0}; p < m_mesh.getNumMPs(); ++p) {
      Index e = m_mesh.getMPelementID(p);
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
        // Update volume
        volume_p[p] *= (1.0 + dStrain_p[p]);
      }
    }
  }

  void resetMesh() {
    m_mesh.setMPCoords(position_p); // Saving the updated MPs' position to Mesh

    mass_n.resetZero();
    momentum_n.resetZero();
    velocity_n.resetZero();
    acceleration_n.resetZero();

    bodyForce_n.resetZero();
    tractionForce_n.resetZero();
    forceExternal_n.resetZero();
    forceInternal_n.resetZero();
    forceTotal.resetZero();

    m_mesh.nodalReset();
  }

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

  void exportResult(const std::string &filename = "mpm1D_results.txt") {}
  void applyBC() {}
  void timeIntegration() {}
};

////////////////////////////////////////////
////////////////::: 2D ::://////////////////
////////////////////////////////////////////

template <typename T, T gridLength, T gridHeight, Index nx, Index ny, T MP_size>
class MPM2D {
private:
  // Kinetic, Potential, Initial Total and Dissipation energies
  T m_kinEnergy{}, m_potentEnergy{}, m_totalEnergy0{}, m_dissiEnergy{};
  // Material properties
  T m_E{100000.0}; // Module Young
  T m_nu{1.0};     // Poisson ratio
  T m_rho{1000};   // Density
  T m_phi{30.0};   // Internal Friction angle
  T m_mu{30.0};    // Material-Boundary Friction
  T m_c{1.0};      // cohesion
  T m_K0{0.5};     // Initial earth pressure coefficient

  // MP domain
  Vector<T> minCorner{}, maxCorner{};
  T m_pLength{}, m_pHeight{};
  T m_volume{1.0}, m_volume0{1.0}; // Volume
  T m_mass{1000};                  // mass
  T m_G{0.0};                      // Gravity Acceleration
  Vector<T> m_v0{T{}, T{}};        // Initial velocity

  // Simulation properties
  T m_currentTime{0.0};  // Current time
  T m_dt{0.0};           // Time step
  T m_duration{10.0};    // Duration of simulation
  Index m_nSteps{0};     // Number of steps
  Index m_interval{100}; // Output interval

  // // Behavior law
  // std::function<T(T)> m_law;

  //  Mesh
  Mesh2D<T> m_mesh{};

  // Nodes n
  Vector<T> mass_n{};
  // Vector<T> position_n{}; // Don't get used
  Vector<Vector<T>> velocity_n{};
  Vector<Vector<T>> acceleration_n{};
  Vector<Vector<T>> momentum_n{};
  Vector<Vector<T>> displacement_n{}; // Displacement
  Vector<Matrix<T, 2, 2>> stress_n{}; // Displacement

  Vector<Vector<char>> velocityConstrained_n{};
  Vector<Vector<char>> accelerationConstrained_n{};
  Vector<Vector<char>> momentumConstrained_n{};
  Vector<Vector<char>> forceConstrained_n{};

  // Constraint values (so constraints can be set anytime and enforced later)
  Vector<Vector<T>> velocityConstraintValue_n{};
  Vector<Vector<T>> accelerationConstraintValue_n{};
  Vector<Vector<T>> momentumConstraintValue_n{};
  Vector<Vector<T>> forceConstraintValue_n{};
  Vector<Vector<T>> bodyForce_n{}, tractionForce_n{};
  // Nodal external forces
  Vector<Vector<T>> forceExternal_n{}, forceInternal_n{}, forceTotal{};

  // Material Points p
  Vector<T> minCorner, maxCorner{};
  Vector<T> volume_p{};           // Volume
  Vector<T> mass_p{};             // Mass
  Vector<Vector<T>> position_p{}; // Position
  Vector<Vector<T>> velocity_p{}; // Velocity
  Vector<Vector<T>> momentum_p{}; // Momentum

  Vector<Matrix<T, 2, 2>> stress_p{}, strain_p{}, deformGradient_p{},
      dStrain_p{}; // Stress + strain

public:
  // Constructor
  MPM2D(T E, T nu, T rho, T mu, T phi, T c, T K0, const Vector<T> &minCorner,
        const Vector<T> &maxCorner, T dt, T duration, T v0 = T{})
      : m_E{E}, m_nu{nu}, m_rho{rho}, m_mu{mu}, m_phi{phi}, m_c{c}, m_K0{K0},
        m_v0{v0}, m_dt{dt}, m_duration{duration},
        m_nSteps{static_cast<Index>(duration / dt)},
        m_pLength{maxCorner.x() - minCorner.x()},
        m_pHeight{maxCorner.y() - minCorner.y()},
        m_volume{constexpr_fabs(m_pLength * m_pHeight)}, m_mass{rho * m_volume},
        m_volume0{m_volume}, m_mesh{Mesh2D<T>{{gridLength, gridHeight},
                                              {nx, ny},
                                              minCorner,
                                              maxCorner,
                                              MP_size}} {

    // Check critical time
    T c = std::sqrt(m_E / rho);
    T dt_crit = (m_pLength < m_pHeight ? m_pLength : m_pHeight) / c;
    assert((dt_crit / 10.0) >= dt &&
           "dt doesn't satisfied CFL condition (too big)");

    // Initialize nodal vectors
    mass_n.resize(nNodes);
    // position_n.resize(nNodes, Vector<T>{});
    velocity_n.resize(nNodes, Vector<T>{0, 0});
    acceleration_n.resize(nNodes, Vector<T>{0, 0});
    momentum_n.resize(nNodes, Vector<T>{0, 0});
    displacement_n.resize(nNodes, Vector<T>{0, 0});
    stress_n.resize(nMPs, Matrix<T, 2, 2>::zero);

    velocityConstrained_n.resize(nNodes, Vector<char>{0, 0});
    accelerationConstrained_n.resize(nNodes, Vector<char>{0, 0});
    momentumConstrained_n.resize(nNodes, Vector<char>{0, 0});
    forceConstrained_n.resize(nNodes, Vector<char>{0, 0});
    velocityConstraintValue_n.resize(nNodes, Vector<T>{0, 0});
    accelerationConstraintValue_n.resize(nNodes, Vector<T>{0, 0});
    momentumConstraintValue_n.resize(nNodes, Vector<T>{0, 0});
    forceConstraintValue_n.resize(nNodes, Vector<T>{0, 0});

    bodyForce_n.resize(nNodes, Vector<T>{0, 0});
    tractionForce_n.resize(nNodes, Vector<T>{0, 0});
    forceExternal_n.resize(nNodes, Vector<T>{0, 0});
    forceInternal_n.resize(nNodes, Vector<T>{0, 0});
    forceTotal.resize(nNodes, Vector<T>{0, 0});

    // Initialize MP vectors
    Index nMPs = m_mesh.getNumMPs();
    volume_p.resize(nMPs, m_volume / nMPs);
    mass_p.resize(nMPs, m_mass / nMPs); // MP's mass is constant
    // position_p.resize(nMPs);
    position_p = m_mesh.getMPCoords();        // Initialize position of MPs
    momentum_p.resize(nMPs, Vector<T>{0, 0}); // Compute later
    velocity_p.resize(nMPs, m_v0);            // Initial velocity

    stress_p.resize(nMPs, Matrix<T, 2, 2>::zero);
    strain_p.resize(nMPs, Matrix<T, 2, 2>::zero);
    deformGradient_p.resize(nMPs, Matrix<T, 2, 2>::identity);
    dStrain_p.resize(nMPs, Matrix<T, 2, 2>::zero);

    for (Index p{0}; p < nMPs; ++p) {
      m_totalEnergy0 += m_G * position_p[p].y() * mass_p[p];
    }

    m_mesh.print();
  };

  // Other defaults
  MPM2D() = default;
  MPM2D(const MPM2D &) = default;
  MPM2D(MPM2D &&) = default;
  MPM2D &operator=(const MPM2D &) = default;
  MPM2D &operator=(MPM2D &&) = default;
  ~MPM2D() = default;

  // Getters
  T getE() const { return m_E; }
  T getG() const { return m_G; }
  T getRho() const { return m_rho; }
  T getMass() const { return m_mass; }
  T getVolume() const { return m_volume; }
  T getCurrentTime() const { return m_currentTime; }
  T getTimeStep() const { return m_dt; }
  T getDuration() const { return m_duration; }
  T getNumSteps() const { return m_nSteps; }
  T getIniVelo() const { return m_v0; }

  const Mesh2D<T> &getMesh() const { return m_mesh; }
  Index getNumNodes() const { return m_mesh.getNumNodes(); }
  Index getNumElements() const { return m_mesh.getNumElements(); }
  Index getNumMps() const { return m_mesh.getNumMPs(); }

  T getNodalMass(Index i) const { return mass_n[i]; }
  T getNodalVelocity(Index i) const { return velocity_n[i]; }
  T getNodalMomentum(Index i) const { return momentum_n[i]; }
  T getNodalExtForce(Index i) const { return forceExternal_n[i]; }
  T getNodalIntForce(Index i) const { return forceInternal_n[i]; }
  T getNodalTotalForce(Index i) const { return forceTotal[i]; }

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
  void setCurrentTime(T value) { m_currentTime = value; }
  void setE(T E) { m_E = E; }
  void setG(T G) { m_G = G; } // Set G if considering gravity

  void setComportmentLaw(std::function<T(T)> law) { m_law = law; }
  void setMPvelocity(Index p, T value) { velocity_p[p] = value; }

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
        forceTotal[i] = forceConstraintValue_n[i];
      }
    }
  }

  void applyFrictionalBC(const Vector<Index>& nodeIDs, Index dir_n, Index signDir_n) {
    dir_t = 1 - dir_n; //tangential diretion
    //Normal and tangential acceleration
    T acc_n = acceleration_n[nodeIDs];
  }

  void initializeStress() {
    T ymax = m_pHeight - MP_size / 2;
    for (Index p{0}; p < this->getNumMps(); ++p) {
      stress_p[p].yy() = -m_G * m_rho * (ymax - position_p[p].y());
      stress_p[p].xx() = m_K0 * stress_p[p].yy();
    }
  }

  void computeEnergy() {
    for (Index p{0}; p < m_mesh.getNumMPs(); ++p) {
      m_potentEnergy += m_G * mass_p[p] * position_p[p].y();
      m_kinEnergy +=
          T{0.5} * mass_p[p] * dotProduct(velocity_p[p], velocity_p[p]);
    }
    m_dissiEnergy = m_totalEnergy0 - m_potentEnergy - m_kinEnergy;
  }

  void setupMP() {
    m_mesh.activateNodes();
    m_mesh.activateElements();
    initializeStress();
    computeEnergy();
  }

  void p2n() {

    for (Index p{0}; p < m_mesh.getNumMPs(); ++p) {
      momentum_p[p].x() = mass_p[p] * velocity_p[p].x();
      momentum_p[p].y() = mass_p[p] * velocity_p[p].y();
    }

    // Map nodal mass + momentum
    for (Index e{0}; e < m_mesh.getNumElements(); ++e) {
      if (!m_mesh.isActiveElement(e)) {
        continue;
      }

      const auto &mps = m_mesh.getMPsInElement(e);
      if (mps.empty()) {
        continue;
      }

      const auto &conn = m_mesh.getEleConnectivity(e);
      const Index n1 = conn[0];
      const Index n2 = conn[1];
      const Index n3 = conn[2];
      const Index n4 = conn[3];

      const auto [x_nodes, y_nodes] = m_mesh.getElementNodes(e);

      for (const Index p : mps) {
        const T x_p = position_p[p].x();
        const T y_p = position_p[p].y();
        const auto [xi, eta] = parentCoor(x_p, y_p, x_nodes, y_nodes);

        mass_n[n1] += N1_ref(xi, eta) * mass_p[p];
        mass_n[n2] += N2_ref(xi, eta) * mass_p[p];
        mass_n[n3] += N3_ref(xi, eta) * mass_p[p];
        mass_n[n4] += N4_ref(xi, eta) * mass_p[p];

        momentum_n[n1] += N1_ref(xi, eta) * momentum_p[p];
        momentum_n[n2] += N2_ref(xi, eta) * momentum_p[p];
        momentum_n[n3] += N3_ref(xi, eta) * momentum_p[p];
        momentum_n[n4] += N4_ref(xi, eta) * momentum_p[p];
      }
    }

    applyNodalMomentumConstraint();
  }

  void nodalEquilibrium() {
    // f^ext = b + t

    // Map nodal mass + momentum
    for (Index e{0}; e < m_mesh.getNumElements(); ++e) {
      if (!m_mesh.isActiveElement(e)) {
        continue;
      }
      const auto &mps = m_mesh.getMPsInElement(e);
      if (mps.empty()) {
        continue;
      }

      const auto &conn = m_mesh.getEleConnectivity(e);
      const Index n1 = conn[0];
      const Index n2 = conn[1];
      const Index n3 = conn[2];
      const Index n4 = conn[3];

      const auto [x_nodes, y_nodes] = m_mesh.getElementNodes(e);

      for (const Index p : mps) {
        const T x_p = position_p[p].x();
        const T y_p = position_p[p].y();
        const auto [xi, eta] = parentCoor(x_p, y_p, x_nodes, y_nodes);

        // Gravity body force (y-direction) : G should be negative!
        bodyForce_n[n1].y() += m_G * N1_ref(xi, eta) * mass_p[p];
        bodyForce_n[n2].y() += m_G * N2_ref(xi, eta) * mass_p[p];
        bodyForce_n[n3].y() += m_G * N3_ref(xi, eta) * mass_p[p];
        bodyForce_n[n4].y() += m_G * N4_ref(xi, eta) * mass_p[p];

        // Traction force t_i (to be implemented)

        // Internal force
        const auto [dN_dx, dN_dy] = dNdxdy(xi, eta, x_nodes, y_nodes);
        forceInternal_n[n1] -=
            volume_p[p] * (stress_p[p] * Vector<T>{dN_dx[0], dN_dy[0]});
        forceInternal_n[n2] -=
            volume_p[p] * (stress_p[p] * Vector<T>{dN_dx[1], dN_dy[1]});
        forceInternal_n[n3] -=
            volume_p[p] * (stress_p[p] * Vector<T>{dN_dx[2], dN_dy[2]});
        forceInternal_n[n4] -=
            volume_p[p] * (stress_p[p] * Vector<T>{dN_dx[3], dN_dy[3]});
      }
    };

    for (Index i{0}; i < m_mesh.getNumNodes(); ++i) {
      if (m_mesh.isActiveNode(i)) {
        forceExternal_n[i] = bodyForce_n[i] + tractionForce_n[i];
        forceTotal[i] = forceExternal_n[i] + forceInternal_n[i];
      }
    }

    // Enforce any stored nodal force constraints after assembly
    applyNodalForceConstraint();
  }

  void n2p() {
    // Update momentum at nodes
    for (Index i{0}; i < getNumNodes(); ++i) {
      if (m_mesh.isActiveNode(i)) {
        momentum_n[i] += forceTotal[i] * m_dt;
        acceleration_n[i] = forceTotal[i] / mass_n[i];
        // velocity_n[i] += acceleration_n[i] * m_dt;
        // applyNodalVeloConstraint();
        // position_n[i] += velocity_n[i] * m_dt; //Grid nodes don't change
        // position?
      }
    }
    applyNodalAccConstraint(); // No need to apply AccConstraint here as Force
                               // Constraint is already applied
    applyFrictionalBC();

    // Update particle position and velocity

    // Map back to MPs
    for (Index p{0}; p < m_mesh.getNumMPs(); ++p) {
      Index e = m_mesh.getMPelementID(p);
      if (e != -1) {
        T x_p = position_p[p];
        Index n1 = m_mesh.getEleConnectivity(e)[0];
        Index n2 = m_mesh.getEleConnectivity(e)[1];
        auto [x1, x2] = m_mesh.getElementNodes(e);
        T xi = parentCoor(x_p, x1, x2);
        // Update velocity and position using FLIP style
        velocity_p[p] += (N1_ref(xi) * acceleration_n[n1] +
                          N2_ref(xi) * acceleration_n[n2]) *
                         m_dt;
        // Hybrid
        // T v_pic = N1 * velocity_n[n1] + N2 * velocity_n[n2];
        // T v_flip = velocity_p[p] + (N1 * a_n1 + N2 * a_n2) * dt;
        // velocity_p[p] = alpha * v_pic + '(1-alpha)' * v_flip;
        // position_p[p] += velocity_p[p] * m_dt;
        position_p[p] += (N1_ref(xi) * momentum_n[n1] / mass_n[n1] +
                          N2_ref(xi) * momentum_n[n2] / mass_n[n2]) *
                         m_dt;
        momentum_p[p] = mass_p[p] * velocity_p[p];
      }
    }

    for (Index p{0}; p < m_mesh.getNumMPs(); ++p) {
      Index e = m_mesh.getMPelementID(p);
      if (e != -1) {
        T x_p = position_p[p];
        Index n1 = m_mesh.getEleConnectivity(e)[0];
        Index n2 = m_mesh.getEleConnectivity(e)[1];
        auto [x1, x2] = m_mesh.getElementNodes(e);
        T xi = parentCoor(x_p, x1, x2);
        velocity_n[n1] += momentum_p[p] * N1_ref(xi) / mass_n[n1];
        velocity_n[n2] += momentum_p[p] * N2_ref(xi) / mass_n[n2];
      }
    }
    applyNodalVeloConstraint();

    for (Index p{0}; p < m_mesh.getNumMPs(); ++p) {
      Index e = m_mesh.getMPelementID(p);
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
        // Update volume
        volume_p[p] *= (1.0 + dStrain_p[p]);
      }
    }
  }

  void resetMesh() {
    m_mesh.setMPCoords(position_p); // Saving the updated MPs' position to
    Mesh

        mass_n.resetZero();
    momentum_n.resetZero();
    velocity_n.resetZero();
    acceleration_n.resetZero();

    bodyForce_n.resetZero();
    tractionForce_n.resetZero();
    forceExternal_n.resetZero();
    forceInternal_n.resetZero();
    forceTotal.resetZero();

    for (auto s : stress_n) {
      s = Matrix<T, 2, 2>::zero();
    }

    m_mesh.nodalReset();
  }

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

  void exportResult(const std::string &filename = "mpm1D_results.txt") {}
  void applyBC() {}
  void timeIntegration() {}
};

#endif // MATERIAL_POINT_METHOD_H
