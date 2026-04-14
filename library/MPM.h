#ifndef MATERIAL_POINT_METHOD_H
#define MATERIAL_POINT_METHOD_H

#include "Matrix.h"
#include "Mesh.h"
#include "Vector.h"
#include "elasticity.h"
#include "gaussQuadrature.h"
#include "parentElement.h"
#include "physicConstants.h"
#include "signFunction.h"
#include <cassert>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <system_error>
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
  DynamicVector<T> n_mass{};
  // DynamicVector<T> position_n{}; // Don't get used
  DynamicVector<T> n_velocity{};
  DynamicVector<T> n_acceleration{};
  DynamicVector<T> n_momentum{};
  DynamicVector<char> n_velocityConstrained{};
  DynamicVector<char> n_accelerationConstrained{};
  DynamicVector<char> n_momentumConstrained{};
  DynamicVector<char> n_forceConstrained{};

  // Constraint values (so constraints can be set anytime and enforced later)
  DynamicVector<T> n_velocityConstraintValue{};
  DynamicVector<T> n_accelerationConstraintValue{};
  DynamicVector<T> n_momentumConstraintValue{};
  DynamicVector<T> n_forceConstraintValue{};
  DynamicVector<T> n_bodyForce{}, n_tractionForce{};
  // Nodal external forces
  DynamicVector<T> n_forceExternal{}, n_forceInternal{}, n_forceTotal{};

  // Material Points p
  DynamicVector<Particle1D<T>> p_particles{};

  static bool endsWith(const std::string &s, const std::string &suffix) {
    return s.size() >= suffix.size() &&
           s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
  }

  static std::string withExtensionReplaced(const std::string &filename,
                                           const std::string &fromExt,
                                           const std::string &toExt) {
    if (endsWith(filename, fromExt))
      return filename.substr(0, filename.size() - fromExt.size()) + toExt;
    return filename + toExt;
  }

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
    n_mass.resize(nNodes, T{});
    // position_n.resize(nNodes, T{});
    n_velocity.resize(nNodes, T{});
    n_acceleration.resize(nNodes, T{});
    n_momentum.resize(nNodes, T{});

    n_velocityConstrained.resize(nNodes, 0);
    n_accelerationConstrained.resize(nNodes, 0);
    n_momentumConstrained.resize(nNodes, 0);
    n_forceConstrained.resize(nNodes, 0);
    n_velocityConstraintValue.resize(nNodes, T{});
    n_accelerationConstraintValue.resize(nNodes, T{});
    n_momentumConstraintValue.resize(nNodes, T{});
    n_forceConstraintValue.resize(nNodes, T{});

    n_bodyForce.resize(nNodes, T{});
    n_tractionForce.resize(nNodes, T{});
    n_forceExternal.resize(nNodes, T{});
    n_forceInternal.resize(nNodes, T{});
    n_forceTotal.resize(nNodes, T{});

    // Initialize MP vectors
    Index nMPs = m_mesh.getNumMPs();
    p_particles = m_mesh.getMPs();
    for (Index p{0}; p < nMPs; ++p) {
      p_particles[p].volume = m_volume / nMPs;
      p_particles[p].mass = m_mass / nMPs; // MP's mass is constant
      p_particles[p].vel = m_v0;
      p_particles[p].acc = T{};
      p_particles[p].force = T{};
      p_particles[p].momentum = T{};
      p_particles[p].stress = T{};
      p_particles[p].strain = T{};
      p_particles[p].strainRate = T{};
      p_particles[p].dStrain = T{};
    }
    m_mesh.setMPs(p_particles);

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

  T getNodalMass(Index i) const { return n_mass[i]; }
  T getNodalVelocity(Index i) const { return n_velocity[i]; }
  T getNodalMomentum(Index i) const { return n_momentum[i]; }
  T getNodalExtForce(Index i) const { return n_forceExternal[i]; }
  T getNodalIntForce(Index i) const { return n_forceInternal[i]; }
  T getNodalTotalForce(Index i) const { return n_forceTotal[i]; }

  T getMPvolume(Index p) const { return p_particles[p].volume; }
  T getMPmass(Index p) const { return p_particles[p].mass; }
  T getMPvelocity(Index p) const { return p_particles[p].vel; }
  T getMPposition(Index p) const { return p_particles[p].pos; }
  T getMPmomentum(Index p) const { return p_particles[p].momentum; }
  T getMPstrain(Index p) const { return p_particles[p].strain; }
  T getMPstrainRate(Index p) const { return p_particles[p].strainRate; }
  T getMPdStrain(Index p) const { return p_particles[p].dStrain; }
  T getMPstress(Index p) const { return p_particles[p].stress; }

  // Setters
  void setCurrentTime(T value) { m_currentTime = value; }
  void setE(T E) { m_E = E; }
  void setG(T G) { m_G = G; } // Set G if considering gravity
  void setAnalyticSolution(std::function<T(T)> sol) {
    m_analyticSolution = sol;
  }
  void setComportmentLaw(std::function<T(T)> law) { m_law = law; }
  void setMPvelocity(Index p, T value) { p_particles[p].vel = value; }

  // Setters: store constraint values + mark constrained
  void setNodalVeloConstraint(Index i, T value) {
    n_velocityConstraintValue[i] = value;
    n_velocityConstrained[i] = 1;
  }
  void setNodalAccConstraint(Index i, T value) {
    n_accelerationConstraintValue[i] = value;
    n_accelerationConstrained[i] = 1;
  }
  void setNodalMomentumConstraint(Index i, T value) {
    n_momentumConstraintValue[i] = value;
    n_momentumConstrained[i] = 1;
  }
  void setNodalForceConstraint(Index i, T value) {
    n_forceConstraintValue[i] = value;
    n_forceConstrained[i] = 1;
  }

  // Apply stored constraints to current nodal state
  void applyNodalVeloConstraint() {
    for (Index i{0}; i < getNumNodes(); ++i) {
      if (n_velocityConstrained[i] != 0) {
        n_velocity[i] = n_velocityConstraintValue[i];
      }
    }
  }
  void applyNodalAccConstraint() {
    for (Index i{0}; i < getNumNodes(); ++i) {
      if (n_accelerationConstrained[i] != 0) {
        n_acceleration[i] = n_accelerationConstraintValue[i];
      }
    }
  }
  void applyNodalMomentumConstraint() {
    for (Index i{0}; i < getNumNodes(); ++i) {
      if (n_momentumConstrained[i] != 0) {
        n_momentum[i] = n_momentumConstraintValue[i];
      }
    }
  }
  void applyNodalForceConstraint() {
    for (Index i{0}; i < getNumNodes(); ++i) {
      if (n_forceConstrained[i] != 0) {
        n_forceTotal[i] = n_forceConstraintValue[i];
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
      p_particles[p].momentum = p_particles[p].mass * p_particles[p].vel;
      Index e = m_mesh.getMPelementID(p);
      if (e != -1) {
        T x_p = p_particles[p].pos;
        Index n1 = m_mesh.getEleConnectivity(e)[0];
        Index n2 = m_mesh.getEleConnectivity(e)[1];
        auto [x1, x2] = m_mesh.getElementNodes(e);
        T xi = parentCoor(x_p, x1, x2);
        n_mass[n1] += N1_ref(xi) * p_particles[p].mass;
        n_mass[n2] += N2_ref(xi) * p_particles[p].mass;
        n_momentum[n1] += N1_ref(xi) * p_particles[p].momentum;
        n_momentum[n2] += N2_ref(xi) * p_particles[p].momentum;
      }
    }

    applyNodalMomentumConstraint();
  }

  void nodalEquilibrium() {
    // f^ext = b + t
    for (Index p{0}; p < m_mesh.getNumMPs(); ++p) {
      Index e = m_mesh.getMPelementID(p);
      if (e != -1) {
        T x_p = p_particles[p].pos;
        Index n1 = m_mesh.getEleConnectivity(e)[0];
        Index n2 = m_mesh.getEleConnectivity(e)[1];
        auto [x1, x2] = m_mesh.getElementNodes(e);
        T xi = parentCoor(x_p, x1, x2);
        n_bodyForce[n1] += m_G * N1_ref(xi) * p_particles[p].mass;
        n_bodyForce[n2] += m_G * N2_ref(xi) * p_particles[p].mass;
        // Traction force t_i (to be implemented)
        n_forceInternal[n1] -=
            p_particles[p].volume * dN1_dx(x1, x2) * p_particles[p].stress;
        n_forceInternal[n2] -=
            p_particles[p].volume * dN2_dx(x1, x2) * p_particles[p].stress;
      }
    }
    n_forceExternal = n_bodyForce + n_tractionForce;
    n_forceTotal = n_forceExternal + n_forceInternal;

    // Enforce any stored nodal force constraints after assembly
    applyNodalForceConstraint();
  }

  void n2p() {
    // Update momentum at nodes
    for (Index i{0}; i < getNumNodes(); ++i) {
      if (m_mesh.isActiveNode(i)) {
        n_momentum[i] += n_forceTotal[i] * m_dt;
        n_acceleration[i] = n_forceTotal[i] / n_mass[i];

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
        T x_p = p_particles[p].pos;
        Index n1 = m_mesh.getEleConnectivity(e)[0];
        Index n2 = m_mesh.getEleConnectivity(e)[1];
        auto [x1, x2] = m_mesh.getElementNodes(e);
        T xi = parentCoor(x_p, x1, x2);
        // Update velocity and position using FLIP style
        p_particles[p].acc =
            N1_ref(xi) * n_acceleration[n1] + N2_ref(xi) * n_acceleration[n2];
        p_particles[p].vel += p_particles[p].acc * m_dt;
        // Hybrid
        // T v_pic = N1 * velocity_n[n1] + N2 * velocity_n[n2];
        // T v_flip = velocity_p[p] + (N1 * a_n1 + N2 * a_n2) * dt;
        // velocity_p[p] = alpha * v_pic + '(1-alpha)' * v_flip;
        // position_p[p] += velocity_p[p] * m_dt;
        p_particles[p].pos += (N1_ref(xi) * n_momentum[n1] / n_mass[n1] +
                               N2_ref(xi) * n_momentum[n2] / n_mass[n2]) *
                              m_dt;
        p_particles[p].momentum = p_particles[p].mass * p_particles[p].vel;
      }
    }

    for (Index p{0}; p < m_mesh.getNumMPs(); ++p) {
      Index e = m_mesh.getMPelementID(p);
      if (e != -1) {
        T x_p = p_particles[p].pos;
        Index n1 = m_mesh.getEleConnectivity(e)[0];
        Index n2 = m_mesh.getEleConnectivity(e)[1];
        auto [x1, x2] = m_mesh.getElementNodes(e);
        T xi = parentCoor(x_p, x1, x2);
        n_velocity[n1] += p_particles[p].momentum * N1_ref(xi) / n_mass[n1];
        n_velocity[n2] += p_particles[p].momentum * N2_ref(xi) / n_mass[n2];
      }
    }
    applyNodalVeloConstraint();

    for (Index p{0}; p < m_mesh.getNumMPs(); ++p) {
      Index e = m_mesh.getMPelementID(p);
      if (e != -1) {
        T x_p = p_particles[p].pos;
        Index n1 = m_mesh.getEleConnectivity(e)[0];
        Index n2 = m_mesh.getEleConnectivity(e)[1];
        auto [x1, x2] = m_mesh.getElementNodes(e);
        T xi = parentCoor(x_p, x1, x2);
        // Attention: x1-x2 belong to nodes (their positions don't get
        // updated), while the velocity is measured at MPs (updated at t+dt)
        p_particles[p].strainRate =
            dN1_dx(x1, x2) * n_velocity[n1] + dN2_dx(x1, x2) * n_velocity[n2];
        p_particles[p].dStrain = p_particles[p].strainRate * m_dt;
        p_particles[p].strain += p_particles[p].dStrain;

        // Constitutive law:
        if (m_law) {
          p_particles[p].stress += m_law(p_particles[p].dStrain);
        } else {
          p_particles[p].stress +=
              m_E * p_particles[p].dStrain; // Default linear elastic
        }
        // Update volume
        p_particles[p].volume *= (T{1} + p_particles[p].dStrain);
      }
    }
  }

  void resetMesh() {
    m_mesh.setMPs(p_particles); // Saving updated MP state back to mesh

    n_mass.resetZero();
    n_momentum.resetZero();
    n_velocity.resetZero();
    n_acceleration.resetZero();

    n_bodyForce.resetZero();
    n_tractionForce.resetZero();
    n_forceExternal.resetZero();
    n_forceInternal.resetZero();
    n_forceTotal.resetZero();

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
    T min_dist = std::abs(p_particles[0].pos - m_xloc);
    Index closest_p = 0;
    for (Index p = 1; p < m_mesh.getNumMPs(); ++p) {
      T dist = std::abs(p_particles[p].pos - m_xloc);
      if (dist < min_dist) {
        min_dist = dist;
        closest_p = p;
      }
    }
    T x_numerical = p_particles[closest_p].pos;

    T x_exact = m_analyticSolution(m_currentTime);
    T error = std::abs(x_numerical - x_exact);

    std::cout << std::setw(10) << std::fixed << std::setprecision(4)
              << m_currentTime << std::setw(15) << std::setprecision(6)
              << x_numerical << std::setw(15) << x_exact << std::setw(15)
              << std::scientific << error << '\n';
  }

  void exportResult(const std::string &filename = "mpm1D_results.vtk") {
    const std::string vtkFile =
        endsWith(filename, ".vtk")
            ? filename
            : withExtensionReplaced(filename, ".txt", ".vtk");

    // Also export the background grid mesh (nodes + elements) to a sibling
    // file. Default naming:
    //   - if vtkFile contains "particles_" => replace it by "mesh_"
    //   - else => append "_mesh" before the extension
    std::string meshFile = vtkFile;
    if (const std::size_t pos = meshFile.rfind("particles_");
        pos != std::string::npos) {
      meshFile.replace(pos, std::string("particles_").size(), "mesh_");
    } else if (endsWith(meshFile, ".vtk")) {
      meshFile = meshFile.substr(0, meshFile.size() - 4) + "_mesh.vtk";
    } else {
      meshFile += "_mesh.vtk";
    }

    std::ofstream vtk(vtkFile);
    if (!vtk)
      throw std::runtime_error("MPM1D::exportResult: cannot open file: " +
                               vtkFile);

    const Index nPoints = p_particles.size();

    vtk << "# vtk DataFile Version 3.0\n";
    vtk << "MPM1D particles\n";
    vtk << "ASCII\n";
    vtk << "DATASET POLYDATA\n";
    vtk << "POINTS " << nPoints << " double\n";
    for (Index p = 0; p < nPoints; ++p) {
      vtk << p_particles[p].pos << " 0 0\n";
    }

    vtk << "VERTICES " << nPoints << " " << (nPoints * 2) << "\n";
    for (Index p = 0; p < nPoints; ++p) {
      vtk << "1 " << p << "\n";
    }

    vtk << "POINT_DATA " << nPoints << "\n";
    vtk << "VECTORS velocity double\n";
    for (Index p = 0; p < nPoints; ++p) {
      vtk << p_particles[p].vel << " 0 0\n";
    }

    vtk << "SCALARS mass double 1\n";
    vtk << "LOOKUP_TABLE default\n";
    for (Index p = 0; p < nPoints; ++p) {
      vtk << p_particles[p].mass << "\n";
    }

    vtk << "SCALARS volume double 1\n";
    vtk << "LOOKUP_TABLE default\n";
    for (Index p = 0; p < nPoints; ++p) {
      vtk << p_particles[p].volume << "\n";
    }

    // ---- Mesh grid export (VTK legacy, UNSTRUCTURED_GRID) ----
    std::ofstream mesh(meshFile);
    if (!mesh)
      throw std::runtime_error("MPM1D::exportResult: cannot open mesh file: " +
                               meshFile);

    const Index numNodes = m_mesh.getNumNodes();
    const Index numElems = m_mesh.getNumElements();

    mesh << "# vtk DataFile Version 3.0\n";
    mesh << "MPM1D background grid\n";
    mesh << "ASCII\n";
    mesh << "DATASET UNSTRUCTURED_GRID\n";
    mesh << "POINTS " << numNodes << " double\n";
    for (Index i = 0; i < numNodes; ++i) {
      mesh << m_mesh.nodeCoords()[i] << " 0 0\n";
    }

    // Each line cell: "2 n1 n2" => 3 integers per element
    mesh << "CELLS " << numElems << " " << (numElems * 3) << "\n";
    for (Index e = 0; e < numElems; ++e) {
      const auto &conn = m_mesh.getEleConnectivity(e);
      mesh << "2 " << conn[0] << " " << conn[1] << "\n";
    }

    // VTK cell type for line = 3
    mesh << "CELL_TYPES " << numElems << "\n";
    for (Index e = 0; e < numElems; ++e) {
      mesh << "3\n";
    }

    mesh << "POINT_DATA " << numNodes << "\n";

    mesh << "SCALARS active int 1\n";
    mesh << "LOOKUP_TABLE default\n";
    const auto &active = m_mesh.getActiveNodes();
    for (Index i = 0; i < numNodes; ++i) {
      const int a = (i < active.size() && active[i] != 0) ? 1 : 0;
      mesh << a << "\n";
    }

    mesh << "SCALARS mass double 1\n";
    mesh << "LOOKUP_TABLE default\n";
    for (Index i = 0; i < numNodes; ++i) {
      const T m = (i < n_mass.size()) ? n_mass[i] : T{0};
      mesh << m << "\n";
    }

    mesh << "VECTORS velocity double\n";
    for (Index i = 0; i < numNodes; ++i) {
      const T v = (i < n_velocity.size()) ? n_velocity[i] : T{0};
      mesh << v << " 0 0\n";
    }

    mesh << "VECTORS acceleration double\n";
    for (Index i = 0; i < numNodes; ++i) {
      const T a = (i < n_acceleration.size()) ? n_acceleration[i] : T{0};
      mesh << a << " 0 0\n";
    }

    mesh << "VECTORS force_total double\n";
    for (Index i = 0; i < numNodes; ++i) {
      const T f = (i < n_forceTotal.size()) ? n_forceTotal[i] : T{0};
      mesh << f << " 0 0\n";
    }

    mesh << "VECTORS momentum double\n";
    for (Index i = 0; i < numNodes; ++i) {
      const T p = (i < n_momentum.size()) ? n_momentum[i] : T{0};
      mesh << p << " 0 0\n";
    }
  }

  void exportVTKFrame(const std::filesystem::path &outputDir, Index frame,
                      int padWidth = 6) {
    std::error_code ec;
    std::filesystem::create_directories(outputDir, ec);

    std::ostringstream name;
    name << "particles_" << std::setw(padWidth) << std::setfill('0')
         << static_cast<long long>(frame) << ".vtk";
    exportResult((outputDir / name.str()).string());
  }

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
  std::pair<T, T> m_minCorner{};
  std::pair<T, T> m_maxCorner{};
  T m_pLength{}, m_pHeight{};
  T m_volume{1.0}, m_volume0{1.0}; // Volume
  T m_mass{1000};                  // mass
  T m_G{0.0};                      // Gravity Acceleration
  DynamicVector<T> m_v0{T{}, T{}}; // Initial velocity

  // Simulation properties
  T m_currentTime{0.0}; // Current time
  T m_dt{0.0};          // Time step
  T m_duration{10.0};   // Duration of simulation
  Index m_nSteps{0};    // Number of steps
  // Index m_interval{10}; // Output interval

  // Behavior law (stress increment from strain increment)
  std::function<Matrix<T, 2, 2>(const Matrix<T, 2, 2> &)> m_law;

  //  Mesh
  Mesh2D<T> m_mesh{};

  // Nodes n
  DynamicVector<T> n_mass{};
  // DynamicVector<T> position_n{}; // Don't get used
  DynamicVector<DynamicVector<T>> n_velocity{};
  DynamicVector<DynamicVector<T>> n_acceleration{};
  DynamicVector<DynamicVector<T>> n_momentum{};
  DynamicVector<DynamicVector<T>> n_displacement{}; // Displacement
  DynamicVector<Matrix<T, 2, 2>> n_stress{};        // Nodal stress (optional)

  DynamicVector<DynamicVector<char>> n_velocityConstrained{};
  DynamicVector<DynamicVector<char>> n_accelerationConstrained{};
  DynamicVector<DynamicVector<char>> n_momentumConstrained{};
  DynamicVector<DynamicVector<char>> n_forceConstrained{};

  // Constraint values (so constraints can be set anytime and enforced later)
  DynamicVector<DynamicVector<T>> n_velocityConstraintValue{};
  DynamicVector<DynamicVector<T>> n_accelerationConstraintValue{};
  DynamicVector<DynamicVector<T>> n_momentumConstraintValue{};
  DynamicVector<DynamicVector<T>> n_forceConstraintValue{};
  DynamicVector<DynamicVector<T>> n_bodyForce{}, n_tractionForce{};
  // Nodal external forces
  DynamicVector<DynamicVector<T>> n_forceExternal{}, n_forceInternal{},
      n_forceTotal{};

  // Material Points p
  DynamicVector<Particle2D<T>> p_particles{};
  DynamicVector<T> p_volume0{}; // Initial volume (for large-strain update)
  DynamicVector<DynamicVector<T>> p_momentum{}; // Momentum

  DynamicVector<Matrix<T, 2, 2>> p_stress{}, p_strain{}, p_strainRate{},
      p_deformGradient{}, p_dStrain{}; // Stress + strain

  static constexpr Index dir_x = 0;
  static constexpr Index dir_y = 1;
  static constexpr Index directionPositive = 1;
  static constexpr Index directionNegative = -1;
  static constexpr Index dimensions = 2;

  static bool endsWith(const std::string &s, const std::string &suffix) {
    return s.size() >= suffix.size() &&
           s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
  }

  static std::string withExtensionReplaced(const std::string &filename,
                                           const std::string &fromExt,
                                           const std::string &toExt) {
    if (endsWith(filename, fromExt))
      return filename.substr(0, filename.size() - fromExt.size()) + toExt;
    return filename + toExt;
  }

public:
  // Constructor
  MPM2D(T E, T nu, T rho, T mu, T phi, T c, T K0,
        const std::pair<T, T> &minCorner, const std::pair<T, T> &maxCorner,
        T dt, T duration, T v0 = T{})
      : m_E{E}, m_nu{nu}, m_rho{rho}, m_mu{mu}, m_phi{phi}, m_c{c}, m_K0{K0},
        m_minCorner{minCorner}, m_maxCorner{maxCorner}, m_v0{v0, T{}}, m_dt{dt},
        m_duration{duration}, m_nSteps{static_cast<Index>(duration / dt)},
        m_pLength{maxCorner.first - minCorner.first},
        m_pHeight{maxCorner.second - minCorner.second},
        m_volume{constexpr_fabs(m_pLength * m_pHeight)}, m_mass{rho * m_volume},
        m_volume0{m_volume}, m_mesh{Mesh2D<T>{{gridLength, gridHeight},
                                              {nx, ny},
                                              minCorner,
                                              maxCorner,
                                              MP_size}} {

    // Check critical time
    T c_wave = std::sqrt(m_E / rho);
    T dt_crit = (m_pLength < m_pHeight ? m_pLength : m_pHeight) / c_wave;
    assert((dt_crit / T{10}) >= dt &&
           "dt doesn't satisfied CFL condition (too big)");

    const Index nNodes = m_mesh.getNumNodes();
    const Index nMPs = m_mesh.getNumMPs();

    // Initialize nodal vectors
    n_mass.resize(nNodes, T{});
    // position_n.resize(nNodes, DynamicVector<T>{});
    n_velocity.resize(nNodes, DynamicVector<T>{0, 0});
    n_acceleration.resize(nNodes, DynamicVector<T>{0, 0});
    n_momentum.resize(nNodes, DynamicVector<T>{0, 0});
    n_displacement.resize(nNodes, DynamicVector<T>{0, 0});
    n_stress.resize(nNodes, Matrix<T, 2, 2>::zero());

    n_velocityConstrained.resize(nNodes, DynamicVector<char>{0, 0});
    n_accelerationConstrained.resize(nNodes, DynamicVector<char>{0, 0});
    n_momentumConstrained.resize(nNodes, DynamicVector<char>{0, 0});
    n_forceConstrained.resize(nNodes, DynamicVector<char>{0, 0});
    n_velocityConstraintValue.resize(nNodes, DynamicVector<T>{0, 0});
    n_accelerationConstraintValue.resize(nNodes, DynamicVector<T>{0, 0});
    n_momentumConstraintValue.resize(nNodes, DynamicVector<T>{0, 0});
    n_forceConstraintValue.resize(nNodes, DynamicVector<T>{0, 0});

    n_bodyForce.resize(nNodes, DynamicVector<T>{0, 0});
    n_tractionForce.resize(nNodes, DynamicVector<T>{0, 0});
    n_forceExternal.resize(nNodes, DynamicVector<T>{0, 0});
    n_forceInternal.resize(nNodes, DynamicVector<T>{0, 0});
    n_forceTotal.resize(nNodes, DynamicVector<T>{0, 0});

    // Initialize MP vectors
    p_particles = m_mesh.getMPs();
    const T mpVolume = m_volume / nMPs;
    const T mpMass = m_mass / nMPs;
    for (Index p{0}; p < nMPs; ++p) {
      p_particles[p].volume = mpVolume;
      p_particles[p].mass = mpMass;
      p_particles[p].vel = m_v0;
      p_particles[p].acc.resetZero();
      p_particles[p].force.resetZero();
    }
    m_mesh.setMPs(p_particles);

    p_volume0.resize(nMPs, m_volume / nMPs);
    p_momentum.resize(nMPs, DynamicVector<T>{0, 0}); // Compute later

    p_stress.resize(nMPs, Matrix<T, 2, 2>::zero());
    p_strain.resize(nMPs, Matrix<T, 2, 2>::zero());
    p_strainRate.resize(nMPs, Matrix<T, 2, 2>::zero());
    p_deformGradient.resize(nMPs, Matrix<T, 2, 2>::identity());
    p_dStrain.resize(nMPs, Matrix<T, 2, 2>::zero());

    for (Index p{0}; p < nMPs; ++p) {
      m_totalEnergy0 += m_G * p_particles[p].pos.y() * p_particles[p].mass;
    }

    // m_mesh.print();
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

  T getNodalMass(Index i) const { return n_mass[i]; }
  T getNodalVelocity(Index i) const { return n_velocity[i].x(); }
  T getNodalMomentum(Index i) const { return n_momentum[i].x(); }
  T getNodalExtForce(Index i) const { return n_forceExternal[i].x(); }
  T getNodalIntForce(Index i) const { return n_forceInternal[i].x(); }
  T getNodalTotalForce(Index i) const { return n_forceTotal[i].x(); }

  T getMPvolume(Index p) const { return p_particles[p].volume; }
  T getMPmass(Index p) const { return p_particles[p].mass; }
  T getMPvelocity(Index p) const { return p_particles[p].vel.x(); }
  T getMPposition(Index p) const { return p_particles[p].pos.x(); }
  T getMPmomentum(Index p) const { return p_momentum[p].x(); }
  T getMPstrain(Index p) const { return p_strain[p].xx(); }
  T getMPstrainRate(Index p) const { return p_strainRate[p].xx(); }
  T getMPdStrain(Index p) const { return p_dStrain[p].xx(); }
  T getMPstress(Index p) const { return p_stress[p].xx(); }

  // Setters
  void setCurrentTime(T value) { m_currentTime = value; }
  void setE(T E) { m_E = E; }
  void setG(T G) { m_G = G; } // Set G if considering gravity

  void setComportmentLaw(
      const std::function<Matrix<T, 2, 2>(const Matrix<T, 2, 2> &)> &law) {
    m_law = law;
  }
  void setComportmentLaw(const std::function<T(T)> &law) {
    if (!law) {
      m_law = {};
      return;
    }
    m_law = [law](const Matrix<T, 2, 2> &dEps) {
      return Matrix<T, 2, 2>{law(dEps.xx()), law(dEps.xy()), law(dEps.yx()),
                             law(dEps.yy())};
    };
  }
  void setMPvelocity(Index p, const DynamicVector<T> &value) {
    p_particles[p].vel.x() = value.x();
    p_particles[p].vel.y() = value.y();
  }
  void setMPvelocity(Index p, T value) {
    p_particles[p].vel.x() = value;
    p_particles[p].vel.y() = T{};
  }

  // Setters: store constraint values + mark constrained
  void setNodalVeloConstraint(Index i, const DynamicVector<T> &value) {
    n_velocityConstraintValue[i] = value;
    n_velocityConstrained[i] = DynamicVector<char>{1, 1};
  }
  void setNodalVeloConstraint(Index i, T value) {
    setNodalVeloConstraint(i, DynamicVector<T>{value, value});
  }

  void setNodalAccConstraint(Index i, const DynamicVector<T> &value) {
    n_accelerationConstraintValue[i] = value;
    n_accelerationConstrained[i] = DynamicVector<char>{1, 1};
  }
  void setNodalAccConstraint(Index i, T value) {
    setNodalAccConstraint(i, DynamicVector<T>{value, value});
  }

  void setNodalMomentumConstraint(Index i, const DynamicVector<T> &value) {
    n_momentumConstraintValue[i] = value;
    n_momentumConstrained[i] = DynamicVector<char>{1, 1};
  }
  void setNodalMomentumConstraint(Index i, T value) {
    setNodalMomentumConstraint(i, DynamicVector<T>{value, value});
  }

  void setNodalForceConstraint(Index i, const DynamicVector<T> &value) {
    n_forceConstraintValue[i] = value;
    n_forceConstrained[i] = DynamicVector<char>{1, 1};
  }
  void setNodalForceConstraint(Index i, T value) {
    setNodalForceConstraint(i, DynamicVector<T>{value, value});
  }

  // Apply stored constraints to current nodal state
  void applyNodalVeloConstraint() {
    for (Index i{0}; i < getNumNodes(); ++i) {
      for (Index d{0}; d < dimensions; ++d) {
        if (n_velocityConstrained[i][d] != 0) {
          n_velocity[i][d] = n_velocityConstraintValue[i][d];
        }
      }
    }
  }
  void applyNodalAccConstraint() {
    for (Index i{0}; i < getNumNodes(); ++i) {
      for (Index d{0}; d < dimensions; ++d) {
        if (n_accelerationConstrained[i][d] != 0) {
          n_acceleration[i][d] = n_accelerationConstraintValue[i][d];
        }
      }
    }
  }
  void applyNodalMomentumConstraint() {
    for (Index i{0}; i < getNumNodes(); ++i) {
      for (Index d{0}; d < dimensions; ++d) {
        if (n_momentumConstrained[i][d] != 0) {
          n_momentum[i][d] = n_momentumConstraintValue[i][d];
        }
      }
    }
  }
  void applyNodalForceConstraint() {
    for (Index i{0}; i < getNumNodes(); ++i) {
      for (Index d{0}; d < dimensions; ++d) {
        if (n_forceConstrained[i][d] != 0) {
          n_forceTotal[i][d] = n_forceConstraintValue[i][d];
        }
      }
    }
  }

  // Source: https://www.geoelements.org/LearnMPM/mpm2d-column-collapse.html
  void frictionalBC(const DynamicVector<Index> &nodeIDs, Index dir_n,
                    Index signDir_n) {
    Index dir_t = 1 - dir_n; // tangential diretion
    // Normal and tangential acceleration
    const Index nNodes = nodeIDs.size();
    DynamicVector<T> acc_n(nNodes), acc_t(nNodes), vel_t(nNodes);
    DynamicVector<char> moveTowardBoundary(nNodes);
    for (Index i{0}; i < nNodes; ++i) {
      const Index nodeID = nodeIDs[i];
      acc_n[i] = n_acceleration[nodeID][dir_n];
      acc_t[i] = n_acceleration[nodeID][dir_t];
      moveTowardBoundary[i] =
          (acc_n[i] * static_cast<T>(signDir_n)) > T{0} ? 1 : 0;
      vel_t[i] = n_momentum[nodeID][dir_t] / n_mass[nodeID];
    }

    // Apply frictional boundary condition
    for (Index i{0}; i < nNodes; ++i) {
      if (moveTowardBoundary[i] == 1) {
        // Determine static or kinetic friction
        if (!approximatelyEqualAbsRel(vel_t[i], T{0})) // kinetic friction
        {
          // Compute tangential velocity at next time step
          const T vel_net = m_dt * acc_t[i] + vel_t[i];
          const T vel_frictional = m_dt * m_mu * constexpr_fabs(acc_n[i]);
          if (constexpr_fabs(vel_net) <= vel_frictional) {
            // friction stops the particle
            acc_t[i] = -vel_t[i] / m_dt; // vel_net = 0
          } else {
            // friction reduces the tangential acceleration
            acc_t[i] -= sgn(vel_net) * m_mu * constexpr_fabs(acc_n[i]);
          }
        } else // static friction
        {
          if (constexpr_fabs(acc_t[i]) <= m_mu * constexpr_fabs(acc_n[i])) {
            acc_t[i] = T{0};
          } else {
            acc_t[i] -= sgn(acc_t[i]) * m_mu * constexpr_fabs(acc_n[i]);
          }
        }
      }
      const Index nodeID = nodeIDs[i];
      // Update tangential acceleration
      n_acceleration[nodeID][dir_t] = acc_t[i];
      // Update nodal force at this node
      n_forceTotal[nodeID] = n_acceleration[nodeID] * n_mass[nodeID];
    }
  }

  void applyFrictionalBC() {
    // Update acceleration and force on boundary
    frictionalBC(m_mesh.leftActiveNodes(), dir_x, directionNegative);
    frictionalBC(m_mesh.rightActiveNodes(), dir_x, directionPositive);
    frictionalBC(m_mesh.bottomActiveNodes(), dir_y, directionNegative);
    frictionalBC(m_mesh.topActiveNodes(), dir_y, directionPositive);

    // Constraint momentum and force on boundary
    for (Index i : m_mesh.leftActiveNodes()) {
      n_momentum[i].x() = T{};
      n_forceTotal[i].x() = T{};
    }
    for (Index i : m_mesh.rightActiveNodes()) {
      n_momentum[i].x() = T{};
      n_forceTotal[i].x() = T{};
    }
    for (Index i : m_mesh.bottomActiveNodes()) {
      n_momentum[i].y() = T{};
      n_forceTotal[i].y() = T{};
    }
    for (Index i : m_mesh.topActiveNodes()) {
      n_momentum[i].y() = T{};
      n_forceTotal[i].y() = T{};
    }
  }

  void initializeStress() {
    T ymax = m_pHeight - MP_size / 2;
    for (Index p{0}; p < this->getNumMps(); ++p) {
      p_stress[p].yy() = -m_G * m_rho * (ymax - p_particles[p].pos.y());
      p_stress[p].xx() = m_K0 * p_stress[p].yy();
    }
  }

  void computeEnergy() {
    m_potentEnergy = T{};
    m_kinEnergy = T{};
    for (Index p{0}; p < m_mesh.getNumMPs(); ++p) {
      m_potentEnergy += m_G * p_particles[p].mass * p_particles[p].pos.y();
      m_kinEnergy += T{0.5} * p_particles[p].mass *
                     dotProduct(p_particles[p].vel, p_particles[p].vel);
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
      p_momentum[p] =
          DynamicVector<T>{p_particles[p].mass * p_particles[p].vel.x(),
                           p_particles[p].mass * p_particles[p].vel.y()};
    }

    // Map nodal mass + momentum
    for (Index e{0}; e < m_mesh.getNumElements(); ++e) {
      if (!m_mesh.isActiveElement(e)) {
        continue;
      }

      const auto &mps = m_mesh.getMPsInElement(e);
      if (mps.size() == 0) {
        continue;
      }

      const auto &conn = m_mesh.getEleConnectivity(e);
      const Index n1 = conn[0];
      const Index n2 = conn[1];
      const Index n3 = conn[2];
      const Index n4 = conn[3];

      const auto [x_nodes, y_nodes] = m_mesh.getElementNodes(e);

      for (const Index p : mps) {
        const T x_p = p_particles[p].pos.x();
        const T y_p = p_particles[p].pos.y();
        const auto [xi, eta] = parentCoor(x_p, y_p, x_nodes, y_nodes);

        n_mass[n1] += N1_ref(xi, eta) * p_particles[p].mass;
        n_mass[n2] += N2_ref(xi, eta) * p_particles[p].mass;
        n_mass[n3] += N3_ref(xi, eta) * p_particles[p].mass;
        n_mass[n4] += N4_ref(xi, eta) * p_particles[p].mass;

        n_momentum[n1] += N1_ref(xi, eta) * p_momentum[p];
        n_momentum[n2] += N2_ref(xi, eta) * p_momentum[p];
        n_momentum[n3] += N3_ref(xi, eta) * p_momentum[p];
        n_momentum[n4] += N4_ref(xi, eta) * p_momentum[p];
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
      if (mps.size() == 0) {
        continue;
      }

      const auto &conn = m_mesh.getEleConnectivity(e);
      const Index n1 = conn[0];
      const Index n2 = conn[1];
      const Index n3 = conn[2];
      const Index n4 = conn[3];

      const auto [x_nodes, y_nodes] = m_mesh.getElementNodes(e);

      for (const Index p : mps) {
        const T x_p = p_particles[p].pos.x();
        const T y_p = p_particles[p].pos.y();
        const auto [xi, eta] = parentCoor(x_p, y_p, x_nodes, y_nodes);

        // Gravity body force (y-direction) : G should be negative!
        n_bodyForce[n1].y() += N1_ref(xi, eta) * m_G * p_particles[p].mass;
        n_bodyForce[n2].y() += N2_ref(xi, eta) * m_G * p_particles[p].mass;
        n_bodyForce[n3].y() += N3_ref(xi, eta) * m_G * p_particles[p].mass;
        n_bodyForce[n4].y() += N4_ref(xi, eta) * m_G * p_particles[p].mass;

        // Traction force t_i (to be implemented)

        // Internal force (both x and y)
        const auto [dN_dx, dN_dy] = gradientN(xi, eta, x_nodes, y_nodes);
        n_forceInternal[n1] -=
            p_particles[p].volume *
            (p_stress[p] * DynamicVector<T>{dN_dx[0], dN_dy[0]});
        n_forceInternal[n2] -=
            p_particles[p].volume *
            (p_stress[p] * DynamicVector<T>{dN_dx[1], dN_dy[1]});
        n_forceInternal[n3] -=
            p_particles[p].volume *
            (p_stress[p] * DynamicVector<T>{dN_dx[2], dN_dy[2]});
        n_forceInternal[n4] -=
            p_particles[p].volume *
            (p_stress[p] * DynamicVector<T>{dN_dx[3], dN_dy[3]});
      }
    };

    for (Index i{0}; i < m_mesh.getNumNodes(); ++i) {
      if (m_mesh.isActiveNode(i)) {
        n_forceExternal[i] = n_bodyForce[i] + n_tractionForce[i];
        n_forceTotal[i] = n_forceExternal[i] + n_forceInternal[i];
      }
    }

    // Enforce any stored nodal force constraints after assembly
    applyNodalForceConstraint();

    // Compute nodal acceleration from assembled forces
    for (Index i{0}; i < getNumNodes(); ++i) {
      if (m_mesh.isActiveNode(i) && !approximatelyEqualAbsRel(n_mass[i], T{})) {
        n_acceleration[i] = n_forceTotal[i] / n_mass[i];
      }
    }

    // Acc constraints first, then friction updates (acc + force)
    applyNodalAccConstraint();
    applyFrictionalBC();

    // Update momentum at nodes (after frictional BC)
    for (Index i{0}; i < getNumNodes(); ++i) {
      if (m_mesh.isActiveNode(i)) {
        n_momentum[i] += n_forceTotal[i] * m_dt;
      }
    }
    applyNodalMomentumConstraint();

    // Update nodal velocity from momentum and enforce velocity constraints
    for (Index i{0}; i < getNumNodes(); ++i) {
      if (m_mesh.isActiveNode(i) && !approximatelyEqualAbsRel(n_mass[i], T{})) {
        n_velocity[i] = n_momentum[i] / n_mass[i];
      }
    }
    applyNodalVeloConstraint();
  }

  void n2p() {
    // Map back to MPs
    for (Index p{0}; p < m_mesh.getNumMPs(); ++p) {
      const Index e = m_mesh.getMPelementID(p);
      if (e == -1) {
        continue;
      }

      const auto &conn = m_mesh.getEleConnectivity(e);
      const Index n1 = conn[0];
      const Index n2 = conn[1];
      const Index n3 = conn[2];
      const Index n4 = conn[3];

      const auto [x_nodes, y_nodes] = m_mesh.getElementNodes(e);
      const T x_p = p_particles[p].pos.x();
      const T y_p = p_particles[p].pos.y();
      const auto [xi, eta] = parentCoor(x_p, y_p, x_nodes, y_nodes);

      const T N1 = N1_ref(xi, eta);
      const T N2 = N2_ref(xi, eta);
      const T N3 = N3_ref(xi, eta);
      const T N4 = N4_ref(xi, eta);

      // Update particles' velocity , position and momentum
      const DynamicVector<T> a_next = N1_ref(xi, eta) * n_acceleration[n1] +
                                      N2_ref(xi, eta) * n_acceleration[n2] +
                                      N3_ref(xi, eta) * n_acceleration[n3] +
                                      N4_ref(xi, eta) * n_acceleration[n4];

      const DynamicVector<T> v_next =
          N1_ref(xi, eta) * n_velocity[n1] + N2_ref(xi, eta) * n_velocity[n2] +
          N3_ref(xi, eta) * n_velocity[n3] + N4_ref(xi, eta) * n_velocity[n4];

      p_particles[p].acc = a_next;
      p_particles[p].vel += StaticVector<T, 2>{a_next * m_dt};
      p_particles[p].pos += StaticVector<T, 2>{v_next * m_dt};
      p_momentum[p] =
          DynamicVector<T>{p_particles[p].mass * p_particles[p].vel.x(),
                           p_particles[p].mass * p_particles[p].vel.y()};

      // Update stress and strain
      const auto [dN_dx, dN_dy] = gradientN(xi, eta, x_nodes, y_nodes);
      Matrix<T, 2, 2> L = Matrix<T, 2, 2>::zero();

      for (Index a{0}; a < 4; ++a) {
        const Index nodeID = conn[a];
        L += tensorProduct<T, 2, 2>(n_velocity[nodeID],
                                    DynamicVector<T>{dN_dx[a], dN_dy[a]});
      }
      // // Strain rate from nodal velocities (small strain)
      // p_strainRate[p] = Matrix<T, 2, 2>::zero();
      // p_strainRate[p].xx() = L.xx();
      // p_strainRate[p].yy() = L.yy();
      // const T shear = T{0.5} * (L.xy() + L.yx());
      // p_strainRate[p].xy() = shear;
      // p_strainRate[p].yx() = shear;

      // p_dStrain[p] = p_strainRate[p] * m_dt;
      // p_strain[p] += p_dStrain[p];

      // if (m_law) {
      //   p_stress[p] += m_law(p_dStrain[p]);
      // } else {
      //   p_stress[p] += m_E * p_dStrain[p];
      // }

      // // 2D volume/area update (small strain): J ≈ 1 + tr(dε)
      // p_volume[p] *= (T{1} + p_dStrain[p].xx() + p_dStrain[p].yy());

      // Large-strain kinematics:
      //   F_{n+1} = F_n (I + L dt)
      //   V = det(F) V0
      p_deformGradient[p] =
          p_deformGradient[p] * (Matrix<T, 2, 2>::identity() + L * m_dt);
      const T J = det(p_deformGradient[p]);
      p_particles[p].volume = J * p_volume0[p];

      // Strain increment (small-strain measure) from nodal velocities:
      //   dε = sym(L) dt
      // This matches the Python reference (dEps = dt * 0.5 * (Lp + Lp.T)).
      p_strainRate[p] = Matrix<T, 2, 2>::zero();
      p_strainRate[p].xx() = L.xx();
      p_strainRate[p].yy() = L.yy();
      const T shear = T{0.5} * (L.xy() + L.yx());
      p_strainRate[p].xy() = shear;
      p_strainRate[p].yx() = shear;
      p_dStrain[p] = p_strainRate[p] * m_dt;

      if (m_law) {
        // Custom constitutive law: returns stress increment Δσ for given Δε.
        p_strain[p] += p_dStrain[p];
        p_stress[p] += m_law(p_dStrain[p]);
      } else {
        // Default: Drucker–Prager return mapping (same structure as Python).
        const Matrix<T, 3, 3> D = elasticityMatrix(m_E, m_nu, "planeStrain");
        const auto [alpha, k] = druckerPrager(m_phi, m_c);

        const DynamicVector<T> stress_n{p_stress[p].xx(), p_stress[p].yy(),
                                        p_stress[p].xy()};
        const DynamicVector<T> strain_n{p_strain[p].xx(), p_strain[p].yy(),
                                        p_strain[p].xy()};
        const DynamicVector<T> strain_increment{
            p_dStrain[p].xx(), p_dStrain[p].yy(), p_dStrain[p].xy()};

        const auto [stress_updated, strain_updated, _delta_lambda] =
            updateStressStrainDruckerPrager(stress_n, strain_n,
                                            strain_increment, D, alpha, k,
                                            "planeStrain", m_nu);

        p_stress[p].xx() = stress_updated[0];
        p_stress[p].yy() = stress_updated[1];
        p_stress[p].xy() = stress_updated[2];
        p_stress[p].yx() = stress_updated[2];

        p_strain[p].xx() = strain_updated[0];
        p_strain[p].yy() = strain_updated[1];
        p_strain[p].xy() = strain_updated[2];
        p_strain[p].yx() = strain_updated[2];
      }
    }
  }

  void resetMesh() {
    // Save updated MP state back to mesh
    m_mesh.setMPs(p_particles);

    n_mass.resetZero();
    n_momentum.resetZero();
    n_velocity.resetZero();
    n_acceleration.resetZero();
    n_displacement.resetZero();

    n_bodyForce.resetZero();
    n_tractionForce.resetZero();
    n_forceExternal.resetZero();
    n_forceInternal.resetZero();
    n_forceTotal.resetZero();

    for (auto &s : n_stress) {
      s = Matrix<T, 2, 2>::zero();
    }

    m_mesh.nodalReset();
  }

  void exportResult(const std::string &filename = "mpm2D_results.vtk") {
    const std::string vtkFile =
        endsWith(filename, ".vtk")
            ? filename
            : withExtensionReplaced(filename, ".txt", ".vtk");

    // Also export the background grid mesh (nodes + elements) to a sibling
    // file. Default naming:
    //   - if vtkFile contains "particles_" => replace it by "mesh_"
    //   - else => append "_mesh" before the extension
    std::string meshFile = vtkFile;
    if (const std::size_t pos = meshFile.rfind("particles_");
        pos != std::string::npos) {
      meshFile.replace(pos, std::string("particles_").size(), "mesh_");
    } else if (endsWith(meshFile, ".vtk")) {
      meshFile = meshFile.substr(0, meshFile.size() - 4) + "_mesh.vtk";
    } else {
      meshFile += "_mesh.vtk";
    }

    std::ofstream vtk(vtkFile);
    if (!vtk)
      throw std::runtime_error("MPM2D::exportResult: cannot open file: " +
                               vtkFile);

    const Index nPoints = p_particles.size();

    vtk << "# vtk DataFile Version 3.0\n";
    vtk << "MPM2D particles\n";
    vtk << "ASCII\n";
    vtk << "DATASET POLYDATA\n";
    vtk << "POINTS " << nPoints << " double\n";
    for (Index p = 0; p < nPoints; ++p) {
      const T x = p_particles[p].pos.x();
      const T y = p_particles[p].pos.y();
      vtk << x << " " << y << " 0\n";
    }

    vtk << "VERTICES " << nPoints << " " << (nPoints * 2) << "\n";
    for (Index p = 0; p < nPoints; ++p) {
      vtk << "1 " << p << "\n";
    }

    vtk << "POINT_DATA " << nPoints << "\n";

    vtk << "VECTORS velocity double\n";
    for (Index p = 0; p < nPoints; ++p) {
      const T vx = p_particles[p].vel.x();
      const T vy = p_particles[p].vel.y();
      vtk << vx << " " << vy << " 0\n";
    }

    vtk << "SCALARS mass double 1\n";
    vtk << "LOOKUP_TABLE default\n";
    for (Index p = 0; p < nPoints; ++p) {
      vtk << p_particles[p].mass << "\n";
    }

    vtk << "SCALARS volume double 1\n";
    vtk << "LOOKUP_TABLE default\n";
    for (Index p = 0; p < nPoints; ++p) {
      vtk << p_particles[p].volume << "\n";
    }

    vtk << "TENSORS stress double\n";
    for (Index p = 0; p < nPoints; ++p) {
      const auto &s = p_stress[p];
      vtk << s.xx() << " " << s.xy() << " 0\n";
      vtk << s.yx() << " " << s.yy() << " 0\n";
      vtk << "0 0 0\n";
    }

    vtk << "TENSORS strain double\n";
    for (Index p = 0; p < nPoints; ++p) {
      const auto &e = p_strain[p];
      vtk << e.xx() << " " << e.xy() << " 0\n";
      vtk << e.yx() << " " << e.yy() << " 0\n";
      vtk << "0 0 0\n";
    }

    // ---- Mesh grid export (VTK legacy, UNSTRUCTURED_GRID) ----
    std::ofstream mesh(meshFile);
    if (!mesh)
      throw std::runtime_error("MPM2D::exportResult: cannot open mesh file: " +
                               meshFile);

    const Index nNodes = m_mesh.getNumNodes();
    const Index nElems = m_mesh.getNumElements();

    mesh << "# vtk DataFile Version 3.0\n";
    mesh << "MPM2D background grid\n";
    mesh << "ASCII\n";
    mesh << "DATASET UNSTRUCTURED_GRID\n";
    mesh << "POINTS " << nNodes << " double\n";
    for (Index i = 0; i < nNodes; ++i) {
      const auto [x, y] = m_mesh.getNodeCoor(i);
      mesh << x << " " << y << " 0\n";
    }

    // Each quad cell line: "4 n1 n2 n3 n4" => 5 integers per element
    mesh << "CELLS " << nElems << " " << (nElems * 5) << "\n";
    for (Index e = 0; e < nElems; ++e) {
      const auto &conn = m_mesh.getEleConnectivity(e);
      mesh << "4 " << conn[0] << " " << conn[1] << " " << conn[2] << " "
           << conn[3] << "\n";
    }

    // VTK cell type for quad = 9
    mesh << "CELL_TYPES " << nElems << "\n";
    for (Index e = 0; e < nElems; ++e) {
      mesh << "9\n";
    }

    mesh << "POINT_DATA " << nNodes << "\n";

    mesh << "SCALARS active int 1\n";
    mesh << "LOOKUP_TABLE default\n";
    const auto &active = m_mesh.getActiveNodes();
    for (Index i = 0; i < nNodes; ++i) {
      const int a = (i < active.size() && active[i] != 0) ? 1 : 0;
      mesh << a << "\n";
    }

    mesh << "SCALARS mass double 1\n";
    mesh << "LOOKUP_TABLE default\n";
    for (Index i = 0; i < nNodes; ++i) {
      const T m = (i < n_mass.size()) ? n_mass[i] : T{0};
      mesh << m << "\n";
    }

    mesh << "VECTORS velocity double\n";
    for (Index i = 0; i < nNodes; ++i) {
      const T vx = (i < n_velocity.size() && n_velocity[i].size() > 0)
                       ? n_velocity[i][0]
                       : T{0};
      const T vy = (i < n_velocity.size() && n_velocity[i].size() > 1)
                       ? n_velocity[i][1]
                       : T{0};
      mesh << vx << " " << vy << " 0\n";
    }

    mesh << "VECTORS acceleration double\n";
    for (Index i = 0; i < nNodes; ++i) {
      const T ax = (i < n_acceleration.size() && n_acceleration[i].size() > 0)
                       ? n_acceleration[i][0]
                       : T{0};
      const T ay = (i < n_acceleration.size() && n_acceleration[i].size() > 1)
                       ? n_acceleration[i][1]
                       : T{0};
      mesh << ax << " " << ay << " 0\n";
    }

    mesh << "VECTORS force_total double\n";
    for (Index i = 0; i < nNodes; ++i) {
      const T fx = (i < n_forceTotal.size() && n_forceTotal[i].size() > 0)
                       ? n_forceTotal[i][0]
                       : T{0};
      const T fy = (i < n_forceTotal.size() && n_forceTotal[i].size() > 1)
                       ? n_forceTotal[i][1]
                       : T{0};
      mesh << fx << " " << fy << " 0\n";
    }

    mesh << "VECTORS momentum double\n";
    for (Index i = 0; i < nNodes; ++i) {
      const T px = (i < n_momentum.size() && n_momentum[i].size() > 0)
                       ? n_momentum[i][0]
                       : T{0};
      const T py = (i < n_momentum.size() && n_momentum[i].size() > 1)
                       ? n_momentum[i][1]
                       : T{0};
      mesh << px << " " << py << " 0\n";
    }
  }

  void exportVTKFrame(const std::filesystem::path &outputDir, Index frame,
                      int padWidth = 6) {
    std::error_code ec;
    std::filesystem::create_directories(outputDir, ec);

    std::ostringstream name;
    name << "particles_" << std::setw(padWidth) << std::setfill('0')
         << static_cast<long long>(frame) << ".vtk";
    exportResult((outputDir / name.str()).string());
  }
  void timeIntegration() {}
};

#endif // MATERIAL_POINT_METHOD_H
