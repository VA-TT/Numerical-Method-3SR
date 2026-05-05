#ifndef MATERIAL_POINT_METHOD_H
#define MATERIAL_POINT_METHOD_H

#include "ContinuumMechanic.h"
#include "Matrix.h"
#include "Mesh.h"
#include "ParentElement.h"
#include "Vector.h"
#include "elasticity.h"
#include "gaussQuadrature.h"
#include "physicConstants.h"
#include "signFunction.h"
#include <cassert>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <span>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

template <typename T, T length, Index nNodes, Index nMPperEle> class MPM1D {
private:
  // Physical properties
  T m_E{};      // Module Young
  T m_volume{}; // Volume
  T m_rho{};    // Density
  T m_mass{};   // mass
  T m_G{};      // Gravity Acceleration

  // Simulation properties
  T m_currentTime{}; // Current time
  T m_dt{};          // Time step
  T m_duration{};    // Duration of simulation
  Index m_nStep{};   // Number of steps

  //  Mesh: included nodes and material points
  Mesh1D<T> m_mesh{};
  // Inner solvers:
  std::span<Node1D<T>> m_nodes;   // Use in main loops for performance & when
  std::span<Particle1D<T>> m_MPs; // index is already controlled by the group

  // Analytical solution (if available)
  std::function<T(T)> m_analyticSolution{};
  // Behavior law
  std::function<T(T)> m_law{};

  // Helpers for saving vtk output
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
  MPM1D(T rho, T E, T dt, T duration, T v0 = {})
      : m_rho{rho}, m_E{E}, m_dt{dt}, m_duration{duration},
        m_volume{length * 1.0}, m_mass{rho * m_volume}, // 1D
        m_mesh{length, nNodes, nMPperEle},
        m_nodes{std::span(m_mesh.getAllNodes())},
        m_MPs{std::span(m_mesh.getAllMPs())} {
    // Check critical time
    T c = std::sqrt(m_E / rho);
    T dt_crit = length / c;
    assert((static_cast<double>(dt_crit) / 10.0) >= static_cast<double>(dt) &&
           "Time step isn't satisfied CFL condition (too big)");
    m_nStep = static_cast<Index>(duration / dt);

    // Distribute equally mass and volume to MPs
    Index nMPs = m_mesh.getNumMPs();
    for (Index p{0}; p < nMPs; ++p) {
      m_MPs[p].V = m_volume / nMPs;
      m_MPs[p].m = m_mass / nMPs; // MP's mass is constant
      m_MPs[p].v = v0;
    }

    m_mesh.print();
  }

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
  T getVolume() const { return m_volume; }
  T getCurrentTime() const { return m_currentTime; }
  T getTimeStep() const { return m_dt; }
  T getDuration() const { return m_duration; }
  T getNumSteps() const { return m_nStep; }

  const Mesh1D<T> &getMesh() const { return m_mesh; }
  Index getNumNodes() const { return m_mesh.getNumNodes(); }
  Index getNumElements() const { return m_mesh.getNumElements(); }
  Index getNumMPs() const { return m_mesh.getNumMPs(); }

  T getNodalMass(Index i) const { return m_mesh.getNode(i).m; }
  T getNodalVelocity(Index i) const { return m_mesh.getNode(i).v; }
  T getNodalMomentum(Index i) const { return m_mesh.getNode(i).p; }
  T getNodalExtForce(Index i) const { return m_mesh.getNode(i).extF; }
  T getNodalIntForce(Index i) const { return m_mesh.getNode(i).intF; }
  T getNodalTotalForce(Index i) const { return m_mesh.getNode(i).totF; }

  T getMPvolume(Index p) const { return m_mesh.getMP(p).V; }
  T getMPmass(Index p) const { return m_mesh.getMP(p).m; }
  T getMPvelocity(Index p) const { return m_mesh.getMP(p).v; }
  T getMPposition(Index p) const { return m_mesh.getMP(p).pos; }
  T getMPmomentum(Index p) const { return m_mesh.getMP(p).p; }
  T getMPstrain(Index p) const { return m_mesh.getMP(p).eps; }
  T getMPstrainRate(Index p) const { return m_mesh.getMP(p).epsDot; }
  T getMPdStrain(Index p) const { return m_mesh.getMP(p).dEps; }
  T getMPstress(Index p) const { return m_mesh.getMP(p).sig; }

  // Setters
  void setE(T E) { m_E = E; }
  void setG(T G) { m_G = G; } // Set G if considering gravity
  void setAnalyticSolution(std::function<T(T)> sol) {
    m_analyticSolution = sol;
  }
  void setComportmentLaw(std::function<T(T)> law) { m_law = law; }
  void setMPvelocity(Index p, T value) { m_MPs[p].v = value; }

  // Setters: store constraint values + mark constrained
  void setNodalVeloConstraint(Index i, T value) {
    m_mesh.getNode(i).vConVal = value;
    m_mesh.getNode(i).vCon = 1;
  }
  void setNodalAccConstraint(Index i, T value) {
    m_mesh.getNode(i).aConVal = value;
    m_mesh.getNode(i).aCon = 1;
  }
  void setNodalMomentumConstraint(Index i, T value) {
    m_mesh.getNode(i).pConVal = value;
    m_mesh.getNode(i).pCon = 1;
  }
  void setNodalForceConstraint(Index i, T value) {
    m_mesh.getNode(i).fConVal = value;
    m_mesh.getNode(i).fCon = 1;
  }

  // Apply stored constraints to current nodal state
  void applyNodalVeloConstraint() {
    for (Index i{0}; i < getNumNodes(); ++i) {
      if (m_nodes[i].vCon) {
        m_nodes[i].v = m_nodes[i].vConVal;
      }
    }
  }
  void applyNodalAccConstraint() {
    for (Index i{0}; i < getNumNodes(); ++i) {
      if (m_nodes[i].aCon) {
        m_nodes[i].a = m_nodes[i].aConVal;
      }
    }
  }
  void applyNodalMomentumConstraint() {
    for (Index i{0}; i < getNumNodes(); ++i) {
      if (m_nodes[i].pCon) {
        m_nodes[i].P = m_nodes[i].pConVal;
      }
    }
  }
  void applyNodalForceConstraint() {
    for (Index i{0}; i < getNumNodes(); ++i) {
      if (m_nodes[i].fCon) {
        m_nodes[i].totF = m_nodes[i].fConVal;
      }
    }
  }

  void setupMP() { m_mesh.activateNodesAndElements(); }

  void p2n() {
    Index nMPs = m_mesh.getNumMPs();
    // Map nodal mass + momentum
    for (Index p{0}; p < nMPs; ++p) {
      m_MPs[p].P = m_MPs[p].m * m_MPs[p].v;
      Index e = m_mesh.getMP(p).eleID;
      if (e != idError) {
        T x_p = m_MPs[p].pos;
        auto conn = m_mesh.getEleConnectivity(e);
        auto ele = m_mesh.getElement(e);
        Index n1 = ele.n1;
        Index n2 = ele.n2;
        T xi = ele.parentCoord(x_p);
        m_nodes[n1].m += ele.N1_ref(xi) * m_MPs[p].m;
        m_nodes[n2].m += ele.N2_ref(xi) * m_MPs[p].m;
        m_nodes[n1].P += ele.N1_ref(xi) * m_MPs[p].P;
        m_nodes[n2].P += ele.N2_ref(xi) * m_MPs[p].P;
      }
    }
    applyNodalMomentumConstraint();
  }

  void nodalEquilibrium() {
    // f^ext = b + t
    for (Index p{0}; p < getNumMPs(); ++p) {
      Index e = m_mesh.getMP(p).eleID;
      if (e != idError) {
        T x_p = m_MPs[p].pos;
        auto conn = m_mesh.getEleConnectivity(e);
        auto ele = m_mesh.getElement(e);
        Index n1 = ele.n1;
        Index n2 = ele.n2;
        T xi = ele.parentCoord(x_p);
        m_nodes[n1].bodyF += m_G * ele.N1_ref(xi) * m_MPs[p].m;
        m_nodes[n2].bodyF += m_G * ele.N2_ref(xi) * m_MPs[p].m;
        // Traction force t_i (to be implemented)
        m_nodes[n1].intF -= m_MPs[p].V * ele.dN1_dx() * m_MPs[p].sig;
        m_nodes[n2].intF -= m_MPs[p].V * ele.dN2_dx() * m_MPs[p].sig;
      }
    }
    for (Index i{0}; i < getNumNodes(); ++i) {
      m_nodes[i].extF = m_nodes[i].bodyF + m_nodes[i].tracF;
      m_nodes[i].totF = m_nodes[i].extF + m_nodes[i].intF;
    }
    // Enforce any stored nodal force constraints after assembly
    applyNodalForceConstraint();
  }

  void n2p() {
    // Update momentum at nodes
    for (Index i{0}; i < getNumNodes(); ++i) {
      if (m_nodes[i].isActive) {
        m_nodes[i].P += m_nodes[i].totF * m_dt;
        m_nodes[i].a = m_nodes[i].totF / m_nodes[i].m;

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
    for (Index p{0}; p < getNumMPs(); ++p) {
      Index e = m_mesh.getMP(p).eleID;
      if (e != idError) {
        T x_p = m_MPs[p].pos;
        ElementL2<T> ele = m_mesh.getElement(e);
        Index n1 = ele.n1;
        Index n2 = ele.n2;
        T xi = ele.parentCoord(x_p);
        // Update velocity and position using FLIP style
        m_MPs[p].a =
            ele.N1_ref(xi) * m_nodes[n1].a + ele.N2_ref(xi) * m_nodes[n2].a;
        m_MPs[p].v += m_MPs[p].a * m_dt;
        // Hybrid
        // T v_pic = N1 * velocity_n[n1] + N2 * velocity_n[n2];
        // T v_flip = velocity_p[p] + (N1 * a_n1 + N2 * a_n2) * dt;
        // velocity_p[p] = alpha * v_pic + '(1-alpha)' * v_flip;
        // position_p[p] += velocity_p[p] * m_dt;
        m_MPs[p].pos += (ele.N1_ref(xi) * m_nodes[n1].P / m_nodes[n1].m +
                         ele.N2_ref(xi) * m_nodes[n2].P / m_nodes[n2].m) *
                        m_dt;
        m_MPs[p].P = m_MPs[p].m * m_MPs[p].v;
      }
    }

    for (Index p{0}; p < getNumMPs(); ++p) {
      Index e = m_mesh.getMP(p).eleID;
      if (e != idError) {
        T x_p = m_MPs[p].pos;
        ElementL2<T> ele = m_mesh.getElement(e);
        Index n1 = ele.n1;
        Index n2 = ele.n2;
        T xi = ele.parentCoord(x_p);
        m_nodes[n1].v += m_MPs[p].P * ele.N1_ref(xi) / m_nodes[n1].m;
        m_nodes[n2].v += m_MPs[p].P * ele.N2_ref(xi) / m_nodes[n2].m;
      }
    }
    applyNodalVeloConstraint();

    for (Index p{0}; p < m_mesh.getNumMPs(); ++p) {
      Index e = m_mesh.getMP(p).eleID;
      if (e != idError) {
        T x_p = m_MPs[p].pos;
        ElementL2<T> ele = m_mesh.getElement(e);
        Index n1 = ele.n1;
        Index n2 = ele.n2;
        T xi = ele.parentCoord(x_p);
        // Attention: x1-x2 belong to nodes (their positions don't get
        // updated), while the velocity is measured at MPs (updated at t+dt)
        m_MPs[p].epsDot =
            ele.dN1_dx() * m_nodes[n1].v + ele.dN2_dx() * m_nodes[n2].v;
        m_MPs[p].dEps = m_MPs[p].epsDot * m_dt;
        m_MPs[p].eps += m_MPs[p].dEps;
      }
      // Constitutive law:
      if (m_law) {
        m_MPs[p].sig += m_law(m_MPs[p].dEps);
      } else {
        m_MPs[p].sig += m_E * m_MPs[p].dEps; // Default linear elastic
      }
      // Update volume
      m_MPs[p].V *= (T{1} + m_MPs[p].dEps);
    }
  }

  void resetMesh() {
    for (Index i{0}; i < nNodes; ++i) {
      auto &node = m_nodes[i];
      node.pos = node.posInit;
      node.m = node.v = node.a = node.P = T{};
      node.bodyF = node.tracF = node.extF = node.intF = node.totF = T{};
    }
  } // Excluding MPs

  // void compareAnalytic(T xloc) {
  //   if (!m_analyticSolution) {
  //     std::cout << "\nNo analytical solution provided.\n";
  //     return;
  //   }

  //   std::cout << "\n=== Comparison with Analytical Solution ===\n";
  //   std::cout << std::setw(10) << "Time" << std::setw(15) << "x_numerical"
  //             << std::setw(15) << "x_exact" << std::setw(15) << "Error" <<
  //             '\n';
  //   std::cout << std::string(55, '-') << '\n';

  //   // Find MP closest to m_xloc
  //   T min_dist = std::abs(m_MPs[0].pos - xloc);
  //   Index closest_p = 0;
  //   for (Index p = 1; p < m_mesh.getNumMPs(); ++p) {
  //     T dist = std::abs(m_MPs[p].pos - xloc);
  //     if (dist < min_dist) {
  //       min_dist = dist;
  //       closest_p = p;
  //     }
  //   }
  //   T x_numerical = m_MPs[closest_p].pos;

  //   T x_exact = m_analyticSolution(m_currentTime);
  //   T error = std::abs(x_numerical - x_exact);

  //   std::cout << std::setw(10) << std::fixed << std::setprecision(4)
  //             << m_currentTime << std::setw(15) << std::setprecision(6)
  //             << x_numerical << std::setw(15) << x_exact << std::setw(15)
  //             << std::scientific << error << '\n';
  // }

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

    const Index numMPs = getNumMPs();
    const Index numNodes = getNumNodes();
    const Index numEles = getNumElements();

    vtk << "# vtk DataFile Version 3.0\n";
    vtk << "MPM1D particles\n";
    vtk << "ASCII\n";
    vtk << "DATASET POLYDATA\n";
    vtk << "POINTS " << numMPs << " double\n";
    for (Index p = 0; p < numMPs; ++p) {
      vtk << m_MPs[p].pos << " 0 0\n";
    }

    vtk << "VERTICES " << numMPs << " " << (numMPs * 2) << "\n";
    for (Index p = 0; p < numMPs; ++p) {
      vtk << "1 " << p << "\n";
    }

    vtk << "POINT_DATA " << numMPs << "\n";
    vtk << "VECTORS velocity double\n";
    for (Index p = 0; p < numMPs; ++p) {
      vtk << m_MPs[p].v << " 0 0\n";
    }

    vtk << "SCALARS mass double 1\n";
    vtk << "LOOKUP_TABLE default\n";
    for (Index p = 0; p < numMPs; ++p) {
      vtk << m_MPs[p].m << "\n";
    }

    vtk << "SCALARS volume double 1\n";
    vtk << "LOOKUP_TABLE default\n";
    for (Index p = 0; p < numMPs; ++p) {
      vtk << m_MPs[p].V << "\n";
    }

    // ---- Mesh grid export (VTK legacy, UNSTRUCTURED_GRID) ----
    std::ofstream mesh(meshFile);
    if (!mesh)
      throw std::runtime_error("MPM1D::exportResult: cannot open mesh file: " +
                               meshFile);

    mesh << "# vtk DataFile Version 3.0\n";
    mesh << "MPM1D background grid\n";
    mesh << "ASCII\n";
    mesh << "DATASET UNSTRUCTURED_GRID\n";
    mesh << "POINTS " << numNodes << " double\n";
    for (Index i = 0; i < numNodes; ++i) {
      mesh << m_nodes[i].pos << " 0 0\n";
    }

    // Each line cell: "2 n1 n2" => 3 integers per element
    mesh << "CELLS " << numEles << " " << (numEles * 3) << "\n";
    for (Index e = 0; e < numEles; ++e) {
      const auto &conn = m_mesh.getEleConnectivity(e);
      mesh << "2 " << conn[0] << " " << conn[1] << "\n";
    }

    // VTK cell type for line = 3
    mesh << "CELL_TYPES " << numEles << "\n";
    for (Index e = 0; e < numEles; ++e) {
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
      const T mass = m_nodes[i].m;
      mesh << mass << "\n";
    }

    mesh << "VECTORS velocity double\n";
    for (Index i = 0; i < numNodes; ++i) {
      const T v = m_nodes[i].v;
      mesh << v << " 0 0\n";
    }

    mesh << "VECTORS acceleration double\n";
    for (Index i = 0; i < numNodes; ++i) {
      const T a = m_nodes[i].a;
      mesh << a << " 0 0\n";
    }

    mesh << "VECTORS force_total double\n";
    for (Index i = 0; i < numNodes; ++i) {
      const T f = m_nodes[i].totF;
      mesh << f << " 0 0\n";
    }

    mesh << "VECTORS momentum double\n";
    for (Index i = 0; i < numNodes; ++i) {
      const T p = m_nodes[i].P;
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

  void timeIntegration() {}
};

//================================================================================//

//                                    2D

//================================================================================//

template <typename T, T gridLength, T gridHeight, Index nx, Index ny, T MP_size>
class MPM2D {
private:
  // Particles, nodes and elements
  Mesh2D<T> m_mesh{};
  std::span<Node2D<T>> m_nodes;
  std::span<Particle2D<T>> m_MPs;

  // Material properties
  T m_rho{}; // Density
  T m_E{};   // Module Young
  T m_nu{};  // Poisson ratio
  T m_phi{}; // Internal Friction angle
  T m_mu{};  // Material-Boundary Friction
  T m_c{};   // cohesion
  T m_K0{};  // Initial earth pressure coefficient

  // MP domain
  T m_pLength{}, m_pHeight{};
  T m_volume{}, m_volume0{}, m_mass{};

  // Kinetic, Potential, Initial Total and Dissipation energies
  T m_kinEnergy{}, m_potentEnergy{}, m_totalEnergy0{}, m_dissiEnergy{};

  T m_G{}; // Gravity Acceleration

  // Time advancement
  Index m_nSteps{};
  T m_dt{}, m_duration{}, m_currentTime{};

  // Constitutive law
  std::function<Matrix<T, 2, 2>(const Matrix<T, 2, 2> &)> m_law{};

  // Boundary contact
  static constexpr Index dimensions = 2;
  static constexpr Index dir_x = 0;
  static constexpr Index dir_y = 1;
  static constexpr Index directionPositive = 1;
  static constexpr Index directionNegative = -1;

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

  void applyNodalVeloConstraint() {
    for (Index i{0}; i < getNumNodes(); ++i) {
      for (Index d{0}; d < dimensions; ++d) {
        if (m_nodes[i].vCon[d]) {
          m_nodes[i].v[d] = m_nodes[i].vConVal[d];
          // m_nodes[i].P[d] = m_nodes[i].v[d] * m_nodes[i].mass;
        }
      }
    }
  }

  void applyNodalAccConstraint() {
    for (Index i{0}; i < getNumNodes(); ++i) {
      for (Index d{0}; d < dimensions; ++d) {
        if (m_nodes[i].aCon[d]) {
          m_nodes[i].a[d] = m_nodes[i].aConVal[d];
          // m_nodes[i].totF[d] = m_nodes[i].a[d] * m_nodes[i].mass;
        }
      }
    }
  }

  void applyNodalMomentumConstraint() {
    for (Index i{0}; i < getNumNodes(); ++i) {
      for (Index d{0}; d < dimensions; ++d) {
        if (m_nodes[i].pCon[d]) {
          m_nodes[i].P[d] = m_nodes[i].pConVal[d];
        }
      }
    }
  }

  void applyNodalForceConstraint() {
    for (Index i{0}; i < getNumNodes(); ++i) {
      for (Index d{0}; d < dimensions; ++d) {
        if (m_nodes[i].fCon[d]) {
          m_nodes[i].totF[d] = m_nodes[i].fConVal[d];
        }
      }
    }
  }

  // Source: https://www.geoelements.org/LearnMPM/mpm2d-column-collapse.html
  void frictionalBC(const DynamicVector<Index> &nodeIDs, Index dir_n,
                    Index signDir_n) {
    const Index dir_t = 1 - dir_n; // tangential direction

    for (const Index nodeID : nodeIDs) {
      auto &node = m_nodes[nodeID];
      if (approximatelyEqualAbsRel(node.mass, T{})) {
        continue;
      }

      T acc_n = node.a[dir_n];
      T acc_t = node.a[dir_t];
      const T vel_t = node.P[dir_t] / node.mass;
      const bool moveTowardBoundary = (acc_n * static_cast<T>(signDir_n)) > T{};

      if (moveTowardBoundary) {
        if (!approximatelyEqualAbsRel(vel_t, T{})) {
          const T vel_net = m_dt * acc_t + vel_t;
          const T vel_frictional = m_dt * m_mu * std::abs(acc_n);
          if (std::abs(vel_net) <= vel_frictional) {
            acc_t = -vel_t / m_dt;
          } else {
            acc_t -= sgn(vel_net) * m_mu * std::abs(acc_n);
          }
        } else {
          if (std::abs(acc_t) <= m_mu * std::abs(acc_n)) {
            acc_t = T{};
          } else {
            acc_t -= sgn(acc_t) * m_mu * std::abs(acc_n);
          }
        }
      }

      node.a[dir_t] = acc_t;
      node.totF = node.a * node.mass;
    }
  }

  void applyFrictionalBC() {
    frictionalBC(m_mesh.leftActiveNodes(), dir_x, directionNegative);
    frictionalBC(m_mesh.rightActiveNodes(), dir_x, directionPositive);
    frictionalBC(m_mesh.bottomActiveNodes(), dir_y, directionNegative);
    frictionalBC(m_mesh.topActiveNodes(), dir_y, directionPositive);

    for (Index i : m_mesh.leftActiveNodes()) {
      m_nodes[i].P.x() = T{};
      m_nodes[i].totF.x() = T{};
    }
    for (Index i : m_mesh.rightActiveNodes()) {
      m_nodes[i].P.x() = T{};
      m_nodes[i].totF.x() = T{};
    }
    for (Index i : m_mesh.bottomActiveNodes()) {
      m_nodes[i].P.y() = T{};
      m_nodes[i].totF.y() = T{};
    }
    for (Index i : m_mesh.topActiveNodes()) {
      m_nodes[i].P.y() = T{};
      m_nodes[i].totF.y() = T{};
    }
  }

  void initializeStress() {
    const T ymax = m_pHeight - MP_size / T{2};
    for (Index p{0}; p < getNumMPs(); ++p) {
      m_MPs[p].sig.yy() = -m_G * m_rho * (ymax - m_MPs[p].pos.y());
      m_MPs[p].sig.xx() = m_K0 * m_MPs[p].sig.yy();
    }
  }

  void computeEnergy() {
    m_potentEnergy = T{};
    m_kinEnergy = T{};
    for (Index p{0}; p < getNumMPs(); ++p) {
      m_potentEnergy += m_G * m_MPs[p].m * m_MPs[p].pos.y();
      m_kinEnergy += T{0.5} * m_MPs[p].m * dotProduct(m_MPs[p].v, m_MPs[p].v);
    }
    m_dissiEnergy = m_totalEnergy0 - m_potentEnergy - m_kinEnergy;
  }

public:
  MPM2D(T rho, T E, T nu, T phi, T c, T mu, const StaticVector<T, 2> &minCorner,
        const StaticVector<T, 2> &maxCorner, T dt, T duration, T v0 = T{})
      : m_E{E}, m_nu{nu}, m_rho{rho}, m_mu{mu}, m_phi{phi}, m_c{c},
        m_K0{nu / (T{1} - nu)}, m_pLength{maxCorner.x() - minCorner.x()},
        m_pHeight{maxCorner.y() - minCorner.y()},
        m_volume{m_pLength * m_pHeight}, m_mass{rho * m_volume}, m_dt{dt},
        m_duration{duration}, m_nSteps{static_cast<Index>(duration / dt)},
        m_mesh{gridLength,
               gridHeight,
               nx,
               ny,
               minCorner.x(),
               minCorner.y(),
               maxCorner.x(),
               maxCorner.y(),
               MP_size},
        m_nodes{std::span(m_mesh.getAllNodes())},
        m_MPs{std::span(m_mesh.getAllMPs())} {
    const T cWave = std::sqrt(m_E / m_rho);
    // const T dx = gridLength / static_cast<T>(nx - 1);
    // const T dy = gridHeight / static_cast<T>(ny - 1);
    const T dt_crit = (m_pLength < m_pHeight ? m_pLength : m_pHeight) / cWave;
    assert((dt_crit / T{10}) >= dt &&
           "dt doesn't satisfied CFL condition (too big)");

    const Index nMPs = m_mesh.getNumMPs();
    for (Index p{0}; p < nMPs; ++p) {
      m_MPs[p].V = m_volume / static_cast<T>(nMPs);
      m_MPs[p].V0 = m_MPs[p].V;
      m_MPs[p].m = m_mass / static_cast<T>(nMPs);
      m_MPs[p].v.x() = v0;
      m_MPs[p].v.y() = T{};
      m_MPs[p].F = Matrix<T, 2, 2>::identity();
    }

    for (Index p{0}; p < nMPs; ++p) {
      m_totalEnergy0 += m_G * m_MPs[p].m * m_MPs[p].pos.y();
    }
  }

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
  T getTimeStep() const { return m_dt; }
  T getDuration() const { return m_duration; }
  Index getNumSteps() const { return m_nSteps; }

  const Mesh2D<T> &getMesh() const { return m_mesh; }
  Mesh2D<T> &getMesh() { return m_mesh; }
  Index getNumNodes() const { return m_mesh.getNumNodes(); }
  Index getNumElements() const { return m_mesh.getNumElements(); }
  Index getNumMPs() const { return m_mesh.getNumMPs(); }

  T getMPvolume(Index p) const { return m_mesh.getMP(p).V; }
  T getMPmass(Index p) const { return m_mesh.getMP(p).m; }
  StaticVector<T, 2> getMPvelocity(Index p) const { return m_mesh.getMP(p).v; }
  StaticVector<T, 2> getMPposition(Index p) const {
    return m_mesh.getMP(p).pos;
  }
  T getMPstrain(Index p) const { return m_mesh.getMP(p).eps.xx(); }
  T getMPstress(Index p) const { return m_mesh.getMP(p).sig.xx(); }

  // Setters
  void setE(T E) { m_E = E; }
  void setG(T G) { m_G = G; }

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
      Matrix<T, 2, 2> dsig{};
      dsig.xx() = law(dEps.xx());
      dsig.yy() = law(dEps.yy());
      dsig.xy() = law(dEps.xy());
      dsig.yx() = law(dEps.yx());
      return dsig;
    };
  }

  void setNodalVeloConstraint(Index i, const StaticVector<T, 2> &value) {
    m_nodes[i].vCon = StaticVector<char, 2>{1, 1};
    m_nodes[i].vConVal = value;
  }
  void setNodalVeloConstraint(Index i, T value) {
    setNodalVeloConstraint(i, StaticVector<T, 2>{value, value});
  }
  void setNodalAccConstraint(Index i, const StaticVector<T, 2> &value) {
    m_nodes[i].aCon = StaticVector<char, 2>{1, 1};
    m_nodes[i].aConVal = value;
  }
  void setNodalAccConstraint(Index i, T value) {
    setNodalAccConstraint(i, StaticVector<T, 2>{value, value});
  }
  void setNodalMomentumConstraint(Index i, const StaticVector<T, 2> &value) {
    m_nodes[i].pCon = StaticVector<char, 2>{1, 1};
    m_nodes[i].pConVal = value;
  }
  void setNodalMomentumConstraint(Index i, T value) {
    setNodalMomentumConstraint(i, StaticVector<T, 2>{value, value});
  }
  void setNodalForceConstraint(Index i, const StaticVector<T, 2> &value) {
    m_nodes[i].fCon = StaticVector<char, 2>{1, 1};
    m_nodes[i].fConVal = value;
  }
  void setNodalForceConstraint(Index i, T value) {
    setNodalForceConstraint(i, StaticVector<T, 2>{value, value});
  }

  void setupMP() {
    m_mesh.activateNodesAndElements();
    m_mesh.updateAllMasks();
    initializeStress();
    computeEnergy();
  }

  void p2n() {
    for (Index p{0}; p < getNumMPs(); ++p) {
      auto &mp = m_MPs[p];
      mp.P = mp.m * mp.v;
      const Index e = mp.eleID;
      if (e == idError)
        continue;

      const auto &ele = m_mesh.getElement(e);
      const auto conn = ele.getConnectivity();
      const auto parent = ele.parentCoord(mp.pos.x(), mp.pos.y());
      const auto N = ele.N_Q4(parent.x(), parent.y());
      for (Index q{0}; q < 4; ++q) {
        const Index nodeID = conn[q];
        m_nodes[nodeID].mass += N[q] * mp.m;
        m_nodes[nodeID].P += N[q] * mp.P;
      }
    }
    applyNodalMomentumConstraint();
  }
  // //Element centric
  //   const auto &eles = m_mesh.getAllElements();
  //   for (Index e{0}; e < getNumElements(); ++e) {
  //     const auto &ele = eles[e];
  //     if (ele.isActive) {
  //       const auto conn = ele.getConnectivity();
  //       for (auto &mp : m_MPs) {
  //         if (mp.eleID != e)
  //           continue;
  //         const auto parent = ele.parentCoord(mp.pos.x(), mp.pos.y());
  //         const auto grad = ele.gradientN(parent.x(), parent.y());
  //         const auto N = ele.N_Q4(parent.x(), parent.y());
  //         const auto &dNdx = grad[0];
  //         const auto &dNdy = grad[1];
  //         for (Index a{0}; a < 4; ++a) {
  //           const Index nodeID = conn[a];
  //           m_nodes[nodeID].bodyF.y() += N[a] * m_G * mp.m;
  //           // Traction force t_i (to be implemented)
  //           const StaticVector<T, 2> gradNa{dNdx[a], dNdy[a]};
  //           m_nodes[nodeID].intF -= mp.V * (mp.sig * gradNa);
  //         }
  //       }
  //     }
  //   }
  void nodalEquilibrium() {
    for (Index p{0}; p < getNumMPs(); ++p) {
      auto &mp = m_MPs[p];
      const Index e = mp.eleID;
      if (e == idError)
        continue;

      const auto &ele = m_mesh.getElement(e);
      const auto conn = ele.getConnectivity();
      const auto parent = ele.parentCoord(mp.pos.x(), mp.pos.y());
      const auto N = ele.N_Q4(parent.x(), parent.y());
      const auto grad = ele.gradientN(parent.x(), parent.y());
      const auto &dNdx = grad[0];
      const auto &dNdy = grad[1];

      for (Index q{0}; q < 4; ++q) {
        const Index nodeID = conn[q];
        m_nodes[nodeID].bodyF.y() += N[q] * m_G * mp.m;
        // Traction force t_i (to be implemented)
        const StaticVector<T, 2> gradNq{dNdx[q], dNdy[q]};
        m_nodes[nodeID].intF -= mp.V * (mp.sig * gradNq);
      }
    }

    // assembly forces
    for (Index i{0}; i < getNumNodes(); ++i) {
      if (!m_nodes[i].isActive)
        continue;
      m_nodes[i].extF = m_nodes[i].bodyF + m_nodes[i].tracF;
      m_nodes[i].totF = m_nodes[i].extF + m_nodes[i].intF;
    }
    applyNodalForceConstraint();

    // Compute nodal acceleration
    for (Index i{0}; i < getNumNodes(); ++i) {
      if (m_nodes[i].isActive &&
          !approximatelyEqualAbsRel(m_nodes[i].mass, T{})) {
        m_nodes[i].a = m_nodes[i].totF / m_nodes[i].mass;
      }
    }
    applyNodalAccConstraint();
    applyFrictionalBC();

    // Update momentum at nodes (after frictional BC)
    for (Index i{0}; i < getNumNodes(); ++i) {
      if (m_nodes[i].isActive)
        m_nodes[i].P += m_nodes[i].totF * m_dt;
    }
    applyNodalMomentumConstraint();

    // Update nodal velocity
    for (Index i{0}; i < getNumNodes(); ++i) {
      if ((m_nodes[i].isActive) &&
          (!approximatelyEqualAbsRel(m_nodes[i].mass, T{}))) {
        m_nodes[i].v = m_nodes[i].P / m_nodes[i].mass;
      }
    }

    applyNodalVeloConstraint();
  }

  void n2p() {
    // Map back to MPs
    for (Index p{0}; p < getNumMPs(); ++p) {
      auto &mp = m_MPs[p];
      const Index e = mp.eleID;
      if (e == idError)
        continue;

      const auto &ele = m_mesh.getElement(e);
      const auto conn = ele.getConnectivity();
      const auto parent = ele.parentCoord(mp.pos.x(), mp.pos.y());
      const auto N = ele.N_Q4(parent.x(), parent.y());

      StaticVector<T, 2> aNext{};
      StaticVector<T, 2> vNext{};
      for (Index q{0}; q < 4; ++q) {
        const Index nodeID = conn[q];
        aNext += N[q] * m_nodes[nodeID].a;
        vNext += N[q] * m_nodes[nodeID].v;
      }

      // Foward Euler
      mp.a = aNext;
      mp.v += aNext * m_dt;
      mp.pos += vNext * m_dt;
      mp.P = mp.m * mp.v;

      // Update stress and strain
      const auto grad = ele.gradientN(parent.x(), parent.y());
      const auto &dNdx = grad[0];
      const auto &dNdy = grad[1];
      Matrix<T, 2, 2> L = Matrix<T, 2, 2>::zero();
      for (Index q{0}; q < 4; ++q) {
        const Index nodeID = conn[q];
        const StaticVector<T, 2> gradNq{dNdx[q], dNdy[q]};
        L += tensorProduct<T, 2, 2>(m_nodes[nodeID].v, gradNq);
      }

      mp.F = mp.F * (Matrix<T, 2, 2>::identity() + L * m_dt);
      mp.V = det(mp.F) * mp.V0; // J = det(F)

      mp.epsDot = Matrix<T, 2, 2>::zero();
      mp.epsDot.xx() = L.xx();
      mp.epsDot.yy() = L.yy();
      const T shear = T{0.5} * (L.xy() + L.yx());
      mp.epsDot.xy() = shear;
      mp.epsDot.yx() = shear;
      mp.dEps = mp.epsDot * m_dt;

      if (m_law) {
        mp.eps += mp.dEps;
        mp.sig += m_law(mp.dEps);
      } else {
        const Matrix<T, 3, 3> D = elasticityMatrix(m_E, m_nu, "planeStrain");
        const auto [alpha, k] = druckerPrager(m_phi, m_c);
        const DynamicVector<T> stressN{mp.sig.xx(), mp.sig.yy(), mp.sig.xy()};
        const DynamicVector<T> strainN{mp.eps.xx(), mp.eps.yy(), mp.eps.xy()};
        const DynamicVector<T> strainInc{mp.dEps.xx(), mp.dEps.yy(),
                                         mp.dEps.xy()};
        const auto [stressUpdated, strainUpdated, deltaLambda] =
            updateStressStrainDruckerPrager(stressN, strainN, strainInc, D,
                                            alpha, k, "planeStrain", m_nu);
        (void)deltaLambda;
        mp.sig.xx() = stressUpdated[0];
        mp.sig.yy() = stressUpdated[1];
        mp.sig.xy() = stressUpdated[2];
        mp.sig.yx() = stressUpdated[2];
        mp.eps.xx() = strainUpdated[0];
        mp.eps.yy() = strainUpdated[1];
        mp.eps.xy() = strainUpdated[2];
        mp.eps.yx() = strainUpdated[2];
      }

      if (!std::isfinite(static_cast<double>(mp.sig.xx())) ||
          !std::isfinite(static_cast<double>(mp.sig.yy())) ||
          !std::isfinite(static_cast<double>(mp.sig.xy())) ||
          !std::isfinite(static_cast<double>(mp.V))) {
        mp.sig = Matrix<T, 2, 2>::zero();
        mp.eps = Matrix<T, 2, 2>::zero();
        mp.epsDot = Matrix<T, 2, 2>::zero();
        mp.dEps = Matrix<T, 2, 2>::zero();
        mp.F = Matrix<T, 2, 2>::identity();
        mp.V = mp.V0;
      }
    }
  }

  void resetMesh() {
    for (Index i{0}; i < getNumNodes(); ++i) {
      auto &node = m_nodes[i];
      node.pos = node.posInit;
      node.mass = T{};
      node.v.resetZero();
      node.a.resetZero();
      node.P.resetZero();
      node.bodyF.resetZero();
      node.tracF.resetZero();
      node.extF.resetZero();
      node.intF.resetZero();
      node.totF.resetZero();
    }
  }

  void exportResult(const std::string &filename = "mpm2D_results.vtk") {
    const std::string vtkFile =
        endsWith(filename, ".vtk")
            ? filename
            : withExtensionReplaced(filename, ".txt", ".vtk");

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

    const Index nMPs = getNumMPs();
    vtk << "# vtk DataFile Version 3.0\n";
    vtk << "MPM2D particles\n";
    vtk << "ASCII\n";
    vtk << "DATASET POLYDATA\n";
    vtk << "POINTS " << nMPs << " double\n";
    for (Index p{0}; p < nMPs; ++p)
      vtk << m_MPs[p].pos.x() << " " << m_MPs[p].pos.y() << " 0\n";

    vtk << "VERTICES " << nMPs << " " << (nMPs * 2) << "\n";
    for (Index p{0}; p < nMPs; ++p)
      vtk << "1 " << p << "\n";

    vtk << "POINT_DATA " << nMPs << "\n";
    vtk << "VECTORS velocity double\n";
    for (Index p{0}; p < nMPs; ++p)
      vtk << m_MPs[p].v.x() << " " << m_MPs[p].v.y() << " 0\n";

    vtk << "SCALARS mass double 1\nLOOKUP_TABLE default\n";
    for (Index p{0}; p < nMPs; ++p)
      vtk << m_MPs[p].m << "\n";

    vtk << "SCALARS volume double 1\nLOOKUP_TABLE default\n";
    for (Index p{0}; p < nMPs; ++p)
      vtk << m_MPs[p].V << "\n";

    vtk << "TENSORS stress double\n";
    for (Index p{0}; p < nMPs; ++p) {
      const auto &s = m_MPs[p].sig;
      vtk << s.xx() << " " << s.xy() << " 0\n";
      vtk << s.yx() << " " << s.yy() << " 0\n";
      vtk << "0 0 0\n";
    }

    std::ofstream mesh(meshFile);
    if (!mesh)
      throw std::runtime_error("MPM2D::exportResult: cannot open mesh file: " +
                               meshFile);

    const Index nNodes = getNumNodes();
    const Index nElems = getNumElements();
    mesh << "# vtk DataFile Version 3.0\n";
    mesh << "MPM2D background grid\n";
    mesh << "ASCII\n";
    mesh << "DATASET UNSTRUCTURED_GRID\n";
    mesh << "POINTS " << nNodes << " double\n";
    for (Index i{0}; i < nNodes; ++i)
      mesh << m_nodes[i].pos.x() << " " << m_nodes[i].pos.y() << " 0\n";

    mesh << "CELLS " << nElems << " " << (nElems * 5) << "\n";
    for (Index e{0}; e < nElems; ++e) {
      const auto conn = m_mesh.getEleConnectivity(e);
      mesh << "4 " << conn[0] << " " << conn[1] << " " << conn[2] << " "
           << conn[3] << "\n";
    }

    mesh << "CELL_TYPES " << nElems << "\n";
    for (Index e{0}; e < nElems; ++e)
      mesh << "9\n";

    mesh << "POINT_DATA " << nNodes << "\n";
    mesh << "SCALARS active int 1\nLOOKUP_TABLE default\n";
    for (Index i{0}; i < nNodes; ++i)
      mesh << (m_nodes[i].isActive ? 1 : 0) << "\n";

    mesh << "SCALARS mass double 1\nLOOKUP_TABLE default\n";
    for (Index i{0}; i < nNodes; ++i)
      mesh << m_nodes[i].mass << "\n";

    mesh << "VECTORS velocity double\n";
    for (Index i{0}; i < nNodes; ++i)
      mesh << m_nodes[i].v.x() << " " << m_nodes[i].v.y() << " 0\n";
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
