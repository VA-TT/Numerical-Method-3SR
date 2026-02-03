#ifndef FINITE_ELEMENT_METHOD_H
#define FINITE_ELEMENT_METHOD_H

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

template <typename T, Index nNodes> class MPM1D {
private:
  // Physical properties
  T m_EA{10.0};      // Axial stiffness
  T m_length{1.0};   // Domain length
  T m_V{1.0};        // Volume
  T m_rho{1000};     // Density
  T m_stress{0.0};   // Strain
  T m_strain{0.0};   // Stress
  T m_velocity{0.0}; // Velocity
  T m_position{0.0}; // Velocity

  // Mesh
  Mesh1D<T> m_mesh{};

  // FEM matrices
  Matrix<T, nNodes, nNodes> m_K{}, m_K_original{};
  Matrix<T, nNodes, 1> m_U{};
  Matrix<T, nNodes, 1> m_F{}, m_F_original{};
  Matrix<T, nNodes, 1> m_R{};

  // Distributed force: f(x)
  std::function<T(T)> m_rhsFunction;

  // Analytical solution (if available)
  std::function<T(T)> m_analyticSolution;

public:
  // Constructor
  FEM1D(T EA, T length)
      : m_EA{EA}, m_length{length}, m_mesh{Mesh1D<T>{length, nNodes}} {
    // Default distributed load: f(x) = 0
    m_rhsFunction = [](T x) { return T(0); };
    m_mesh.print();
  };

  // Other defaults
  FEM1D() = default;
  FEM1D(const FEM1D &) = default;
  FEM1D(FEM1D &&) = default;
  FEM1D &operator=(const FEM1D &) = default;
  FEM1D &operator=(FEM1D &&) = default;
  ~FEM1D() = default;

  // Getters
  T getEA() const { return m_EA; }
  T getLength() const { return m_length; }
  Index getNumNodes() const { return m_mesh.getNumNodes(); }
  Index getNumElements() const { return m_mesh.getNumElements(); }
  const Mesh1D<T> &getMesh() const { return m_mesh; }
  const Matrix<T, nNodes, nNodes> &getK() const { return m_K; }
  const Matrix<T, nNodes, 1> &getU() const { return m_U; }
  const Matrix<T, nNodes, 1> &getF() const { return m_F; }
  const Matrix<T, nNodes, 1> &getR() const { return m_R; }
  T getDisplacement(Index node) const {
    assert(node >= 0 && node < nNodes && "Invalid node index");
    return m_U(node, 0);
  }

  // Setters
  void setEA(T EA) { m_EA = EA; }
  void setDistributedLoad(std::function<T(T)> f) {
    m_rhsFunction = f;
  } // f(x) != const
  void setDistributedLoad(T constant_load) {
    m_rhsFunction = [constant_load](T x) {
      return constant_load;
    }; // f(x) = const
  }
  void setAnalyticSolution(std::function<T(T)> sol) {
    m_analyticSolution = sol;
  }

  // FEM Assembly
  void assembleKF() {
    // Building global matrix via stiffness matrix: need connectitivty matrix
    // C in order to assemble (2x2) into (nNodes x nNodes)
    // Avoiding use of connectivity matrix C for efficiency
    m_K.resetZero();
    m_F.resetZero();
    for (Index e{0}; e < m_mesh.getNumElements(); ++e) {
      // Elementary stiffness matrix assembly
      T h = m_mesh.getLengthEle(e);
      T k = m_EA / h;
      m_K(e, e) += k;
      m_K(e, e + 1) += -k;
      m_K(e + 1, e) += -k;
      m_K(e + 1, e + 1) += k;

      // Elementary force matrix assembly: F_e = int_element N^T * f(x) dx
      auto [x1, x2] = m_mesh.getElementNodes(e);
      auto integrand0 = [=, this](T xi) {
        return m_rhsFunction(physicCoor(xi, x1, x2)) * N1_r(xi);
      };
      auto integrand1 = [=, this](T xi) {
        return m_rhsFunction(physicCoor(xi, x1, x2)) * N2_r(xi);
      };
      m_F[e] += integrationGauss1D_ref(x1, x2, integrand0, 2);
      m_F[e + 1] += integrationGauss1D_ref(x1, x2, integrand1, 2);
    }
    // K is singular before BC applied (no constraints)
  }

  // Boundary Conditions
  void applyNeumanCondition(Index node, T value) {
    // Neumann's condition: EA * u'(x) = value
    assert(node >= 0 && node < nNodes && "Invalid node index");
    m_F[node] += value;
  }

  // Essential BC: at least 1 condition required
  void applyDirichletCondition(Index node, T value) {
    // Back up to calculate reactions later
    m_K_original = m_K;
    m_F_original = m_F;

    // Direct/Elimination method
    m_F[node] = value;
    for (Index i = 0; i < nNodes; ++i) {
      if (i == node)
        continue;
      m_F[i] -= value * m_K(i, node);
    }
    for (Index i = 0; i < nNodes; ++i) {
      m_K(i, node) = 0.0;
    }
    for (Index j = 0; j < nNodes; ++j) {
      m_K(node, j) = 0.0;
    }
    m_K(node, node) = 1.0;
  }

  // Solve
  void solveFEM() {
    m_U = solveLinearSystem(m_K, m_F);
    std::cout << "Displacement vector U:\n";
    for (Index i = 0; i < nNodes; ++i) {
      std::cout << "  u[" << i << "] = " << std::fixed << std::setprecision(6)
                << m_U(i, 0) << '\n';
    }
  }

  // Post-processing
  void calculateReaction() {
    // Compute reaction vector R = K_original * U - F_original
    m_R = -(m_F_original - m_K_original * m_U);
    std::cout << "Reaction force vector R: \n" << Vector<T>(m_R) << std::endl;
  }

  void compareAnalytic() {
    if (!m_analyticSolution) {
      std::cout << "\nNo analytical solution provided.\n";
      return;
    }

    std::cout << "\n=== Comparison with Analytical Solution ===\n";
    std::cout << std::setw(8) << "Node" << std::setw(12) << "x" << std::setw(15)
              << "u_FEM" << std::setw(15) << "u_Exact" << std::setw(15)
              << "Error" << '\n';
    std::cout << std::string(65, '-') << '\n';

    T max_error = 0.0;
    T sum_sq_error = 0.0;

    for (Index i = 0; i < nNodes; ++i) {
      T x = m_mesh.nodeCoords()[i];
      T u_fem = m_U(i, 0);
      T u_exact = m_analyticSolution(x);
      T error = std::abs(u_fem - u_exact);

      max_error = std::max(max_error, error);
      sum_sq_error += error * error;

      std::cout << std::setw(8) << i << std::setw(12) << std::fixed
                << std::setprecision(4) << x << std::setw(15)
                << std::setprecision(6) << u_fem << std::setw(15) << u_exact
                << std::setw(15) << std::scientific << error << '\n';
    }

    T rms_error = std::sqrt(sum_sq_error / nNodes);
    std::cout << "\nMax Error: " << std::scientific << max_error << '\n';
    std::cout << "RMS Error: " << rms_error << '\n';
  }

  void exportResult(const std::string &filename = "fem1D_results.txt") {
    std::cout << "\n=== Exporting Results ===\n";

    // Compute element axial forces: N_e = EA * (u[j] - u[i]) / h_e
    // Store element forces, then assign to nodes
    Vector<T> element_forces(m_mesh.getNumElements());
    for (Index e = 0; e < m_mesh.getNumElements(); ++e) {
      T h = m_mesh.getLengthEle(e);
      element_forces[e] = m_EA * (m_U[e + 1] - m_U[e]) / h;
    }

    // Assign element forces to nodes
    Vector<T> axial_force(nNodes);
    axial_force[0] = element_forces[0]; // First node: force from first element
    for (Index i = 1; i < nNodes - 1; ++i) {
      // Internal nodes: average of adjacent elements (equilibrium)
      axial_force[i] = 0.5 * (element_forces[i - 1] + element_forces[i]);
    }
    axial_force[nNodes - 1] =
        element_forces[m_mesh.getNumElements() - 1]; // Last node: last element

    // Export to file
    std::ofstream file(filename);
    if (!file.is_open()) {
      std::cerr << "Error: Cannot open file " << filename << '\n';
      return;
    }

    file << "# FEM 1D Results\n";
    file << "# EA = " << m_EA << ", L = " << m_length << '\n';
    file << "# Nodes = " << nNodes << ", Elements = " << nNodes - 1 << '\n';
    file << "# Columns: Node, x, u(x), N(x)\n";
    file << std::setw(8) << "Node" << std::setw(15) << "x" << std::setw(18)
         << "u(x)" << std::setw(18) << "N(x)" << '\n';

    for (Index i = 0; i < nNodes; ++i) {
      T x = m_mesh.nodeCoords()[i];
      file << std::setw(8) << i << std::setw(15) << std::fixed
           << std::setprecision(6) << x << std::setw(18) << std::setprecision(8)
           << m_U[i] << std::setw(18) << axial_force[i] << '\n';
    }

    file.close();
    std::cout << "Results exported to: " << filename << '\n';
  }
};

#endif // FEM_H
