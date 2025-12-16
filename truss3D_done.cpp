#include "library/DualDiffrentiation.h"
#include "library/Matrix.h" //Approximative Comparsion
#include "library/Vector.h"
#include "library/clock.h"
#include "vector"
#include <cassert> // for assert
#include <fstream> //working with files
#include <iomanip> //tab
#include <iostream>
#include <numbers>     // for std::numbers::pi
#include <stdexcept>   //throw exception
#include <type_traits> // precision

// Model
#if 0
 4      5
 o------o
 | \    |
 | 2 \  | 3
 o----\-o
 | \    |
 |   \  |
 o 0   \o 1  
///    ///
#endif

// Small deformation truss

// Model Parameters
namespace modelParameters {
// Problem dimension:
constexpr Index d{3};

// Unit vectors
Vector<double> i1{1.0, 0.0, 0.0};
Vector<double> i2{0.0, 1.0, 0.0};
Vector<double> i3{0.0, 0.0, 1.0};

// Geometry of the truss
constexpr Index nNodes{8}; // number of nodes
constexpr Index nBars{10}; // number of bars

Vector<Vector<double>> nodes{
    {0.0, 0.0, 0.0},     // Node 0
    {10.0, 0.0, 0.0},    // Node 1
    {0.0, 0.0, 10.0},    // Node 2
    {10.0, 0.0, 10.0},   // Node 3
    {0.0, -10.0, 0.0},   // Node 4
    {10.0, -10.0, 0.0},  // Node 5
    {0.0, -10.0, 10.0},  // Node 6
    {10.0, -10.0, 10.0}, // Node 7

};
// Bar connectivity: store node indices for each bar's origin and end
Vector<Index> barOrigin{0, 2, 3, 2, 3, 4, 6, 7, 3, 2};
Vector<Index> barEnd{2, 3, 1, 6, 7, 6, 7, 5, 6, 7};
Vector<Vector<double>> vectorBars(nBars), unitVectorBars(nBars);
Vector<double> lengthBars(nBars);
Vector<double> length0Bars(nBars);
Vector<Index> nodeImposed{0, 1, 4, 5};
Vector<Vector<double>> displacementImposed{
    {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
Vector<Index> nodeFree(nNodes - nodeImposed.size());

// Section dimension
double youngModulus{25e9}; // Module Young
double b{0.02}, h{0.05};
double A{b * h};
double alpha{youngModulus * A};

//  External loads at nodes - Direct 3D force vectors
// External Force applying
Vector<Vector<double>> externalForce{
    {0.0, 0.0, 0.0},      // Node 0: no force
    {0.0, 0.0, 0.0},      // Node 1: no force
    {0.0, 0.0, -15000.0}, // Node 2: 15kN downward (z-direction)
    {0.0, 0.0, -15000.0}, // Node 3: 15kN downward
    {0.0, 0.0, 0.0},      // Node 4: no force
    {0.0, 0.0, 0.0},      // Node 5: no force
    {0.0, 0.0, -15000.0}, // Node 6: 15kN downward
    {0.0, 0.0, -15000.0}  // Node 7: 15kN downward
};
// Displacement Vector
Vector<double> U(nNodes *d), UImposed(nNodes *d), UFree(nNodes *d);
Vector<double> totalDispalcement(nNodes *d);
// Axial Force inside the Bars and applying on Nodes
Vector<Vector<double>> internalForceBar(nBars);
Vector<Vector<double>> internalForceNodes(nNodes);
std::vector<int> iteration_array;
std::vector<Vector<double>> deltaU_array;
std::vector<double> incrementNorm_array;
std::vector<double> residualNorm_array;
Vector<double> k(nBars);
// Choosing non-linear Saint-Venant Kirchhoff law
int law{1};

// tolerance for Newton's method
double epsilon{1e-8};
int max_iteration{100};
} // namespace modelParameters

// Constitutive law
template <typename T>
inline auto constitutiveLaw(int law, double alpha, T l, double l0) {
  switch (law) {
  case 0:
    return (alpha * (l - l0) / l0); // Linear law
  case 1:
    return (alpha * (l * l - l0 * l0) /
            (2 * l0 * l0)); // Saint-Venant Kirchhoff law
  case 2:
    return (alpha * log(l / l0)); // Logarithmic law
  default:
    return (alpha * (l - l0) / l0); // Linear law by default
  }
}

void validateInput() {
  using namespace modelParameters;
  assert(nodes.size() == nNodes && "numbers of nodes must be consistent!");
  assert(barOrigin.size() == nBars && barEnd.size() == nBars &&
         "numbers of bars must be consistent!");

  // Verify external forces
  assert(externalForce.size() == nNodes &&
         "Number of external forces must match number of nodes");
  for (Index i = 0; i < nNodes; ++i) {
    assert(externalForce[i].size() == d &&
           "Each force vector must have d components");
  }
}

void setUp() {
  using namespace modelParameters;
  // Identify free nodes space Nf
  Index nf{0};
  bool isImposed = false;
  for (Index n{0}; n < nNodes; ++n) {
    isImposed = false;
    for (Index k{0}; k < nodeImposed.size(); ++k) {
      if (n == nodeImposed[k]) {
        isImposed = true;
        break;
      }
    }
    if (!isImposed) {
      nodeFree[nf] = n;
      ++nf;
    }
  }
  assert(nf == nodeFree.size() && "Mismatch in free node count");
  iteration_array.reserve(max_iteration);
  deltaU_array.reserve(max_iteration);
}

int main() {
  Timer t;
  using namespace modelParameters;
  validateInput();
  setUp();
  // At(l) = N(l)/l. We wrap l0 via the lambda capture when differentiating.
  auto At = [=](double l, double l0) {
    return constitutiveLaw(law, alpha, l, l0) / l;
  };
  // kt = dA/dl as a function
  auto kt = [=](double l, double l0) {
    auto func = [=](auto x) { return constitutiveLaw(law, alpha, x, l0) / x; };
    return automaticDiff(func, l);
  };

  int iteration{0};

  std::cout << "=== 3D TRUSS PROBLEM ===" << std::endl;

  while (iteration < max_iteration) {

    // ===== STEP 1: Update bar geometry based on current node positions =====
    for (Index b{0}; b < nBars; ++b) {
      vectorBars[b] = nodes[barEnd[b]] - nodes[barOrigin[b]];
      lengthBars[b] = magnitude(vectorBars[b]);
      // Saving neutral bar's lengths (only at iteration 0)
      if (iteration == 0) {
        length0Bars[b] = lengthBars[b];
      }
      unitVectorBars[b] = vectorBars[b] / lengthBars[b];
    }

    // ===== STEP 2: Compute internal forces at current configuration =====
    for (Index b{0}; b < nBars; ++b) {
      internalForceBar[b] =
          constitutiveLaw(law, alpha, lengthBars[b], length0Bars[b]) *
          unitVectorBars[b];
    }
    // Assemble internal force vector at nodes
    for (Index n{0}; n < nNodes; ++n) {
      internalForceNodes[n] = Vector<double>(d);
      for (Index k = 0; k < d; ++k)
        internalForceNodes[n][k] = 0.0;
    }
    for (Index b{0}; b < nBars; ++b) {
      internalForceNodes[barEnd[b]] += internalForceBar[b];
      internalForceNodes[barOrigin[b]] += -internalForceBar[b];
    }

    // Flatten force vectors
    Vector<double> externalForceFlatten{flatten(externalForce, nNodes, d)};
    Vector<double> internalForceFlatten{flatten(internalForceNodes, nNodes, d)};

    // ===== STEP 3: Compute residual force =====
    Vector<double> residualForce = externalForceFlatten - internalForceFlatten;
    double residualNorm = magnitude(residualForce);

    // Diagnostics: residual on free DOFs (L2), L1 norm, and resultant force
    double residualFreeL2 = 0.0, residualL1 = 0.0;
    Vector<double> resultant(d);
    for (Index k = 0; k < d; ++k)
      resultant[k] = 0.0;

    // Sum resultant and L1 over all DOFs
    for (Index j = 0; j < residualForce.size(); ++j) {
      residualL1 += std::abs(residualForce[j]);
    }
    // Resultant by summing per node components
    for (Index n = 0; n < nNodes; ++n) {
      for (Index k = 0; k < d; ++k) {
        resultant[k] += residualForce[d * n + k];
      }
    }
    // L2 over free DOFs only
    for (Index ii = 0; ii < nodeFree.size(); ++ii) {
      Index n = nodeFree[ii];
      for (Index k = 0; k < d; ++k) {
        double v = residualForce[d * n + k];
        residualFreeL2 += v * v;
      }
    }
    residualFreeL2 = std::sqrt(residualFreeL2);
    double resultantNorm = magnitude(resultant);

    std::cout << "Iteration " << iteration
              << ": |Fk|_L2(all) = " << residualNorm
              << ", |Fk|_L2(free) = " << residualFreeL2
              << ", |Fk|_L1 = " << residualL1
              << ", |Resultant| = " << resultantNorm << std::endl;

    // Print component-wise global force (first 12 DOFs for brevity)
    std::cout << "  F_global (first 12 DOFs): [";
    for (Index i = 0; i < std::min((Index)12, (Index)residualForce.size());
         ++i) {
      std::cout << residualForce[i];
      if (i < std::min((Index)11, (Index)residualForce.size() - 1))
        std::cout << ", ";
    }
    std::cout << "]\n";

    // ===== STEP 4: Check convergence =====
    if (residualNorm < epsilon) {
      std::cout << "Solution converged at iteration " << iteration << ".\n\n";
      std::cout << "Final displacement: " << totalDispalcement << std::endl;
      std::cout << "Final positions:\n";
      for (Index n = 0; n < nNodes; ++n) {
        std::cout << "  Node " << n << ": " << nodes[n] << std::endl;
      }
      std::cout << "Total iterations: " << iteration << std::endl;
      break;
    }

    // If residual is converging slowly but increments are tiny, we can also
    // accept convergence based on increment norm alone:
    // if (iteration > 1 && incrementNorm < epsilon * 1e-5) {

    // ===== STEP 5: Assemble stiffness matrix =====
    // kb
    for (Index b{0}; b < nBars; ++b) {
      k[b] = kt(lengthBars[b], length0Bars[b]);
    }

    Vector<Matrix<double, d, d>> elementaryApplicationK(nBars);
    for (Index b{0}; b < nBars; ++b) {
      elementaryApplicationK[b] =
          At(lengthBars[b], length0Bars[b]) * Matrix<double, d, d>::identity() +
          (1 / lengthBars[b]) * k[b] *
              tensorProduct<d, d>(vectorBars[b], vectorBars[b]);
    }
    //  Connectivity Matrix
    Vector<Matrix<double, d, d * nNodes>> C(nNodes);

    const auto Id = Matrix<double, d, d>::identity();
    Matrix<double, d, d * nNodes> Ci{};
    for (Index n{0}; n < nNodes; ++n) {
      Ci = Matrix<double, d, d * nNodes>{};
      for (Index i = 0; i < d; ++i) {
        for (Index j = 0; j < d; ++j) {
          Ci(i, n * d + j) = Id(i, j);
        }
      }
      C[n] = Ci;
    }

    Matrix<double, d * nNodes, d * nNodes> assemblyStiffnessK{};

    for (Index b{0}; b < nBars; ++b) {
      assemblyStiffnessK += (C[barEnd[b]] - C[barOrigin[b]]).transpose() *
                            elementaryApplicationK[b] *
                            (C[barEnd[b]] - C[barOrigin[b]]);
    }
    // Check the singularity of stiffness matrix K
    // Note: If matrix is singular, the penalty method should regularize it
    // double detK = det(assemblyStiffnessK);
    // if (approximatelyEqualAbsRel(detK, 0.0)) {
    //   std::cout << "Warning: Stiffness matrix is singular before
    //   penalization\n";
    // }

    // Penalization method (encastree)
    // Apply penalty to enforce displacement constraints at imposed nodes
    // The penalty enforces: u_imposed = displacementImposed
    assert(nodeImposed.size() == displacementImposed.size() &&
           "Size of imposed nodes must be consistent!");

    double penalty{1e10}; // Large penalty parameter
    int diagonalIndex{};

    // Copy residual to forceF (will be modified by penalization)
    Vector<double> forceF = residualForce;

    // Apply penalty: Add penalty * (u - u_imposed) term
    // This modifies K[i,i] += penalty and F[i] += penalty * u_imposed
    for (Index idx = 0; idx < nodeImposed.size(); ++idx) {
      Index n = nodeImposed[idx];
      for (int k = 0; k < d; ++k) {
        diagonalIndex = d * n + k;
        assemblyStiffnessK(diagonalIndex, diagonalIndex) += penalty;
        // RHS: F = F - K*u_current + penalty*u_imposed
        // Simplified: forceF[i] = F[i] + penalty*(u_imposed[i] - u_current[i])
        forceF[diagonalIndex] = residualForce[diagonalIndex] +
                                penalty * (displacementImposed[idx][k] -
                                           totalDispalcement[diagonalIndex]);
      }
    }

    // Solve the linear system to find displacement
    Vector<double> deltaU{solveLinearSystem(assemblyStiffnessK, forceF)};
    double incrementNorm = magnitude(deltaU);

    std::cout << "        |dx| = " << incrementNorm << std::endl;

    // Print component-wise dx (first 12 DOFs for brevity)
    std::cout << "  dx (first 12 DOFs): [";
    for (Index i = 0; i < std::min((Index)12, (Index)deltaU.size()); ++i) {
      std::cout << deltaU[i];
      if (i < std::min((Index)11, (Index)deltaU.size() - 1))
        std::cout << ", ";
    }
    std::cout << "]\n";

    // ===== STEP 7: Update node positions and total displacement =====
    nodes += unflatten(deltaU, nNodes, d);
    totalDispalcement += deltaU;

    // Convergence check: if increment norm is below tolerance, solution has
    // converged
    if (iteration > 0 && incrementNorm < epsilon * 1e-2) {
      // Recompute internal forces at the converged configuration
      for (Index b{0}; b < nBars; ++b) {
        vectorBars[b] = nodes[barEnd[b]] - nodes[barOrigin[b]];
        lengthBars[b] = magnitude(vectorBars[b]);
        unitVectorBars[b] = vectorBars[b] / lengthBars[b];
      }
      for (Index b{0}; b < nBars; ++b) {
        internalForceBar[b] =
            constitutiveLaw(law, alpha, lengthBars[b], length0Bars[b]) *
            unitVectorBars[b];
      }
      for (Index n{0}; n < nNodes; ++n) {
        internalForceNodes[n] = Vector<double>(d);
        for (Index k = 0; k < d; ++k)
          internalForceNodes[n][k] = 0.0;
      }
      for (Index b{0}; b < nBars; ++b) {
        internalForceNodes[barEnd[b]] += internalForceBar[b];
        internalForceNodes[barOrigin[b]] += -internalForceBar[b];
      }
      Vector<double> internalForceFlattenFinal{
          flatten(internalForceNodes, nNodes, d)};
      Vector<double> externalForceFlattenFinal{
          flatten(externalForce, nNodes, d)};
      Vector<double> residualForceFinal =
          externalForceFlattenFinal - internalForceFlattenFinal;

      // Calculate physical reactions at imposed nodes
      // Reaction = Internal force - External force at supports (to maintain
      // equilibrium) Or equivalently: Reaction = -Residual at imposed nodes
      Vector<Vector<double>> reactionR(nNodes);
      for (Index ii = 0; ii < nNodes; ++ii) {
        reactionR[ii] = Vector<double>(d);
        for (Index k = 0; k < d; ++k)
          reactionR[ii][k] = 0.0;
      }
      for (Index idx = 0; idx < nodeImposed.size(); ++idx) {
        Index n = nodeImposed[idx];
        for (Index k = 0; k < d; ++k) {
          int diagIdx = d * n + k;
          // Physical reaction = internal - external at support
          reactionR[n][k] = internalForceFlattenFinal[diagIdx] -
                            externalForceFlattenFinal[diagIdx];
        }
      }
      Vector<double> reactionFlatten = flatten(reactionR, nNodes, d);

      // Print residual per free DOF
      double freeL2 = 0.0;
      std::cout << "Residual on free DOFs (final):\n";
      for (Index ii = 0; ii < nodeFree.size(); ++ii) {
        Index n = nodeFree[ii];
        double rx = residualForceFinal[d * n + 0];
        double ry = residualForceFinal[d * n + 1];
        double rz = residualForceFinal[d * n + 2];
        freeL2 += rx * rx + ry * ry + rz * rz;
        std::cout << "  node " << n << ": (" << rx << ", " << ry << ", " << rz
                  << ")\n";
      }
      freeL2 = std::sqrt(freeL2);
      std::cout << "  ||residual||_L2(free) = " << freeL2 << "\n\n";

      // Check equilibrium on imposed DOFs: residual + reaction should be ~0
      double imposedEquilibriumL2 = 0.0;
      std::cout << "Equilibrium check on imposed DOFs (Residual + Reaction):\n";
      for (Index idx = 0; idx < nodeImposed.size(); ++idx) {
        Index n = nodeImposed[idx];
        double ex = residualForceFinal[d * n + 0] + reactionFlatten[d * n + 0];
        double ey = residualForceFinal[d * n + 1] + reactionFlatten[d * n + 1];
        double ez = residualForceFinal[d * n + 2] + reactionFlatten[d * n + 2];
        imposedEquilibriumL2 += ex * ex + ey * ey + ez * ez;
        std::cout << "  node " << n << ": residual + reaction = (" << ex << ", "
                  << ey << ", " << ez << ")\n";
      }
      imposedEquilibriumL2 = std::sqrt(imposedEquilibriumL2);
      std::cout << "  ||residual + reaction||_L2(imposed) = "
                << imposedEquilibriumL2 << "\n\n";

      // Global equilibrium: sum of all (F_external - F_internal + R) should be
      // ~0
      Vector<double> globalEquilibrium = residualForceFinal + reactionFlatten;
      double globalEquilibriumNorm = magnitude(globalEquilibrium);
      std::cout
          << "Global equilibrium check (Residual + Reaction over all DOFs):\n";
      std::cout << "  ||F_ext - F_int + R||_L2 = " << globalEquilibriumNorm
                << "\n\n";

      std::cout << "Solution converged (increment < tolerance) at iteration "
                << iteration << ".\n\n";
      std::cout << "Final displacement: " << totalDispalcement << std::endl;
      std::cout << "Final positions:\n";
      for (Index n = 0; n < nNodes; ++n) {
        std::cout << "  Node " << n << ": " << nodes[n] << std::endl;
      }
      std::cout << "Total iterations: " << iteration << std::endl;
      break;
    }

    // Saving the output
    iteration_array.push_back(iteration);
    deltaU_array.push_back(deltaU);
    incrementNorm_array.push_back(incrementNorm);
    residualNorm_array.push_back(residualNorm);

    // Calculate the reaction
    Vector<Vector<double>> reactionR(nNodes);
    for (Index ii = 0; ii < nNodes; ++ii) {
      reactionR[ii] = Vector<double>(d);
      for (Index k = 0; k < d; ++k)
        reactionR[ii][k] = 0.0;
    }
    for (Index idx = 0; idx < nodeImposed.size(); ++idx) {
      Index n = nodeImposed[idx];
      for (Index k = 0; k < d; ++k) {
        diagonalIndex = d * n + k;
        reactionR[n][k] = penalty * (totalDispalcement[diagonalIndex] -
                                     displacementImposed[idx][k]);
      }
    }

    // Hardening behavoir(to be implemented)

    iteration++;
  }

  // Write final node positions to file
  std::ofstream nodeFile("truss3D_node.dat");
  if (nodeFile.is_open()) {
    nodeFile << "# Node positions after deformation (x, y, z)\n";
    nodeFile << "# Node_ID x y z\n";
    for (Index n = 0; n < nNodes; ++n) {
      nodeFile << n << " " << nodes[n][0] << " " << nodes[n][1] << " "
               << nodes[n][2] << "\n";
    }
    nodeFile.close();
    std::cout << "\nFinal node positions written to truss3D_node.dat\n";
  } else {
    std::cerr << "Error: Could not open truss3D_node.dat for writing\n";
  }

  // Write iteration convergence data to report directory
  std::ofstream plotFile("report/truss3D_plot.dat");
  if (plotFile.is_open()) {
    plotFile << "# Iteration |dx| |F_global|\n";
    for (size_t i = 0; i < iteration_array.size(); ++i) {
      plotFile << iteration_array[i] << " " << incrementNorm_array[i] << " "
               << residualNorm_array[i] << "\n";
    }
    plotFile.close();
    std::cout << "Convergence data written to report/truss3D_plot.dat\n";
  } else {
    std::cerr << "Error: Could not open report/truss3D_plot.dat for writing\n";
  }

  std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";
  return 0;
}