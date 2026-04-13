#include "library/DualDifferentiation.h"
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
DynamicVector<double> i1{1.0, 0.0, 0.0};
DynamicVector<double> i2{0.0, 1.0, 0.0};
DynamicVector<double> i3{0.0, 0.0, 1.0};

// Geometry of the truss
constexpr Index nNodes{8}; // number of nodes
constexpr Index nBars{10}; // number of bars

DynamicVector<DynamicVector<double>> nodes{
    {0.0, 0.0, 0.0},   // Node 0
    {10.0, 0.0, 0.0},  // Node 1
    {0.0, 0.0, 10.0},  // Node 2
    {10.0, 0.0, 10.0}, // Node 3

};
// Bar connectivity: store node indices for each bar's origin and end
DynamicVector<Index> barOrigin{0, 2, 3, 2, 3, 4, 6, 7, 3, 2};
DynamicVector<Index> barEnd{2, 3, 1, 6, 7, 6, 7, 5, 6, 7};
DynamicVector<DynamicVector<double>> vectorBars(nBars), unitVectorBars(nBars);
DynamicVector<double> lengthBars(nBars);
DynamicVector<double> length0Bars(nBars);
DynamicVector<Index> nodeImposed{0, 1, 4, 5};
DynamicVector<DynamicVector<double>> displacementImposed{
    {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
DynamicVector<Index> nodeFree(nNodes - nodeImposed.size());

// Section dimension
double youngModulus{25e9}; // Module Young
double b{0.02}, h{0.05};
double A{b * h};
double alpha{youngModulus * A};

// External Force applying
double forceZ{-1.5e6};
DynamicVector<DynamicVector<double>> externalForce{
    {0.0, 0.0, 0.0},    // Node 0: no force
    {0.0, 0.0, 0.0},    // Node 1: no force
    {0.0, 0.0, forceZ}, // Node 2: 1.5 MN downward (z-direction)
    {0.0, 0.0, forceZ}, // Node 3: 1.5 MN downward
    {0.0, 0.0, 0.0},    // Node 4: no force
    {0.0, 0.0, 0.0},    // Node 5: no force
    {0.0, 0.0, forceZ}, // Node 6: 1.5 MN downward
    {0.0, 0.0, forceZ}  // Node 7: 1.5 MN downward
};
// Displacement Vector
DynamicVector<double> U(nNodes *d), UImposed(nNodes *d), UFree(nNodes *d);
DynamicVector<double> totalDispalcement(nNodes *d);
// Axial Force inside the Bars and applying on Nodes
DynamicVector<DynamicVector<double>> internalForceBar(nBars);
DynamicVector<DynamicVector<double>> internalForceNodes(nNodes);
DynamicVector<DynamicVector<double>> reactionR(nNodes);
Matrix<double, d * nNodes, d * nNodes> assemblyStiffnessK{};
std::vector<int> iteration_array;
std::vector<DynamicVector<double>> deltaU_array;
std::vector<double> incrementNorm_array;
std::vector<double> residualNorm_array;
DynamicVector<double> k(nBars);

// Choosing non-linear Saint-Venant Kirchhoff law
int law{1};

// Penalization method
double penalty{1.0 / 1e-18}; // 1/epsilon

// tolerance for Newton's method
double epsilon{1e-10};
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
  assert(nodeImposed.size() == displacementImposed.size() &&
         "Size of imposed nodes must be consistent!");
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

    // Update bar geometry based on current node positions
    for (Index b{0}; b < nBars; ++b) {
      vectorBars[b] = nodes[barEnd[b]] - nodes[barOrigin[b]];
      lengthBars[b] = magnitude(vectorBars[b]);
      // Saving neutral bar's lengths
      if (iteration == 0) {
        length0Bars[b] = lengthBars[b];
      }
      unitVectorBars[b] = vectorBars[b] / lengthBars[b];
    }

    // Compute internal forces at current configuration:
    // N[b] = A[b] * e[b] = At[b] * r[b]
    for (Index b{0}; b < nBars; ++b) {
      internalForceBar[b] =
          constitutiveLaw(law, alpha, lengthBars[b], length0Bars[b]) *
          unitVectorBars[b];
    }

    // Assemble internal force vector at nodes
    for (Index n{0}; n < nNodes; ++n) {
      internalForceNodes[n] =
          DynamicVector<double>::zero(d); // reset to prepare for accumulation
    }
    for (Index b{0}; b < nBars; ++b) {
      internalForceNodes[barEnd[b]] += internalForceBar[b];
      internalForceNodes[barOrigin[b]] += -internalForceBar[b];
    }

    // Flatten force vectors
    DynamicVector<double> externalresidualForcelatten{
        flatten(externalForce, nNodes, d)};
    DynamicVector<double> internalresidualForcelatten{
        flatten(internalForceNodes, nNodes, d)};

    // Compute residual force : F = Fext - Fint
    DynamicVector<double> residualForce =
        externalresidualForcelatten - internalresidualForcelatten;

    // Reaction at fixed points
    for (Index ii = 0; ii < nNodes; ++ii) {
      reactionR[ii] = DynamicVector<double>::zero(d);
    }

    for (Index idx = 0; idx < nodeImposed.size(); ++idx) {
      Index n = nodeImposed[idx];
      for (Index k = 0; k < d; ++k) {
        // R = F_int - F_ext
        reactionR[n][k] = internalresidualForcelatten[d * n + k] -
                          externalresidualForcelatten[d * n + k];
      }
    }
    DynamicVector<double> reactionFlatten = flatten(reactionR, nNodes, d);

    // Stopping criterion: F_ext - F_int + R = 0
    DynamicVector<double> globalEquilibrium = residualForce + reactionFlatten;
    double globalEquilibriumNorm = magnitude(globalEquilibrium);

    // Check global force equilibrium at nodes
    if (magnitude(globalEquilibrium) < epsilon) {
      std::cout << "Solution converged at iteration " << iteration << ".\n\n";

      std::cout << "Final displacement: " << totalDispalcement << std::endl;
      std::cout << "\nFinal positions:\n";
      for (Index n = 0; n < nNodes; ++n) {
        std::cout << "  Node " << n << ": " << nodes[n] << std::endl;
      }
      std::cout << "\nReactions at supports:\n";
      for (Index idx = 0; idx < nodeImposed.size(); ++idx) {
        Index n = nodeImposed[idx];
        std::cout << "  Node " << n << ": " << reactionR[n] << std::endl;
      }
      break;
    }

    //  Assemble stiffness matrix
    // kb
    for (Index b{0}; b < nBars; ++b) {
      k[b] = kt(lengthBars[b], length0Bars[b]);
    }

    DynamicVector<Matrix<double, d, d>> elementaryApplicationK(nBars);
    for (Index b{0}; b < nBars; ++b) {
      elementaryApplicationK[b] =
          At(lengthBars[b], length0Bars[b]) * Matrix<double, d, d>::identity() +
          (1 / lengthBars[b]) * k[b] *
              tensorProduct<d, d>(vectorBars[b], vectorBars[b]);
    }
    //  Connectivity Matrix
    DynamicVector<Matrix<double, d, d * nNodes>> C(nNodes);
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
    assemblyStiffnessK.resetZero(); // reset for next += step

    for (Index b{0}; b < nBars; ++b) {
      assemblyStiffnessK += (C[barEnd[b]] - C[barOrigin[b]]).transpose() *
                            elementaryApplicationK[b] *
                            (C[barEnd[b]] - C[barOrigin[b]]);
    }
    // Check the singularity of stiffness matrix K

    assert(approximatelyEqualAbsRel(det(assemblyStiffnessK), 0.0));
    assert(assemblyStiffnessK.isSymmetric());

    // Penalization method on fixed
    // Stifness matrix and force adjustment: K' and F'
    int diagonalIndex{};
    for (Index idx = 0; idx < nodeImposed.size(); ++idx) {
      Index n = nodeImposed[idx];
      for (int k = 0; k < d; ++k) {
        diagonalIndex = d * n + k;
        assemblyStiffnessK(diagonalIndex, diagonalIndex) += penalty;
        residualForce[diagonalIndex] = penalty * displacementImposed[idx][k];
      }
    }

    // Solve the linear system to find displacement
    DynamicVector<double> deltaU{
        solveLinearSystem(assemblyStiffnessK, residualForce)};
    double incrementNorm = magnitude(deltaU);

    nodes += unflatten(deltaU, nNodes, d);
    totalDispalcement += deltaU;

    // Saving the output
    iteration_array.push_back(iteration);
    deltaU_array.push_back(deltaU);
    incrementNorm_array.push_back(incrementNorm);
    residualNorm_array.push_back(globalEquilibriumNorm); // Saving global force

    // output to console
    std::cout << "Iteration " << iteration
              << ": |Fk| = " << magnitude(globalEquilibrium)
              << ": |dx| = " << incrementNorm << std::endl;

    // Increase the iteration
    iteration++;
  }

  // Write final node positions to file
  std::ofstream nodeFile("report/truss3D_node.dat");

  nodeFile << "# Node positions after deformation (x, y, z)\n";
  nodeFile << "# Node_ID x y z\n";
  for (Index n = 0; n < nNodes; ++n) {
    nodeFile << n << " " << nodes[n][0] << " " << nodes[n][1] << " "
             << nodes[n][2] << "\n";
  }
  nodeFile.close();

  // Write iteration convergence data to report directory
  std::ofstream plotFile("report/truss3D_plot.dat");
  plotFile << "# Iteration |dx| |F_global|\n";
  for (size_t i = 0; i < iteration_array.size(); ++i) {
    plotFile << iteration_array[i] << " " << incrementNorm_array[i] << " "
             << residualNorm_array[i] << "\n";
  }
  plotFile.close();

  std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";
  return 0;
}
