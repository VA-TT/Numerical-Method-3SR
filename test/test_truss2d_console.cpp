#include "../library/DualDiffrentiation.h"
#include "../library/Matrix.h"
#include "../library/Vector.h"
#include "../library/clock.h"
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <vector>

// Model giống console2D
#if 0
|o-----o
      / |\
 ^i2 /  |  \> F
 |  /    theta
|o/---> i1
#endif

namespace modelParameters {
// Dimension
constexpr Index d{2};

// Geometry (giống console2D)
constexpr Index nNodes{
    3}; // 3 nodes: 0 (origin), 1 (top-left), 2 (C - free node)
constexpr Index nBars{2}; // 2 bars

double a{10.0}; // Distance

// Initial node positions
Vector<Vector<double>> nodes{
    {0.0, 0.0}, // Node 0: origin (fixed)
    {0.0, a},   // Node 1: top-left (fixed)
    {a, a}      // Node 2: C (initial position, will move)
};

// Bar connectivity
Vector<Index> barOrigin{0, 1}; // Bar 0: 0->2, Bar 1: 1->2
Vector<Index> barEnd{2, 2};

Vector<Vector<double>> vectorBars(nBars), unitVectorBars(nBars);
Vector<double> lengthBars(nBars), length0Bars(nBars);

// Fixed nodes (0 and 1)
Vector<Index> nodeImposed{0, 1};
Vector<Vector<double>> displacementImposed{{0.0, 0.0}, {0.0, 0.0}};
Vector<Index> nodeFree{2}; // Only node 2 is free

// Material properties (giống console2D)
double E{25e9};
double b1{0.02}, h1{0.05};
double b2{0.02}, h2{0.05};
double A1{b1 * h1}, A2{b2 * h2};
double alpha1{E * A1}, alpha2{E * A2};

// External force (giống console2D)
double force{1.5e6};
double theta{20.0};
double thetaRadian{theta * std::numbers::pi / 180.0};

Vector<Vector<double>> externalForce{
    {0.0, 0.0}, // Node 0: no force
    {0.0, 0.0}, // Node 1: no force
    {force * std::sin(thetaRadian), -force *std::cos(thetaRadian)} // Node 2
};

// State vectors
Vector<double> totalDisplacement(nNodes *d);
Vector<Vector<double>> internalForceBar(nBars);
Vector<Vector<double>> internalForceNodes(nNodes);
Vector<Vector<double>> reactionR(nNodes);

// Convergence
int law{1}; // SVK like console2D
double penalty{1e25};
double epsilon{1e-8};
int max_iteration{100};

std::vector<int> iteration_array;
std::vector<Vector<double>> deltaU_array;
std::vector<double> incrementNorm_array;
std::vector<double> residualNorm_array;
} // namespace modelParameters

// Constitutive law (heterogeneous bars with different alpha)
template <typename T>
inline auto constitutiveLaw(int law, double alpha, T l, double l0) {
  switch (law) {
  case 0:
    return (alpha * (l - l0) / l0);
  case 1:
    return (alpha * (l * l - l0 * l0) / (2 * l0 * l0)); // SVK
  case 2:
    return (alpha * log(l / l0));
  default:
    return (alpha * (l - l0) / l0);
  }
}

void validateInput() {
  using namespace modelParameters;
  assert(nodes.size() == nNodes);
  assert(barOrigin.size() == nBars && barEnd.size() == nBars);
  assert(externalForce.size() == nNodes);
  assert(nodeImposed.size() == displacementImposed.size());
}

void setUp() {
  using namespace modelParameters;
  iteration_array.reserve(max_iteration);
  deltaU_array.reserve(max_iteration);

  // Initialize total displacement to zero
  for (Index i = 0; i < nNodes * d; ++i) {
    totalDisplacement[i] = 0.0;
  }

  // Initialize internal force nodes with correct size
  for (Index n = 0; n < nNodes; ++n) {
    internalForceNodes[n] = Vector<double>::zero(d);
  }

  // Initialize reaction vectors with correct size
  for (Index n = 0; n < nNodes; ++n) {
    reactionR[n] = Vector<double>::zero(d);
  }

  // Initialize bar vectors
  for (Index b = 0; b < nBars; ++b) {
    vectorBars[b] = Vector<double>::zero(d);
    unitVectorBars[b] = Vector<double>::zero(d);
    internalForceBar[b] = Vector<double>::zero(d);
  }
}

int main() {
  Timer t;
  using namespace modelParameters;

  validateInput();
  setUp();

  // Lambda functions for each bar (different alpha)
  auto At = [=](Index barIdx, double l, double l0) {
    double alpha = (barIdx == 0) ? alpha1 : alpha2;
    return constitutiveLaw(law, alpha, l, l0) / l;
  };

  auto kt = [=](Index barIdx, double l, double l0) {
    double alpha = (barIdx == 0) ? alpha1 : alpha2;
    auto func = [=](auto x) { return constitutiveLaw(law, alpha, x, l0) / x; };
    return automaticDiff(func, l);
  };

  int iteration{0};

  std::cout << "=== 2D TRUSS TEST (Console2D input, Truss frame) ==="
            << std::endl;
  std::cout << "Geometry: 3 nodes, 2 bars" << std::endl;
  std::cout << "Material: E=" << E << ", A1=" << A1 << ", A2=" << A2
            << std::endl;
  std::cout << "Force: F=" << force << "N at theta=" << theta << " degrees"
            << std::endl;
  std::cout << "Law: " << law << " (SVK)\n" << std::endl;

  while (iteration < max_iteration) {
    // Update bar geometry
    for (Index b{0}; b < nBars; ++b) {
      vectorBars[b] = nodes[barEnd[b]] - nodes[barOrigin[b]];
      lengthBars[b] = magnitude(vectorBars[b]);
      if (iteration == 0) {
        length0Bars[b] = lengthBars[b];
      }
      unitVectorBars[b] = vectorBars[b] / lengthBars[b];
    }

    // Compute internal forces (each bar has different alpha)
    for (Index b{0}; b < nBars; ++b) {
      double alpha = (b == 0) ? alpha1 : alpha2;
      internalForceBar[b] =
          constitutiveLaw(law, alpha, lengthBars[b], length0Bars[b]) *
          unitVectorBars[b];
    }

    // Assemble nodal forces
    for (Index n{0}; n < nNodes; ++n) {
      if (internalForceNodes[n].size() == d) {
        internalForceNodes[n].resetZero();
      } else {
        internalForceNodes[n] = Vector<double>::zero(d);
      }
    }
    for (Index b{0}; b < nBars; ++b) {
      internalForceNodes[barEnd[b]] += internalForceBar[b];
      internalForceNodes[barOrigin[b]] += -internalForceBar[b];
    }

    // Flatten forces
    Vector<double> externalForceFlatten{flatten(externalForce, nNodes, d)};
    Vector<double> internalForceFlatten{flatten(internalForceNodes, nNodes, d)};
    Vector<double> residualForce = externalForceFlatten - internalForceFlatten;

    // Calculate reactions
    for (Index ii = 0; ii < nNodes; ++ii) {
      if (reactionR[ii].size() == d) {
        reactionR[ii].resetZero();
      } else {
        reactionR[ii] = Vector<double>::zero(d);
      }
    }
    for (Index idx = 0; idx < nodeImposed.size(); ++idx) {
      Index n = nodeImposed[idx];
      for (Index k = 0; k < d; ++k) {
        reactionR[n][k] =
            internalForceFlatten[d * n + k] - externalForceFlatten[d * n + k];
      }
    }
    Vector<double> reactionFlatten = flatten(reactionR, nNodes, d);

    // Global equilibrium check
    Vector<double> globalEquilibrium = residualForce + reactionFlatten;
    double globalEquilibriumNorm = magnitude(globalEquilibrium);

    if (globalEquilibriumNorm < epsilon) {
      std::cout << "\nSolution converged at iteration " << iteration << ".\n\n";

      std::cout << "=== TRUSS FRAME RESULT ===" << std::endl;
      std::cout << "Final displacement of node 2 (C): ["
                << totalDisplacement[d * 2 + 0] << ", "
                << totalDisplacement[d * 2 + 1] << "]" << std::endl;
      std::cout << "Final position of node 2 (C): " << nodes[2] << std::endl;

      std::cout << "\nReactions at supports:" << std::endl;
      for (Index idx = 0; idx < nodeImposed.size(); ++idx) {
        Index n = nodeImposed[idx];
        std::cout << "  Node " << n << ": " << reactionR[n] << std::endl;
      }

      // Write result for comparison
      std::ofstream outFile("report/truss2d_console_result.dat");
      outFile << std::setprecision(16);
      outFile << "# Truss frame result with console2D input\n";
      outFile << "# displacement_x displacement_y position_x position_y\n";
      outFile << totalDisplacement[d * 2 + 0] << " "
              << totalDisplacement[d * 2 + 1] << " " << nodes[2][0] << " "
              << nodes[2][1] << "\n";
      outFile.close();

      break;
    }

    // Build stiffness matrix
    Vector<double> k(nBars);
    for (Index b{0}; b < nBars; ++b) {
      k[b] = kt(b, lengthBars[b], length0Bars[b]);
    }

    Vector<Matrix<double, d, d>> elementaryApplicationK(nBars);
    for (Index b{0}; b < nBars; ++b) {
      elementaryApplicationK[b] =
          At(b, lengthBars[b], length0Bars[b]) *
              Matrix<double, d, d>::identity() +
          (1 / lengthBars[b]) * k[b] *
              tensorProduct<d, d>(vectorBars[b], vectorBars[b]);
    }

    // Connectivity matrix
    Vector<Matrix<double, d, d * nNodes>> C(nNodes);
    const auto Id = Matrix<double, d, d>::identity();
    for (Index n{0}; n < nNodes; ++n) {
      Matrix<double, d, d * nNodes> Ci{};
      for (Index i = 0; i < d; ++i) {
        for (Index j = 0; j < d; ++j) {
          Ci(i, n * d + j) = Id(i, j);
        }
      }
      C[n] = Ci;
    }

    // Assemble global stiffness
    Matrix<double, d * nNodes, d * nNodes> assemblyStiffnessK{};
    for (Index b{0}; b < nBars; ++b) {
      assemblyStiffnessK += (C[barEnd[b]] - C[barOrigin[b]]).transpose() *
                            elementaryApplicationK[b] *
                            (C[barEnd[b]] - C[barOrigin[b]]);
    }

    // Penalty
    int diagonalIndex{};
    for (Index idx = 0; idx < nodeImposed.size(); ++idx) {
      Index n = nodeImposed[idx];
      for (int k = 0; k < d; ++k) {
        diagonalIndex = d * n + k;
        assemblyStiffnessK(diagonalIndex, diagonalIndex) += penalty;
        residualForce[diagonalIndex] = penalty * displacementImposed[idx][k];
      }
    }

    // Solve
    Vector<double> deltaU{solveLinearSystem(assemblyStiffnessK, residualForce)};
    double incrementNorm = magnitude(deltaU);

    std::cout << "Iteration " << iteration
              << ": |Fk| = " << magnitude(residualForce)
              << ": |F_global| = " << globalEquilibriumNorm
              << ": |dx| = " << incrementNorm << std::endl;

    totalDisplacement += deltaU;
    nodes += unflatten(deltaU, nNodes, d);

    iteration_array.push_back(iteration);
    deltaU_array.push_back(deltaU);
    incrementNorm_array.push_back(incrementNorm);
    residualNorm_array.push_back(globalEquilibriumNorm);

    iteration++;
  }

  if (iteration >= max_iteration) {
    std::cerr << "Warning: Maximum iterations reached!" << std::endl;
  }

  std::cout << "\nTime elapsed: " << t.elapsed() << " seconds\n";

  std::cout << "\n=== Now compare with console2D ===" << std::endl;
  std::cout << "Run: ./console2D" << std::endl;
  std::cout << "Compare displacement of node C (node 2)" << std::endl;

  return 0;
}