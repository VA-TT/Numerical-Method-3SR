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

// Rigidity kt of each bar
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

    for (Index b{0}; b < nBars; ++b) {
      // compute bar's vectors from node coordinates
      vectorBars[b] = nodes[barEnd[b]] - nodes[barOrigin[b]];
      lengthBars[b] = magnitude(vectorBars[b]);
      // Saving neutral bar's lengths
      if (iteration == 0) {
        length0Bars[b] = lengthBars[b];
      }
      unitVectorBars[b] = vectorBars[b] / lengthBars[b];
    }

    // COMPUTING FORCE VECTOR F
    // Internal force calculation: N[b] = A[b] * e[b] = At[b] * r[b]
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

    // Flatting vector F
    Vector<double> externalForceFlatten{flatten(externalForce, nNodes, d)};
    Vector<double> internalForceFlatten{flatten(internalForceNodes, nNodes, d)};
    // Modify the right-hand side: F = F_external - F_internal
    Vector<double> forceF = externalForceFlatten - internalForceFlatten;

    // Saving
    double residualNorm = magnitude(forceF);

    if (residualNorm < epsilon) {
      std::cout << "Solution founded at iteration " << iteration << ".\n\n";

      // Final result
      std::cout << "Final displacement: " << totalDispalcement << std::endl;
      std::cout << "Final positions:\n";
      for (Index n = 0; n < nNodes; ++n) {
        std::cout << "  Node " << n << ": " << nodes[n] << std::endl;
      }
      std::cout << "Total iterations: " << iteration << std::endl;
      break;
    }

    // COMPUTING STIFNESS MATRIX K
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
    assert(approximatelyEqualAbsRel(det(assemblyStiffnessK), 0.0));

    // Penalization method (encastree)
    // Stifness matrix and force adjustment: K' and F'
    assert(nodeImposed.size() == displacementImposed.size() &&
           "Size of imposed nodes must be consistent!");

    double penalty{1.0 / 1e-10}; // 1/epsilon
    int diagonalIndex{};

    // Apply penalization for each imposed node.
    for (Index idx = 0; idx < nodeImposed.size(); ++idx) {
      Index n = nodeImposed[idx];
      for (int k = 0; k < d; ++k) {
        diagonalIndex = d * n + k;
        assemblyStiffnessK(diagonalIndex, diagonalIndex) += penalty;
        forceF[diagonalIndex] = penalty * displacementImposed[idx][k];
      }
    }

    // Solve the linear system to find displacement
    Vector<double> deltaU{solveLinearSystem(assemblyStiffnessK, forceF)};
    double incrementNorm = magnitude(deltaU);

    // IN RA THÔNG TIN HỘI TỤ GIỐNG CODE 2D
    std::cout << "Iteration " << iteration << ": |Fk| = " << residualNorm
              << ": |dx| = " << incrementNorm << std::endl;

    totalDispalcement += deltaU;

    // Saving the output
    iteration_array.push_back(iteration);
    deltaU_array.push_back(deltaU);

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

    // Update the position
    nodes += unflatten(deltaU, nNodes, d);
    // Hardening behavoir(to be implemented)

    iteration++;
  }
  std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";
  return 0;
}