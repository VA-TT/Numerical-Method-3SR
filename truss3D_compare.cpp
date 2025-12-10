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
// Problem dimension: Considering the 2D implementation first
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
// Chuyen ve Matrix
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

template <typename T>
inline auto constitutiveLaw(double alpha, T l, double l0) {
  return (alpha * (l - l0) / l0);
  // return (alpha * (l * l - l0 * l0) / (2 * l0 * l0));
  // return (alpha * std::log(l / l0));
}

// tolerance for Newton's method
double epsilon{1e-8};
int max_iteration{100};
} // namespace modelParameters

int main() {
  Timer t;
  using namespace modelParameters;
  // auto func1 = [](auto x1) { return constitutiveLaw(alpha1, x1, l01); };
  // auto func2 = [](auto x2) { return constitutiveLaw(alpha2, x2, l02); };

  // Checking the input
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

  // Setting up

  Vector<double> U(nNodes * d), UImposed(nNodes * d), UFree(nNodes * d);

  // Identify free nodes
  Vector<Index> nodeFree(nNodes - nodeImposed.size());
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

  int iteration{0};

  std::vector<int> iteration_array;
  std::vector<Vector<double>> deltaU_array;
  iteration_array.reserve(max_iteration);
  deltaU_array.reserve(max_iteration);
  Vector<double> totalDispalcement(nNodes * d);

  std::cout << "=== Newton-Raphson Iteration ===" << std::endl;

  while (iteration < max_iteration) {

    for (Index b{0}; b < nBars; ++b) {
      // compute original bar's vectors from node coordinates
      vectorBars[b] = nodes[barEnd[b]] - nodes[barOrigin[b]];
      lengthBars[b] = magnitude(vectorBars[b]);
      unitVectorBars[b] = vectorBars[b] / lengthBars[b];
    }
    // Rigidity in small deformation configuration: N = k e (u2 - u1)
    Vector<double> k(nBars);

    // Fix the constitutive law
    for (Index b{0}; b < nBars; ++b) {
      k[b] = youngModulus * A / lengthBars[b];
    }

    Vector<Matrix<double, d, d>> elementaryApplicationK(nBars);
    for (Index b{0}; b < nBars; ++b) {
      elementaryApplicationK[b] =
          k[b] * tensorProduct<d, d>(unitVectorBars[b], unitVectorBars[b]);
    }
    //   std::cout << elementaryApplicationK;

    // Matrix<double, 2, 2> connectivityMatrix{1.0, -1.0, -1.0, 1.0};
    // Vector<Matrix<double, d * 2, d * 2>> elementaryK(nBars);
    // for (Index b{0}; b < nBars; ++b) {
    //   elementaryK[b] =
    //       tensorProduct(connectivityMatrix, elementaryApplicationK[b]);
    // }
    // std::cout << elementaryK;

    // Connectivity Matrix
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
    std::cout << std::boolalpha << (det(assemblyStiffnessK) == 0) << std::endl;

    // Apply Dirichlet BC via elimination (fem1D style) instead of penalization
    assert(nodeImposed.size() == displacementImposed.size() &&
           "Size of imposed nodes must be consistent!");

    // Build flattened RHS from external forces (no penalty modifications)
    Vector<double> forceF = flatten(externalForce, nNodes, d);

    // Save original stiffness and RHS for reaction calculation
    Matrix<double, d * nNodes, d * nNodes> K_original = assemblyStiffnessK;
    Matrix<double, d * nNodes, 1> F_original_m{};
    for (Index i = 0; i < d * nNodes; ++i)
      F_original_m(i, 0) = forceF[i];

    // Build list of prescribed DOFs and their values
    std::vector<Index> prescribedDOFs;
    std::vector<double> prescribedVals;
    for (Index idx = 0; idx < nodeImposed.size(); ++idx) {
      Index n = nodeImposed[idx];
      for (Index k = 0; k < d; ++k) {
        prescribedDOFs.push_back(d * n + k);
        prescribedVals.push_back(displacementImposed[idx][k]);
      }
    }

    // Apply elimination to assemblyStiffnessK and forceF
    for (size_t p_i = 0; p_i < prescribedDOFs.size(); ++p_i) {
      Index p = prescribedDOFs[p_i];
      double gval = prescribedVals[p_i];
      for (Index i = 0; i < d * nNodes; ++i) {
        if (i == p)
          continue;
        forceF[i] -= gval * assemblyStiffnessK(i, p);
      }
      for (Index i = 0; i < d * nNodes; ++i) {
        assemblyStiffnessK(i, p) = 0.0;
        assemblyStiffnessK(p, i) = 0.0;
      }
      assemblyStiffnessK(p, p) = 1.0;
      forceF[p] = gval;
    }

    // Solve the linear system to find displacement increment
    Vector<double> deltaU{solveLinearSystem(assemblyStiffnessK, forceF)};

    std::cout << "Increment of displacement vector U: " << deltaU << std::endl;
    totalDispalcement += deltaU;

    // Saving the output
    iteration_array.push_back(iteration);
    deltaU_array.push_back(deltaU);

    // Compute reaction vector from equilibrium: R = K_original * U_full -
    // F_original
    Matrix<double, d * nNodes, 1> U_full_m{};
    for (Index i = 0; i < d * nNodes; ++i)
      U_full_m(i, 0) = totalDispalcement[i];
    Matrix<double, d * nNodes, 1> Rm = K_original * U_full_m - F_original_m;
    Vector<Vector<double>> reactionR(nNodes);
    for (Index ii = 0; ii < nNodes; ++ii) {
      reactionR[ii] = Vector<double>(d);
      for (Index kk = 0; kk < d; ++kk) {
        reactionR[ii][kk] = Rm(ii * d + kk, 0);
      }
    }

    // Update the position
    nodes += unflatten(deltaU, nNodes, d);
    // Update the force

    iteration++;
  }
  std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";
  return 0;
}