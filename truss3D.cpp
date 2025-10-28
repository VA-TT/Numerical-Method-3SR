#include "library/Matrix.h" //Approximative Comparsion
#include "library/Vector.h"
#include "library/clock.h"
#include "library/dualDiffrentiation.h"
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
constexpr Index nNodes{6}; // number of nodes
constexpr Index nBars{8};  // number of bars

Vector<Vector<double>> nodes{
    {0.0, 0.0, 0.0},   // Node 0
    {10.0, 0.0, 0.0},  // Node 1
    {0.0, 0.0, 10.0},  // Node 2
    {10.0, 0.0, 10.0}, // Node 3
    {0.0, 0.0, 20.0},  // Node 4
    {10.0, 0.0, 20.0}  // Node 5
};
// Bar connectivity: store node indices for each bar's origin and end
Vector<Index> barOrigin{0, 2, 1, 2, 2, 3, 3, 4};
Vector<Index> barEnd{2, 1, 3, 3, 4, 4, 5, 5};
Vector<Vector<double>> vectorBars(nBars), unitVectorBars(nBars);
Vector<double> lengthBars(nBars);
Vector<Index> nodeImposed{0, 1};
Vector<Vector<double>> displacementImposed{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
Vector<Index> nodeFree(nNodes - nodeImposed.size());

// Section dimension
double youngModulus{25e9}; // Module Young
double b{0.02}, h{0.05};
double A{b * h};
double alpha{youngModulus * A};

//  External loads at nodes - Direct 3D force vectors
Vector<Vector<double>> externalForce{
    {0.0, 0.0, 0.0},   // Node 0: no force
    {0.0, 0.0, 0.0},   // Node 1: no force
    {0.0, 0.0, -15.0}, // Node 2: 15N downward (z-direction)
    {0.0, 0.0, -15.0}, // Node 3: 15N downward
    {0.0, 0.0, -15.0}, // Node 4: 15N downward
    {0.0, 0.0, -15.0}  // Node 5: 15N downward
};

} // namespace modelParameters

int main() {
  Timer t;
  using namespace modelParameters;
  // Checking the input
  assert(nodes.size() == nNodes && "numbers of nodes must be consistent!");
  assert(barOrigin.size() == nBars && barEnd.size() == nBars &&
         "numbers of bars must be consistent!");
  for (Index b{0}; b < nBars; ++b) {
    // compute bar vector from node coordinates
    vectorBars[b] = nodes[barEnd[b]] - nodes[barOrigin[b]];
    lengthBars[b] = magnitude(vectorBars[b]);
    unitVectorBars[b] = vectorBars[b] / lengthBars[b];
  }
  std::cout << "Unit vectors for each bar:\n" << unitVectorBars << std::endl;

  // Verify external forces
  assert(externalForce.size() == nNodes &&
         "Number of external forces must match number of nodes");
  for (Index i = 0; i < nNodes; ++i) {
    assert(externalForce[i].size() == d &&
           "Each force vector must have d components");
  }
  std::cout << "External forces at nodes:\n" << externalForce << std::endl;

  // Setting up
  int iteration{0};
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

  // Rigidity in small deformation configuration: N = k e (u2 - u1)
  Vector<double> k(nBars);
  for (Index b{0}; b < nBars; ++b) {
    k[b] = youngModulus * A / lengthBars[b];
  }

  std::cout << k;
  Vector<Matrix<double, d, d>> elementaryApplicationK(nBars);
  for (Index b{0}; b < nBars; ++b) {
    elementaryApplicationK[b] =
        k[b] * tensorProduct<d, d>(unitVectorBars[b], unitVectorBars[b]);
  }
  std::cout << elementaryApplicationK;

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
    // std::cout << C[n].getCols() << "  and  " << C[n].getRows() << std::endl;
  }

  Matrix<double, d * nNodes, d * nNodes> assemblyStiffnessK{};

  for (Index b{0}; b < nBars; ++b) {
    assemblyStiffnessK += (C[barEnd[b]] - C[barOrigin[b]]).transpose() *
                          elementaryApplicationK[b] *
                          (C[barEnd[b]] - C[barOrigin[b]]);
  }
  // Check the singularity of stiffness matrix K
  std::cout << std::boolalpha << (det(assemblyStiffnessK) == 0) << std::endl;

  // Penalization method (encastree)
  // Stifness matrix and force adjustment: K' and F'
  assert(nodeImposed.size() == displacementImposed.size() &&
         "Size of imposed nodes must be consistent!");
  double penalty{1.0 / 1e-8}; // 1/epsilon
  std::cout << penalty << std::endl;
  int i{};
  for (auto n : nodeImposed) {
    for (int k = 0; k < d; ++k) // Encastre tous 3 directions
    {
      i = d * n + k;
      assemblyStiffnessK(i, i) += penalty;
      externalForce[d * n][k] = penalty * displacementImposed[n][k];
    }
  }

  std::cout << assemblyStiffnessK << std::endl;

  std::cout << std::boolalpha << (det(assemblyStiffnessK) == 0) << std::endl;

  // Flatting vector F
  Vector<double> forceF(nNodes * d);
  for (Index i{0}; i < forceF.size(); ++i) {
    forceF[i] = externalForce[i / d][i % d];
  }

  // Solve the linear system K'U' = F':
  Vector<double> increment_displacement{
      solveLinearSystem(assemblyStiffnessK, forceF)};

  // // Update the position
  // Vector<double> totalDispalcementU;
  // totalDispalcementU += increment_displacement;

  // // Saving the output
  // iteration_array.push_back(iteration);
  // deltaX_array.push_back(deltaX);

  // Calculate the reaction at the end of the loop
  // Index i{0};
  // Vector<double> reactionR{nodeImposed};
  // for (auto n : nodeImposed) {
  //   for (int k = 0; k < 3; k++) // Encastre tous 3 directions
  //   {
  //     i = 3 * n + k;
  //     reactionR[n][k] =
  //         penalty * (totalDispalcementU[i] - nodeImposed[n][k]);
  //   }
  // }

  // iteration++;
  std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";
  return 0;
}