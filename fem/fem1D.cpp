#include "../library/DualDifferentiation.h"
#include "../library/Matrix.h" //Approximative Comparsion
#include "../library/Vector.h"
#include "../library/clock.h"
#include "../library/gaussQuadrature.h"
#include "../library/interpolate.h"
#include "vector"
#include <cassert> // for assert
#include <fstream> //working with files
#include <iomanip> //tab
#include <iostream>
#include <numbers>     // for std::numbers::pi
#include <stdexcept>   //throw exception
#include <type_traits> // precision

// 1D problem
//  U_xx + x = 0, 0 < x < 1 (a < x < b)
//  u(0) = g = 1 : Dirichlet's Boundary Condition (Must at least 1 condition to
//  prevent rigid movement and ensure unique solution!)
//  u_x(1) = h = 1: Neuman's Boundary Condition

// Exact solution for u'' + x = 0 with u(0)=1 and u'(1)=1
// Solve: u'' = -x -> u' = -x^2/2 + C1. Using u'(1)=1 -> C1 = 1 + 1/2 = 1.5
// u = -x^3/6 + 1.5*x + C2. Using u(0)=1 -> C2 = 1
// prototype for analytic solution (defined after modelParameters so it can
// use those values)
double solution(double x);

namespace modelParameters {
////////////////////////////////////////////////////////////
///////////////////// Input parameters /////////////////////
////////////////////////////////////////////////////////////
// Problem's domain (a,b)
double a{0.0};
double b{1.0};
// Dirichlet and Neumann boundary values
constexpr double g{0.0};  // Dirichlet at x=0 (u(0)=g)
constexpr double h{10.0}; // Neumann at x=1 (N(1)=h)
// Uniform distributed load . If spatially varying load is wanted, replace
// rhsFunction accordingly (e.g. [](auto x){ return x; }).
constexpr double uniform_load{5.0};
constexpr Index nNodes{6};             // Numbers of nodes
constexpr Index nElements{nNodes - 1}; // Numbers of elements
// Start and end of elements, will be generated automatically in main
DynamicVector<Index> eleOrigin;
DynamicVector<Index> eleEnd;
constexpr double EA = 10.0;

// Right hand side function: default is uniform load `uniform_load`.
auto rhsFunction = [](auto x) { return uniform_load; };

////////////////////////////////////////////////////////////
//////////////////////////// END ///////////////////////////
////////////////////////////////////////////////////////////

// Initiate needed containers
DynamicVector<DynamicVector<std::function<double(double)>>> N(nNodes),
    N_x(nNodes);
Matrix<double, nNodes, nNodes> K{};
Matrix<double, nNodes, 1> U{};
Matrix<double, nNodes, 1> F{};
Matrix<double, nNodes, 1> R{};
DynamicVector<double> nodes{};
DynamicVector<DynamicVector<double>> element(nElements);
DynamicVector<double> k(nElements);
DynamicVector<double> length(nElements);
Index first_node{0};         // index of the last node
Index last_node{nNodes - 1}; // index of the last node

} // namespace modelParameters

// Analytic solution for uniform load: u(x) = (1/EA)*( f*(L*x - x^2/2) + F*x )
double solution(double x) {
  using namespace modelParameters;
  double f = uniform_load;
  double L = b;
  double F = h; // Neumann at x = L
  return (1.0 / EA) * (f * (L * x - 0.5 * x * x) + F * x);
}

// Mesh generated function
DynamicVector<double> generateMesh(double a, double b, Index n) {
  DynamicVector<double> nodes;
  for (Index i{0}; i < n; ++i) {
    nodes.push_back((b - a) / (n - 1) * i + a);
  }
  return nodes;
}

// calculating shape functions and its derivative
// (obviously redundance but objectif here is to show the algorithm
// in non-shorcut form and test out implemented libraries)
void shapeFunction() {
  // Equally divied
  using namespace modelParameters;
  int local_i = 0;
  int local_j = 1;
  for (Index e{0}; e < nElements; ++e) {
    int i = eleOrigin[e];
    int j = eleEnd[e];
    double x_i = nodes[i];
    double x_j = nodes[j];
    element[e] = {x_i, x_j};
    length[e] = constexpr_fabs(x_j - x_i);
    k[e] = EA / length[e];

    // MUST set parameter "x" here to type <<auto>> in order to accept Dual
    // class as input for derivative calculating
    auto shapefunction_i = [=](auto x) {
      return basisLagrange(local_i, element[e], x);
    };
    auto dShape_i = [=](double x) { return automaticDiff(shapefunction_i, x); };
    auto shapefunction_j = [=](auto x) {
      return basisLagrange(local_j, element[e], x);
    };
    auto dShape_j = [=](double x) { return automaticDiff(shapefunction_j, x); };
    for (Index kk{0}; kk < nNodes; ++kk) {
      if (kk == i) {
        N[kk].push_back(shapefunction_i);
        N_x[kk].push_back(dShape_i);
      } else if (kk == j) {
        N[kk].push_back(shapefunction_j);
        N_x[kk].push_back(dShape_j);
      } else {
        N[kk].push_back([](double x) { return 0.0; });
        N_x[kk].push_back([](double x) { return 0.0; });
      }
    }
  }
}

// More likely than assemble, calculate each element in the rigidity matrix K
//  and the force vector F
void assembleKF(Index nodeNeuman, double h) {
  // Rigidity Matrix K: assemble element-wise contributions
  using namespace modelParameters;
  int local_i = 0;
  int local_j = 1;
  for (Index i = 0; i < nNodes; ++i) {
    for (Index j = i; j < nNodes; ++j) {
      // sum contributions from each element
      for (Index e = 0; e < nElements; ++e) {
        // for element e, use the e-th local derivative function stored in N_x
        // include EA in the stiffness integrand so K = \int EA N'_i N'_j
        auto integrand_K = [&, i, j, e](double x) {
          return EA * N_x[i][e](x) * N_x[j][e](x);
        };
        K(i, j) += integrationGauss1D(element[e][local_i], element[e][local_j],
                                      integrand_K, 4);
      }
    }
  }
  K.reflect(); // K is symmetric

  // Force Vector
  F = Matrix<double, nNodes, 1>::zero();
  for (Index e = 0; e < nElements; ++e) {
    double x1 = element[e][local_i];
    double x2 = element[e][local_j];
    for (Index node = 0; node < nNodes; ++node) {
      auto integrand_F = [=](double x) {
        return rhsFunction(x) * N[node][e](x);
      };
      F[node] += integrationGauss1D(x1, x2, integrand_F, 2);
    }
  }
  // Neuman's condition h * N3(x=1) = h * 1 (activate at the last node)
  F[nodeNeuman] += h * 1.0;
}
//   Handle Dirichlet's condition
void applyBC(Index node, double g) {
  using namespace modelParameters;
  // Set prescribed DOF value, then remove its contribution from RHS for all
  // other DOFs: F[i] -= g * K(i,node) for i != node
  F[node] = g;
  for (Index i = 0; i < nNodes; ++i) {
    if (i == node)
      continue;
    F[i] -= g * K(i, node);
  }
  for (Index i = 0; i < nNodes; ++i) {
    K(i, node) = 0.0;
  }
  for (Index j = 0; j < nNodes; ++j) {
    K(node, j) = 0.0;
  }
  K(node, node) = 1.0;
}

// Compute reaction vector R = F_original - K_original * U
void calculateReactions(
    const Matrix<double, modelParameters::nNodes, modelParameters::nNodes>
        &K_original,
    const Matrix<double, modelParameters::nNodes, 1> &F_original,
    const Matrix<double, modelParameters::nNodes, 1> &U) {
  using namespace modelParameters;
  R = -(F_original - K_original * U);
}

int main() {
  // Equally divied
  Timer t;
  using namespace modelParameters;

  std::cout << "=== FEM 1D PROBLEM ===" << std::endl;
  nodes = generateMesh(a, b, nNodes);
  std::cout << "Coordinates of nodes x_i: \n" << nodes << '\n';
  // Generate element connectivity automatically from number of nodes
  eleOrigin = DynamicVector<Index>();
  eleEnd = DynamicVector<Index>();
  for (Index e = 0; e < nElements; ++e) {
    eleOrigin.push_back(e);
    eleEnd.push_back(e + 1);
  }
  shapeFunction();
  assembleKF(last_node, h);

  assert(approximatelyEqualAbsRel(det(K), 0.0));

  // Save original K and F before applying BC (for reaction force calculation)
  Matrix<double, nNodes, nNodes> K_original = K;
  Matrix<double, nNodes, 1> F_original = F;

  applyBC(first_node, g); // Apply u = g1 at first node

  // std::cout << K << std::endl;
  // std::cout << F << std::endl;
  U = solveLinearSystem(K, F);
  std::cout << "Displacement vector U: \n"
            << DynamicVector<double>(U) << std::endl;

  // Calculate reaction forces using helper function
  calculateReactions(K_original, F_original, U);

  std::cout << "Reaction force vector R: \n"
            << DynamicVector<double>(R) << std::endl;

  {
    // Compute element axial forces (constant per linear element)
    DynamicVector<double> N_elem(nElements);
    for (Index e = 0; e < nElements; ++e) {
      Index i = eleOrigin[e];
      Index j = eleEnd[e];
      N_elem[e] = EA * (U[j] - U[i]) / length[e];
    }

    // Assign element axial force to its two end nodes (overwrite if shared)
    DynamicVector<double> N_node(nNodes);
    for (Index e = 0; e < nElements; ++e) {
      Index i = eleOrigin[e];
      Index j = eleEnd[e];
      N_node[i] = N_elem[e];
      N_node[j] = N_elem[e];
    }
    std::cout << "Axial force vector N: \n" << N_node << std::endl;
    std::ofstream txtFile("report/fem1D_6nodes.txt");
    // Write three columns: x, u(x), N(x)
    for (Index i = 0; i < modelParameters::nNodes; ++i) {
      txtFile << nodes[i] << " " << U[i] << " " << N_node[i] << "\n";
    }
    txtFile.close();
  }
  std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";
  return 0;
}
