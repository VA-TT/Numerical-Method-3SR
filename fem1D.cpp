#include "library/DualDiffrentiation.h"
#include "library/Matrix.h" //Approximative Comparsion
#include "library/Vector.h"
#include "library/clock.h"
#include "library/gaussQuadrature.h"
#include "library/interpolate.h"
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
//  find unique solution!)
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
// Dirichlet and Neumann boundary values (defaults)
constexpr double g{0.0};  // Dirichlet at x=0 (u(0)=g)
constexpr double h{10.0}; // Neumann at x=1 (N(1)=h)
// Uniform distributed load (default). If you want a spatially varying
// load, replace rhsFunction accordingly (e.g. [](auto x){ return x; }).
constexpr double uniform_load{5.0};
constexpr Index nNodes{6};             // Numbers of nodes
constexpr Index nElements{nNodes - 1}; // Numbers of elements
// Start and end of elements, will be generated automatically in main
Vector<Index> eleOrigin;
Vector<Index> eleEnd;
constexpr double EA = 10.0;

// Right hand side function: default is uniform load `uniform_load`.
auto rhsFunction = [](auto x) { return uniform_load; };

////////////////////////////////////////////////////////////
//////////////////////////// END ///////////////////////////
////////////////////////////////////////////////////////////

// Initiate needed containers
Vector<Vector<std::function<double(double)>>> N(nNodes), N_x(nNodes);
Matrix<double, nNodes, nNodes> K{};
Matrix<double, nNodes, 1> U{};
Matrix<double, nNodes, 1> F{};
Matrix<double, nNodes, 1> R{};
Vector<double> nodes{};
Vector<Vector<double>> element(nElements);
Vector<double> k(nElements);
Vector<double> length(nElements);
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
Vector<double> generateMesh(double a, double b, Index n) {
  Vector<double> nodes;
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

  nodes = generateMesh(a, b, nNodes);
  std::cout << nodes << "'\n";
  // Generate element connectivity automatically from number of nodes
  eleOrigin = Vector<Index>();
  eleEnd = Vector<Index>();
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

  std::cout << K << std::endl;
  std::cout << F << std::endl;
  U = solveLinearSystem(K, F);
  std::cout << U << std::endl;

  // Calculate reaction forces using helper function
  calculateReactions(K_original, F_original, U);
  // Compute internal axial force at x=0 (from first element) and print
  double N0 = 0.0;
  if (length.size() > 0) {
    // axial force N = EA * du/dx approximated by linear element slope
    N0 = EA * (U[1] - U[0]) / length[0];
  }
  std::cout << "\nReaction force at node " << first_node
            << " (x=" << nodes[first_node] << "): R = " << R[first_node]
            << std::endl;
  std::cout << "Internal axial N(0) = " << N0 << "   R[0] = " << R[0]
            << "   N(0)+R[0] = " << (N0 + R[0]) << std::endl;
  std::cout << "All reaction forces: " << R << std::endl;
  // Export numeric nodal data for plotting (Asymptote)
  {
    std::ofstream dataFile("report/fem1D_data.asy");
    dataFile << "// Generated nodal data from fem1D_test.cpp\n";
    dataFile << "real[] xn = {";
    for (Index i = 0; i < modelParameters::nNodes; ++i) {
      dataFile << nodes[i];
      if (i + 1 < modelParameters::nNodes)
        dataFile << ",";
    }
    dataFile << "};\n";
    dataFile << "real[] un = {";
    for (Index i = 0; i < modelParameters::nNodes; ++i) {
      dataFile << U[i];
      if (i + 1 < modelParameters::nNodes)
        dataFile << ",";
    }
    dataFile << "};\n";
    dataFile.close();
  }
  // Also export a simple two-column text file for Asymptote input
  {
    // Compute element axial forces (constant per linear element)
    Vector<double> N_elem;
    for (Index e = 0; e < modelParameters::nElements; ++e) {
      Index i = modelParameters::eleOrigin[e];
      Index j = modelParameters::eleEnd[e];
      double Ne = modelParameters::EA * (U[j] - U[i]) / length[e];
      N_elem.push_back(Ne);
    }

    // Compute nodal axial forces by averaging adjacent element forces
    Vector<double> N_node;
    N_node.resize(modelParameters::nNodes);
    for (Index idx = 0; idx < modelParameters::nNodes; ++idx) {
      if (idx == 0) {
        N_node[idx] = (N_elem.size() > 0) ? N_elem[0] : 0.0;
      } else if (idx == modelParameters::nNodes - 1) {
        N_node[idx] = (N_elem.size() > 0) ? N_elem[N_elem.size() - 1] : 0.0;
      } else {
        double left = N_elem[idx - 1];
        double right = N_elem[idx];
        N_node[idx] = 0.5 * (left + right);
      }
    }

    std::ofstream txtFile("report/fem1D_6nodes.txt");
    // Write three columns: x, u(x), N(x)
    for (Index i = 0; i < modelParameters::nNodes; ++i) {
      txtFile << nodes[i] << " " << U[i] << " " << N_node[i] << "\n";
    }
    txtFile.close();
  }
  // ------------------ Error check against exact solution ------------------
  // Compute per-node absolute error, max error and RMS error
  double max_err = 0.0;
  double sum_sq = 0.0;
  for (Index i = 0; i < nNodes; ++i) {
    double x = nodes[i];
    double u_num = U[i];
    double u_ex = solution(x);
    double err = std::fabs(u_num - u_ex);
    if (err > max_err)
      max_err = err;
    sum_sq += err * err;
    std::cout << "x=" << x << "  U_num=" << u_num << "  U_ex=" << u_ex
              << "  err=" << err << "\n";
  }
  double rms = std::sqrt(sum_sq / static_cast<double>(nNodes));
  std::cout << "max error = " << max_err << "  RMS error = " << rms
            << std::endl;
  std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";
  return 0;
}
