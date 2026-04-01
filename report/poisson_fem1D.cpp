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
//  U_xx + 6x = 0, 0 < x < 1 (a < x < b)
//  u(0) = g1 = 1 : Dirichlet's Boundary Condition (Must at least 1 condition to
//  u(1) = g2 = 0 : Dirichlet's Boundary Condition (Must at least 1 condition to
//  find unique solution!)
//  u_x(1) = h = 1: Neuman's Boundary Condition

// Exact solution
double solution(double x) { return (-std::pow(x, 3) + x); }
// Input parameters
namespace modelParameters {
// Problem's domain (a,b)
double a{0.0};
double b{1.0};
constexpr double g1{0.0};
constexpr double g2{0.0};
constexpr double h{0.0};
constexpr Index nNodes{6};
constexpr Index nElements{nNodes - 1};
DynamicVector<Index> eleOrigin{0, 1, 2, 3, 4};
DynamicVector<Index> eleEnd{1, 2, 3, 4, 5};
constexpr double EA = 1.0;

auto rhsFunction = [](auto x) { return x; };

DynamicVector<DynamicVector<std::function<double(double)>>> N(nNodes), N_x(nNodes);
Matrix<double, nNodes, nNodes> K{};
Matrix<double, nNodes, 1> U{};
Matrix<double, nNodes, 1> F{};
DynamicVector<double> nodes{};
DynamicVector<DynamicVector<double>> element(nElements);
DynamicVector<double> k(nElements);
DynamicVector<double> length(nElements);

} // namespace modelParameters

DynamicVector<double> generateMesh(double a, double b, Index n) {
  DynamicVector<double> nodes;
  for (Index i{0}; i < n; ++i) {
    nodes.push_back((b - a) / (n - 1) * i + a);
  }
  return nodes;
}

// Cluster code
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
    for (Index k{0}; k < nNodes; ++k) {
      if (k == i) {
        N[k].push_back(shapefunction_i);
        N_x[k].push_back(dShape_i);
      } else if (k == j) {
        N[k].push_back(shapefunction_j);
        N_x[k].push_back(dShape_j);
      } else {
        N[k].push_back([](double x) { return 0.0; });
        N_x[k].push_back([](double x) { return 0.0; });
      }
    }
  }
}
void assemble() { // Rigidity Matrix K: assemble element-wise contributions
  using namespace modelParameters;
  int local_i = 0;
  int local_j = 1;
  for (Index i = 0; i < nNodes; ++i) {
    for (Index j = i; j < nNodes; ++j) {
      // sum contributions from each element
      for (Index e = 0; e < nElements; ++e) {
        // for element e, use the e-th local derivative function stored in N_x
        auto integrand_K = [&, i, j, e](double x) {
          return N_x[i][e](x) * N_x[j][e](x);
        };
        K(i, j) += integrationGauss1D(element[e][local_i], element[e][local_j],
                                      integrand_K, 4);
      }
    }
  }
  K.reflect(); // K is symmetric
  assert(det(K) == 0 &&
         "Before applying the boundary condition, K must be singular!");

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
  // Neuman's condition; N3(x=1) = 1 (activate at the last node)
  F[nNodes - 1] += h * 1.0;
}
void applyBC(Index node, double g) { //   Handle Dirichlet's condition
  using namespace modelParameters;
  F[node] = g;
  for (Index i = 1; i < nNodes; ++i) {
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

int main() {
  // Equally divied
  using namespace modelParameters;

  nodes = generateMesh(a, b, nNodes);
  shapeFunction();
  assemble();
  applyBC(g1, 0);
  applyBC(g1, 6);
  std::cout << nodes << "'\n";
  std::cout << K << std::endl;
  std::cout << F << std::endl;
  U = solveLinearSystem(K, F);
  std::cout << U << std::endl;
  // // ------------------ Error check against exact solution ------------------
  // // Compute per-node absolute error, max error and RMS error
  // double max_err = 0.0;
  // double sum_sq = 0.0;
  // for (Index i = 0; i < nNodes; ++i) {
  //   double x = nodes[i];
  //   double u_num = U[i];
  //   double u_ex = solution(x);
  //   double err = std::fabs(u_num - u_ex);
  //   if (err > max_err) max_err = err;
  //   sum_sq += err * err;
  //   std::cout << "x=" << x << "  U_num=" << u_num << "  U_ex=" << u_ex
  //             << "  err=" << err << "\n";
  // }
  // double rms = std::sqrt(sum_sq / static_cast<double>(nNodes));
  // std::cout << "max error = " << max_err << "  RMS error = " << rms <<
  // std::endl;

  return 0;
}
