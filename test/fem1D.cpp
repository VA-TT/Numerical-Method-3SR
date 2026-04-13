#include "library/DualDifferentiation.h"
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

// Exact solution
double solution(double x) { return (-std::pow(x, 3) + 9 * x + 6) / 6; }
// Input parameters
namespace modelParameters {
// Problem's domain (a,b)
double a{0.0};
double b{1.0};
constexpr double g{1.0};
constexpr double h{1.0};
constexpr Index nNodes{6};
constexpr Index nElements{nNodes - 1};
DynamicVector<Index> eleOrigin{0, 1, 2, 3, 4};
DynamicVector<Index> eleEnd{1, 2, 3, 4, 5};
constexpr double EA = 1.0;

auto rhsFunction = [](auto x) { return x; };

DynamicVector<DynamicVector<std::function<double(double)>>> N(nNodes),
    N_x(nNodes);
Matrix<double, nNodes, nNodes> K{};
Matrix<double, nNodes, 1> U{};
Matrix<double, nNodes, 1> F{};

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
void assemble() {}
void applyBC() {}

int main() {
  // Equally divied
  using namespace modelParameters;

  DynamicVector<double> nodes{generateMesh(a, b, nNodes)};
  std::cout << nodes << "'\n";

  // Local indices within the element (0 or 1)
  int local_i = 0;
  int local_j = 1;

  // Elements, shape functions on elements
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

  // Rigidity Matrix K: assemble element-wise contributions
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

  // Force Vector
  std::cout << K << std::endl;

  //   Handle Dirichlet's condition
  F[0] = g;
  for (Index i = 1; i < nNodes; ++i) {
    F[i] -= g * K(i, 0);
  }
  for (Index i = 0; i < nNodes; ++i) {
    K(i, 0) = 0.0;
  }
  for (Index j = 0; j < nNodes; ++j) {
    K(0, j) = 0.0;
  }
  K(0, 0) = 1.0;

  std::cout << F << std::endl;

  // Solve the linear system to obtain deplacement vector U
  U = solveLinearSystem(K, F);
  std::cout << U << std::endl;

  return 0;
}
