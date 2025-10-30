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
//  u(0) = g : Dirichlet's Boundary Condition
//  u_x(0) = h : Neuman's Boundary Condition

// Exact solution
double solution(double x) { return (std::pow(x, 3) + 9 * x + 6) / 6; }
namespace modelParameters {
// Problem's domain (a,b)
double a{0.0};
double b{1.0};
constexpr double g{1.0};
constexpr double h{1.0};
constexpr Index nNodes{3};
constexpr Index nElements{nNodes - 1};

Vector<double> N(nNodes), N_x(nNodes);
Matrix<double, nNodes, nNodes> K{};
Matrix<double, nNodes, 1> U{};
Matrix<double, nNodes, 1> F{};
Vector<double> length(nElements);
Vector<Index> eleOrigin{0, 1};
Vector<Index> eleEnd{1, 2};
constexpr double EA = 1.0;

Vector<Vector<double>> element(nElements);
} // namespace modelParameters

Vector<double> generateMesh(double a, double b, Index n) {
  Vector<double> nodes;
  for (Index i{0}; i < n; ++i) {
    nodes.push_back((b - a) / (n - 1) * i + a);
  }
  return nodes;
}

void assembleK() {}
void applyBC() {}
double rhsFuction(double x) { return x; }

int main() {
  // Equally divied
  using namespace modelParameters;

  Vector<double> nodes{generateMesh(a, b, nNodes)};
  std::cout << nodes << "'\n";

  // Elements
  for (Index e{0}; e < nElements; ++e) {
    element[e].resize(2);
    double x1 = nodes[eleOrigin[e]];
    double x2 = nodes[eleEnd[e]];
    length[e] = EA * constexpr_fabs(x2 - x1);
    // Element stiffness matrix adding to the global stifness matrix
    K(e, e) += 1 / length[e];
    K(e, e + 1) -= 1 / length[e];
    K(e + 1, e) -= 1 / length[e];
    K(e + 1, e + 1) += 1 / length[e];
    // Element force vector adding to the global force vector
    auto integrand1 = [=](auto x) { return (x2 - x) / h * rhsFuction(x); };
    auto integrand2 = [=](auto x) { return (x - x1) / h * rhsFuction(x); };
    F[e] += integrationGauss1D(x1, x2, integrand1);
    F[e + 1] += integrationGauss1D(x1, x2, integrand2);
  }

  //   // Shape function
  //   for (Index i{0}; i < nNodes; ++i)
  //     for (Index j{i}; j < nNodes; ++j) {
  //       // Lambda parametes must be auto in order to take the derivative with
  //       // Dual
  //       auto shape_i = [=](auto x) { return basisLagrange(i, nodes, x); };
  //       auto shape_j = [=](auto x) { return basisLagrange(j, nodes, x); };
  //       auto dShape_i = [=](double x) { return automaticDiff(shape_i, x); };
  //       auto dShape_j = [=](double x) { return automaticDiff(shape_j, x); };
  //       auto integrand = [=](double x) { return dShape_i(x) * dShape_j(x); };
  //       K(i, j) = integrationGauss1D(a, b, integrand, 2);
  //     }
  //   K.reflect();

  //   Matrix<double, 3, 3> K_a = {2.0, -2.0, 0.0, -2.0, 4.0, -2.0, 0.0,
  //   -2.0, 2.0}; std::cout << K << std::endl;

  // Handle Dirichlet's condition
  //   if (g == 0) {
  //     //
  //   } else { ///
  //   }

  //   for (Index j{0}; j < nNodes; ++j) {
  //     auto shape_j = [=](auto x) { return basisLagrange(j, nodes, x) * x; };
  //     auto integrand = [&](double x) { return x * shape_j(x); };
  //     F[j] = shape_j(0.0) + integrationGauss1D(a, b, integrand, 2);
  //   }
  std::cout << K << std::endl;
  std::cout << F << std::endl;

  // Apply Dirichlet condition at node 0 (u(0) = g) if specified
  if (g != 0) {
    Index p = 0; // Dirichlet node index
    double u_d = g;
    // modify RHS: F_i <- F_i - K(i,p)*u_d for i != p
    for (Index i = 0; i < nNodes; ++i) {
      if (i == p)
        continue;
      F[i] -= K(i, p) * u_d;
    }
    // zero out column and row p
    for (Index i = 0; i < nNodes; ++i) {
      K(i, p) = 0.0;
      K(p, i) = 0.0;
    }
    K(p, p) = 1.0;
    F[p] = u_d;
  }
  U = solveLinearSystem(K, F);
  std::cout << U << std::endl;
  return 0;
}
