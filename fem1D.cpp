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
//  u(0) = g : Dirichlet's Boundary Condition (Must at least 1 condition to find
//  unique solution!)
// u_x(0) = h : Neuman's Boundary Condition

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

Vector<Vector<std::function<double(double)>>> N(nNodes), N_x(nNodes);
Matrix<double, nNodes, nNodes> K{};
Matrix<double, nNodes, 1> U{};
Matrix<double, nNodes, 1> F{};
Vector<double> length(nElements);
Vector<Index> eleOrigin{0, 1};
Vector<Index> eleEnd{1, 2};
constexpr double EA = 1.0;
Vector<double> k(nElements);
auto rhsFunction = [](auto x) { return x; };
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

  // Elements, shape functions on elements
  for (Index e{0}; e < nElements; ++e) {
    int i = eleOrigin[e];
    int j = eleEnd[e];
    // element[e].resize(2);
    // element[e][0] = nodes[i];
    // element[e][1] = nodes[j];
    double x_i = nodes[i];
    double x_j = nodes[j];
    element[e] = {x_i, x_j};
    length[e] = constexpr_fabs(x_j - x_i);
    k[e] = EA / length[e];

    // Local indices within the element (0 or 1)
    int local_i = 0;
    int local_j = 1;
    auto shapefunction_i = [=](auto x) {
      return basisLagrange(local_i, element[e], x);
    };
    auto dShape_i = [=](double x) { return automaticDiff(shapefunction_i, x); };
    auto shapefunction_j = [=](auto x) {
      return basisLagrange(local_j, element[e], x);
    };
    auto dShape_j = [=](double x) { return automaticDiff(shapefunction_j, x); };
    N[i].push_back(dShape_i);
    N[j].push_back(dShape_j);
    N_x[i].push_back(dShape_i);
    N_x[j].push_back(dShape_j);
  }

  // Shape function
  for (Index i{0}; i < nNodes; ++i) {
    for (Index j{i}; j < nNodes; ++j) {
      // General integrand: sum all local derivative contributions for node i
      // and j
      auto integrand = [&, i, j](double x) {
        double si = 0.0;
        for (Index p{0}; p < static_cast<Index>(N_x[i].size()); ++p)
          si += N_x[i][p](x);
        double sj = 0.0;
        for (Index q{0}; q < static_cast<Index>(N_x[j].size()); ++q)
          sj += N_x[j][q](x);
        return si * sj;
      };
      K(i, j) = integrationGauss1D(a, b, integrand, 2);
    }
  }
  //   K.reflect();

  Matrix<double, 3, 3> K_a = {2.0, -2.0, 0.0, -2.0, 4.0, -2.0, 0.0, -2.0, 2.0};
  std::cout << K << std::endl;

  //   //   Handle Dirichlet's condition
  //   if (g == 0) {
  //     //
  //   } else { ///
  //   }

  //   for (Index j{0}; j < nNodes; ++j) {
  //     auto shape_j = [=](auto x) { return basisLagrange(j, nodes, x) * x; };
  //     auto integrand = [&](double x) { return x * shape_j(x); };
  //     F[j] = shape_j(0.0) + integrationGauss1D(a, b, integrand, 2);
  //   }
  //   std::cout << K << std::endl;
  //   std::cout << F << std::endl;

  //   // Apply Dirichlet condition at node 0 (u(0) = g) if specified
  //   if (g != 0) {
  //     Index p = 0; // Dirichlet node index
  //     double u_d = g;
  //     // modify RHS: F_i <- F_i - K(i,p)*u_d for i != p
  //     for (Index i = 0; i < nNodes; ++i) {
  //       if (i == p)
  //         continue;
  //       F[i] -= K(i, p) * u_d;
  //     }
  //     // zero out column and row p
  //     for (Index i = 0; i < nNodes; ++i) {
  //       K(i, p) = 0.0;
  //       K(p, i) = 0.0;
  //     }
  //     K(p, p) = 1.0;
  //     F[p] = u_d;
  //   }
  U = solveLinearSystem(K, F);
  std::cout << U << std::endl;
  return 0;
}
