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
constexpr Index nNodes{3};    // Divied elements
constexpr Index nElements{2}; // Divied elements

Vector<Vector<double>> element(nElements);

} // namespace modelParameters

int main() {
  // Equally divied
  using namespace modelParameters;
  Vector<double> nodes;
  for (Index i{0}; i < nNodes; ++i) {
    nodes.push_back((b - a) / (nNodes - 1) * i + a);
  }
  Vector<Index> eleOrigin{0, 1};
  Vector<Index> eleEnd{1, 2};

  // Shape function

  Vector<double> N(nNodes), N_x(nNodes);

  Matrix<double, nNodes, nNodes> K{};
  Matrix<double, nNodes, 1> U{};
  Matrix<double, nNodes, 1> F{};

  for (Index i{0}; i < nNodes; ++i)
    for (Index j{i}; j < nNodes; ++j) {
      // Lambda parametes must be auto in order to take the derivative with Dual
      auto shape_i = [=](auto x) { return basisLagrange(i, nodes, x); };
      auto shape_j = [=](auto x) { return basisLagrange(j, nodes, x); };
      auto dShape_i = [&](double x) { return automaticDiff(shape_i, x); };
      auto dShape_j = [&](double x) { return automaticDiff(shape_j, x); };
      auto integrand = [&](double x) { return dShape_i(x) * dShape_j(x); };
      K(i, j) = integrationGauss1D(a, b, integrand, 2);
    }
  K.reflect();
  static_assert(
      std::is_same_v<decltype(automaticDiff(
                         std::declval<decltype([](auto) { return 0.0; })>(),
                         0.0)),
                     double>,
      "automaticDiff must return double for double-valued integrand");

  Matrix<double, 3, 3> K_a = {2.0, -2.0, 0.0, -2.0, 4.0, -2.0, 0.0, -2.0, 2.0};
  std::cout << K << std::endl;

  // Handle Dirichlet's condition
  if (g == 0) {
    //
  } else { ///
  }

  for (Index j{0}; j < nNodes; ++j) {
    auto shape_j = [=](auto x) { return basisLagrange(j, nodes, x) * x; };
    auto integrand = [&](double x) { return x * shape_j(x); };
    F[j] = shape_j(0.0) + integrationGauss1D(a, b, integrand, 2);
  }
  std::cout << F << std::endl;

  U = solveLinearSystem(K, F);
  std::cout << U << std::endl;
  return 0;
}
