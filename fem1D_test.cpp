#include "library/DualDiffrentiation.h"
#include "library/Matrix.h" //Approximative Comparsion
#include "library/Vector.h"
#include "library/clock.h"
#include "library/gaussQuadrature.h"
#include "library/interpolate.h"
#include "vector"
#include <cassert> // for assert
#include <fstream> //working with files
#include <functional>
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
    double h = x2 - x1;
    length[e] = EA * std::fabs(h);
    // Element stiffness matrix adding to the global stiffness matrix
    K(e, e) += 1.0 / h;
    K(e, e + 1) -= 1.0 / h;
    K(e + 1, e) -= 1.0 / h;
    K(e + 1, e + 1) += 1.0 / h;
    // Element force vector adding to the global force vector (using hat
    // functions)
    auto integrand1 = [=](double x) { return (x2 - x) / h * rhsFuction(x); };
    auto integrand2 = [=](double x) { return (x - x1) / h * rhsFuction(x); };
    F[e] += integrationGauss1D(x1, x2, integrand1, 2);
    F[e + 1] += integrationGauss1D(x1, x2, integrand2, 2);
  }
  // Optional: reflect symmetric part if needed
  K.reflect();

  // Apply boundary conditions similar to the Julia code
  // Neumann BC on the right end: add u_x0 to last DOF
  double u0 = g;         // Dirichlet left value
  double u_x0 = h;       // Neumann value at right
  F[nNodes - 1] += u_x0; // Neumann contribution

  // Dirichlet BC at first node: enforce u(0) = u0
  for (Index i = 0; i < nNodes; ++i) {
    F[i] -= K(i, 0) * u0;
  }
  for (Index j = 0; j < nNodes; ++j)
    K(0, j) = 0.0;
  for (Index i = 0; i < nNodes; ++i)
    K(i, 0) = 0.0;
  K(0, 0) = 1.0;
  F[0] = u0;

  U = solveLinearSystem(K, F);
  std::cout << "K =\n" << K << std::endl;
  std::cout << "F =\n" << F << std::endl;
  std::cout << "U =\n" << U << std::endl;
  return 0;
}
