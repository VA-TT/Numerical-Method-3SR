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
  Vector<double> x_i;
  for (Index i{0}; i < nNodes; ++i) {
    x_i.push_back((b - a) / (nNodes - 1) * i + a);
  }
  for (Index e{0}; e < nElements; ++e) {
    element[e][e] = x_i[e];
    element[e][e + 1] = x_i[e + 1];
  }
  std::cout << element << std::endl;
  return 0;
}
