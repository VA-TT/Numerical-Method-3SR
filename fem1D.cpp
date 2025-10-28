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

// U_xx + x = 0, 0 < x < 1
// u(0) = g : Dirichlet's Boundary Condition
// u_x(0) = 1 : Neuman's Boundary Condition

// Exact solution
double solution(double x) { return (std::pow(x, 3) + 9 * x + 6) / 6; }
namespace modelParameters {
constexpr double g{1.0};
constexpr Index nNodes{3};    // Divied elements
constexpr Index nElements{2}; // Divied elements

} // namespace modelParameters

int main() { return 0; }
