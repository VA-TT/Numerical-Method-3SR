#include "library/DualDiffrentiation.h"
#include "library/Matrix.h" //Approximative Comparsion
#include "library/Vector.h"
#include "library/clock.h"
#include "library/gaussQuadrature.h"
#include "library/interpolate.h"
#include "library/physicConstants.h"
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
constexpr double a{0.0};
constexpr double b{1.0};
constexpr double L{b - a};

// Input
constexpr double E{4 * constants::pi * constants::pi};
constexpr double rho{1.0}; // gonna change
constexpr double v0{0.1};
constexpr double vp{1};
constexpr double xp{0.5};
constexpr double c{constexpr_sqrt(E / rho)};
constexpr double duration{10.0};
constexpr double dt{0.01};
constexpr double nSteps{duration / dt};
constexpr double x_loc{1.0};
constexpr double w{1.0 / L * c};

constexpr Index nNodes{2};             // Numbers of nodes
constexpr Index nElements{nNodes - 1}; // Numbers of elements
constexpr Index nMPs{1};               // Numbers material points

////////////////////////////////////////////////////////////
//////////////////////////// END ///////////////////////////
////////////////////////////////////////////////////////////

// Initiate needed containers
Vector<Vector<std::function<double(double)>>> N(nNodes), N_x(nNodes);

Vector<double> nodes{};
} // namespace modelParameters

// Analytic solution for uniform load: u(x) = (1/EA)*( f*(L*x - x^2/2) + F*x )
double velocity_solution(double t) {
  using namespace modelParameters;
  return v0 * std::cos(w * t);
}

double position_solution(double t) {
  using namespace modelParameters;
  return x_loc * std::exp(v0 / (L * w) * std::sin(w * t));
}

// Mesh generated function
Vector<double> generateMesh(double a, double b, Index n) {
  Vector<double> nodes;
  for (Index i{0}; i < n; ++i) {
    nodes.push_back((b - a) / (n - 1) * i + a);
  }
  return nodes;
}

int main() {
  // Equally divied
  Timer t;
  using namespace modelParameters;
  nodes = generateMesh(a, b, nNodes);
  std::cout << nodes;
  std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";
  return 0;
}
