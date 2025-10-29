#include "library/DualDiffrentiation.h"
#include "library/Matrix.h" //Approximative Comparsion
#include "library/Vector.h"
#include "library/clock.h" // time-measuring
#include <fstream>         //working with files
#include <iomanip>         //tab
#include <iostream>
#include <numbers>     // for std::numbers::pi
#include <stdexcept>   //throw exception
#include <type_traits> // precision
#include <vector>

// Model
#if 0
|o-----o
      / |\
 ^i2 /  |  \> F
 |  /    theta
|o/---> i1
#endif

// Model Parameters
namespace modelParameters {
constexpr Index n{2}; // number of equilibrium equations

// Unit vectors
Vector<double> i1{1.0, 0.0};
Vector<double> i2{0.0, 1.0};

// Bar's parameters
double E{25e9}; // Module Young
double a{10.0}; // Distance between node 0 and 1 = bar length

Vector<double> x0{(i1 * a + i2 * a)};
double l01{magnitude(x0)};          // length bar 1
double l02{magnitude(x0 - a * i2)}; // length bar 2
double force{150e3};                // Imposed Force
double theta{20}; // Inclined angle of the force with relative to vertical
double thetaRadian{theta * std::numbers::pi / 180}; // converte to radian

Vector<double> f1{force * std::sin(thetaRadian) *
                  i1}; // force's horizontal component
Vector<double> f2{-force * std::cos(thetaRadian) *
                  i2}; // force's vertical component
Vector<double> externalForce{f1 + f2};

// Section dimension
double b1{0.02}, h1{0.05};
double b2{0.02}, h2{0.05};
double A1{b1 * h1}, A2{b2 * h2};
double alpha1{E * A1}, alpha2{E * A2};

// tolerance for Newton's method
double epsilon{1e-8};
int max_iteration{100};
} // namespace modelParameters

// Constitutive law
template <typename T>
inline auto constitutiveLaw(double alpha, T l, double l0) {
  return (alpha * (l - l0) / l0);
  // return (alpha * (l * l - l0 * l0) / (2 * l0 * l0));
  // return (alpha * std::log(l / l0));
}

int main() {
  Timer t;
  using namespace modelParameters;

  // Define func as a lambda that matches the expected signature for
  // automaticDiff
  auto func1 = [](auto x1) { return constitutiveLaw(alpha1, x1, l01); };
  auto func2 = [](auto x2) { return constitutiveLaw(alpha2, x2, l02); };

  // Setting up
  int iteration{0};

  std::vector<int> iteration_array;
  std::vector<Vector<double>> deltaX_array;
  iteration_array.reserve(max_iteration);
  deltaX_array.reserve(max_iteration);
  Vector<double> deltaX(2);
  Vector<double> Fk(2);

  std::cout << "=== Newton-Raphson Iteration ===" << std::endl;

  Vector<double> x(2);
  double l1{}, l2{};
  Vector<double> e1(2), e2(2);
  Matrix<double, n, n> nablaF{};
  auto I = Matrix<double, 2, 2>::identity();
  Vector<double> deltaX_increment(2);

  // Use Newton_Raphson method to find the approximation with tolerance
  while (iteration < max_iteration) {
    // Calculate updated informations
    x = x0 + deltaX;
    l1 = magnitude(x);
    l2 = magnitude(x - a * i2);
    e1 = x / l1;
    e2 = (x - a * i2) / l2;

    // Calculate F at iteration k
    Fk = -func1(l1) * e1 - func2(l2) * e2 + externalForce;
    std::cout << "Iteration " << iteration << ": |Fk| = " << magnitude(Fk)
              << std::endl;

    // Use F as the quality of approximation
    if (magnitude(Fk) < epsilon) {
      std::cout << "Solution founded.\n";
      break;
    }

    // Calculate NablaF at itaration k
    nablaF = Matrix<double, n, n>{};
    nablaF += -automaticDiff(func1, l1) * tensorProduct<n, n>(e1, e1);
    nablaF += -automaticDiff(func2, l2) * tensorProduct<n, n>(e2, e2);
    nablaF += -(func1(l1) / l1) * (I - tensorProduct<n, n>(e1, e1));
    nablaF += -(func2(l2) / l2) * (I - tensorProduct<n, n>(e2, e2));

    // Solve the linear system
    std::cout << tensorProduct<n, n>(e1, e1) << '\n';
    std::cout << std::setprecision(12);
    std::cout << "[DEBUG] l1=" << l1 << " l01=" << l01
              << " func1(l1)=" << func1(l1)
              << " expected dfunc1/dl1=" << (alpha1 / l01)
              << " autoDiff(func1,l1)=" << automaticDiff(func1, l1) << '\n';
    std::cout << "nabla F: " << nablaF << '\n';
    std::cout << "Fk: " << Fk << '\n';
    deltaX_increment = solveLinearSystem(nablaF, -Fk);
    deltaX += deltaX_increment;

    // Saving the output
    iteration_array.push_back(iteration);
    deltaX_array.push_back(deltaX);

    iteration++;
  }

  if (iteration >= max_iteration) {
    std::cout << "Failed to converge after " << max_iteration << " iterations"
              << std::endl;
  }

  std::cout << "\nFinal displacement: " << deltaX << std::endl;
  std::cout << "Final position: " << x0 + deltaX << "\n";
  std::cout << "Total iterations: " << iteration << std::endl;
  std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";

  // Luu iteration va gia tri delta x de plot
  return 0;
}
