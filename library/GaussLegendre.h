#ifndef INTEGRATION_GAUSS_LEGENDRE
#define INTEGRATION_GAUSS_LEGENDRE

#include <array>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>

// Note: In this method, we should use (n) Gauss points to approximate a
// polynomials with degree (2n - 1),
double myfunc(double x) { return std::exp(x); }

std::array<int, 4> pointsNumber = {1, 2, 3, 4};

// Gauss Points
std::array<double, 1> xiOne = {0};
std::array<double, 2> xiTwo = {-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
std::array<double, 3> xiThree = {-std::sqrt(3.0 / 5.0), 0.0,
                                 std::sqrt(3.0 / 5.0)};
inline std::array<double, 4> xiFour = {
    -std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
    -std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
    std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
    std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0))};

// Weights
std::array<double, 1> wOne = {2.0};
std::array<double, 2> wTwo = {1.0, 1.0};
std::array<double, 3> wThree = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
inline std::array<double, 4> wFour = {
    (18.0 - std::sqrt(30.0)) / 36.0, (18.0 + std::sqrt(30.0)) / 36.0,
    (18.0 + std::sqrt(30.0)) / 36.0, (18.0 - std::sqrt(30.0)) / 36.0};

// Implemented in range [a,b] for 2 points
double integrationGaussTwoPoints(double a, double b,
                                 std::function<double(double)> f) {
  double I{};
  double c1{(b - a) / 2.0}; // mapping from [-1,1] to [a,b]
  double c2{(b + a) / 2.0};
  for (int i = 0; i < 2; ++i) {
    I += c1 * wTwo[i] * f(c1 * xiTwo[i] + c2);
  }
  return I;
}

// Implemented in range [a,b] for 4 points
double integrationGaussFourPoints(double a, double b,
                                  std::function<double(double)> f) {
  double I{};
  double c1{(b - a) / 2.0}; // mapping from [-1,1] to [a,b]
  double c2{(b + a) / 2.0};
  for (int i = 0; i < 4; ++i) {
    I += c1 * wFour[i] * f(c1 * xiFour[i] + c2);
  }
  return I;
}

#endif

/// Implement ham chung integrationGauss(n = 4, a, b, f) voi n la so GaussPoint
// Implement 2D

// double integrationGauss2D_two_by_two(double ax, double bx, double ay, double
// by,
//                                      std::function<double(double, double)> f)
//                                      {
//   double c1x = (bx - ax) / 2.0, c2x = (bx + ax) / 2.0;
//   double c1y = (by - ay) / 2.0, c2y = (by + ay) / 2.0;
//   double I = 0.0;
//   for (int i = 0; i < 2; ++i)
//     for (int j = 0; j < 2; ++j) {
//       double x = c1x * xiTwo[i] + c2x;
//       double y = c1y * xiTwo[j] + c2y;
//       I += c1x * c1y * wTwo[i] * wTwo[j] * f(x, y);
//     }
//   return I;
// }

// int main()
// {

//     double expQuadrature{integrationGaussFourPoints(-2, 3, myfunc)};
//     double expBi{integrationGaussTwoPoints(-2, 3, myfunc)};

//     // Print with 15 decimal precisions
//     std::cout
//         << std::fixed << std::setprecision(15);
//     std::cout << "Integration [-1,1] of function e^x: I = " << expBi << ' '
//     << expQuadrature << std::endl;

//     return 0;
// }