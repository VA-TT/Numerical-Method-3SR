#ifndef INTEGRATION_GAUSS_LEGENDRE
#define INTEGRATION_GAUSS_LEGENDRE

#include <array>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>

double myfunc(double x) { return std::exp(x); }

std::array<int, 4> pointsNumber = {1, 2, 3, 4};

std::array<double, 1> xiOne = {0};
std::array<double, 2> xiTwo = {-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
std::array<double, 3> xiThree = {-std::sqrt(3.0 / 5.0), 0.0,
                                 std::sqrt(3.0 / 5.0)};
std::array<double, 4> xiFour = {
    -std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
    -std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
    std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
    std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0))};

std::array<double, 1> wOne = {2.0};
std::array<double, 2> wTwo = {1.0, 1.0};
std::array<double, 3> wThree = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
std::array<double, 4> wFour = {
    (18.0 - std::sqrt(30.0)) / 36.0, (18.0 + std::sqrt(30.0)) / 36.0,
    (18.0 + std::sqrt(30.0)) / 36.0, (18.0 - std::sqrt(30.0)) / 36.0};

// Implemented in range [-1,1] only
double integrationGaussTwoPoints(double a, double b,
                                 std::function<double(double)> f) {
  double I{};
  double c1{(b - a) / 2.0};
  double c2{(b + a) / 2.0};
  for (int i = 0; i < 2; ++i) {
    I += c1 * wTwo[i] * f(c1 * xiTwo[i] + c2);
  }
  return I;
}
double integrationGaussFourPoints(double a, double b,
                                  std::function<double(double)> f) {
  double I{};
  double c1{(b - a) / 2.0};
  double c2{(b + a) / 2.0};
  for (int i = 0; i < 4; ++i) {
    I += c1 * wFour[i] * f(c1 * xiFour[i] + c2);
  }
  return I;
}

#endif
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