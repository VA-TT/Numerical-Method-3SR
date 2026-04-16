#include "library/newtonRaphson.h"
#include <array>
#include <cmath>
#include <iostream>

void printPoint(const StaticVector<double, 3> &p) {
  std::cout << '(' << p[0] << ", " << p[1] << ", " << p[2] << ')';
}

int main() {
  const double sigma = 10.0;
  const double r = 28.0;
  const double beta = 8.0 / 3.0;

  std::cout << "Lorenz fixed points (analytical):\n";
  std::cout << "P0 = (0, 0, 0)\n";
  const double a = std::sqrt(beta * (r - 1.0));
  StaticVector<double, 3> pPlus{a, a, r - 1.0};
  StaticVector<double, 3> pMinus{-a, -a, r - 1.0};
  std::cout << "P+ = ";
  printPoint(pPlus);
  std::cout << '\n';
  std::cout << "P- = ";
  printPoint(pMinus);
  std::cout << "\n\n";

  const VectorFunction3D lorenzFixedPoint = [sigma, r, beta](Dual x, Dual y,
                                                             Dual z) {
    return StaticVector<Dual, 3>{sigma * (y - x), r * x - y - x * z,
                                 x * y - beta * z};
  };

  const std::array<StaticVector<double, 3>, 3> guesses{
      StaticVector<double, 3>{1.0, 1.0, 1.0},
      StaticVector<double, 3>{10.0, 10.0, 10.0},
      StaticVector<double, 3>{-10.0, -10.0, -10.0}};

  std::cout << "Newton results:\n";
  for (Index i = 0; i < static_cast<Index>(guesses.size()); ++i) {
    auto root = newtonRaphsonSystemEquations(lorenzFixedPoint, guesses[i]);
    std::cout << "guess " << (i + 1) << " = ";
    printPoint(guesses[i]);
    std::cout << " -> root = ";
    printPoint(root);
    std::cout << '\n';
  }

  return 0;
}
