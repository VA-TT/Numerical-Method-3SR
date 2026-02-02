#include "../library/Matrix.h"
#include "../library/Vector.h"
#include <iomanip>
#include <iostream>

int main() {
  Matrix<double, 2, 2> J;
  J(0, 0) = 1.0;
  J(0, 1) = 0.0;
  J(1, 0) = 0.0;
  J(1, 1) = 0.5;

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "Original J:\n";
  std::cout << "  [" << J(0, 0) << ", " << J(0, 1) << "]\n";
  std::cout << "  [" << J(1, 0) << ", " << J(1, 1) << "]\n";

  Matrix<double, 2, 2> J_inv = J.inverse();
  std::cout << "\nJ_inv:\n";
  std::cout << "  [" << J_inv(0, 0) << ", " << J_inv(0, 1) << "]\n";
  std::cout << "  [" << J_inv(1, 0) << ", " << J_inv(1, 1) << "]\n";

  // Test matrix multiplication
  Matrix<double, 2, 4> dN;
  dN(0, 0) = -0.25;
  dN(0, 1) = 0.25;
  dN(0, 2) = 0.25;
  dN(0, 3) = -0.25;
  dN(1, 0) = -0.25;
  dN(1, 1) = -0.25;
  dN(1, 2) = 0.25;
  dN(1, 3) = 0.25;

  std::cout << "\ndN_parent:\n";
  for (int i = 0; i < 2; ++i) {
    std::cout << "  [";
    for (int j = 0; j < 4; ++j) {
      std::cout << dN(i, j);
      if (j < 3)
        std::cout << ", ";
    }
    std::cout << "]\n";
  }

  Matrix<double, 2, 4> result = J_inv * dN;
  std::cout << "\nJ_inv * dN_parent:\n";
  for (int i = 0; i < 2; ++i) {
    std::cout << "  [";
    for (int j = 0; j < 4; ++j) {
      std::cout << result(i, j);
      if (j < 3)
        std::cout << ", ";
    }
    std::cout << "]\n";
  }

  return 0;
}
