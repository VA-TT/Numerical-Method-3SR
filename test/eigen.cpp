#include "../library/Matrix.h"
#include <iostream>

int main() {
  Matrix<double, 3, 3> A{0.5, 0.75, 0.5, 1.0, 0.5, 0.75, 0.25, 0.25, 0.25};
  std::cout << A.QRdecomposition();
  Matrix<double, 3, 3> B{6.0, 5.5, -1.0, 5.5, 1.0, -2.0, -1.0, -2.0, -3.0};
  //   std::cout << B.QRdecomposition();
  return 0;
}