#include <cmath>
#include <iomanip>
#include <iostream>

double square_root(double a, double x0, double epsilon) {
  double xn = 0.5 * (x0 + a / x0); // x1 ≥ √a theo AM-GM
  while ((xn * xn - a) > epsilon)  // Không cần abs() từ x1 trở đi!
  {
    xn = 0.5 * (xn + a / xn);
  }
  return xn;
}

int main() {
  // Error eps = 1e-12
  double x = square_root(2, 1.2, 1e-12);

  // In với độ chính xác cao
  std::cout << std::fixed << std::setprecision(15);
  std::cout << "Result: " << x << std::endl;

  // So sánh với sqrt() chuẩn
  std::cout << "std::sqrt(2): " << std::sqrt(2) << std::endl;

  // Tính sai số
  double error = std::abs(x - std::sqrt(2));
  std::cout << std::scientific << std::setprecision(3);
  std::cout << "Error: " << error << std::endl;

  return 0;
}