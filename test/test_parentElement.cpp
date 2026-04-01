#include "../library/parentElement.h"
#include <iomanip>

int main() {
  // Test rectangle element: (0,0), (2,0), (2,1), (0,1)
  DynamicVector<double> x_nodes = {0.0, 2.0, 2.0, 0.0};
  DynamicVector<double> y_nodes = {0.0, 0.0, 1.0, 1.0};

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "=== Test parentElement.h Functions ===" << '\n';

  // Test 1D functions
  std::cout << "\n--- 1D Functions ---\n";
  double x1 = 2.0, x2 = 5.0;
  std::cout << "Element: [" << x1 << ", " << x2 << "]\n";
  std::cout << "Jacobian: " << Jacobian(x1, x2) << '\n';
  std::cout << "physicCoor(0): " << physicCoor(0.0, x1, x2) << '\n';
  std::cout << "parentCoor(3.5): " << parentCoor(3.5, x1, x2) << '\n';

  // Test 2D functions
  std::cout << "\n--- 2D Functions ---\n";
  auto [x_center, y_center] = physicCoor2D(0.0, 0.0, x_nodes, y_nodes);
  std::cout << "physicCoor2D(0,0): (" << x_center << ", " << y_center << ")\n";

  Matrix<double, 2, 2> J = Jacobian2D(0.0, 0.0, x_nodes, y_nodes);
  std::cout << "Jacobian2D at center: det(J) = " << det(J) << '\n';

  auto [dN_dx, dN_dy] = dNdxdy(0.0, 0.0, x_nodes, y_nodes);
  std::cout << "dN/dx = [";
  for (int i = 0; i < 4; ++i) {
    std::cout << dN_dx[i];
    if (i < 3)
      std::cout << ", ";
  }
  std::cout << "]\n";

  std::cout << "dN/dy = [";
  for (int i = 0; i < 4; ++i) {
    std::cout << dN_dy[i];
    if (i < 3)
      std::cout << ", ";
  }
  std::cout << "]\n";

  // Test inverse mapping
  auto [xi, eta] = parentCoor2D(1.0, 0.5, x_nodes, y_nodes);
  std::cout << "\nparentCoor2D(1.0, 0.5): (xi=" << xi << ", eta=" << eta
            << ")\n";

  // Test integration
  auto f_const = [](double xi, double eta) { return 1.0; };
  double area = integrationGauss2D_ref(x_nodes, y_nodes, f_const, 2);
  std::cout << "\nArea (integrationGauss2D_ref): " << area
            << " (expected: 2.0)\n";

  std::cout << "\n✓ All tests completed!" << '\n';
  return 0;
}
