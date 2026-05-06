#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "../library/Matrix.h"
#include "../library/Vector.h"

constexpr double TOLERANCE = 1e-8;

bool approximatelyEqual(double a, double b, double tol = TOLERANCE) {
  return std::abs(a - b) < tol;
}

// ============================================================================
// Test Case 1: Linear Regression (y = 2x + 3)
// ============================================================================
void test_leastSquare_linear() {
  std::cout << "=== Test 1: Linear Regression y = 2x + 3 ===" << std::endl;

  // Generate data points: y = 2*x + 3
  StaticVector<double, 5> x{0.0, 1.0, 2.0, 3.0, 4.0};
  StaticVector<double, 5> y{3.0, 5.0, 7.0, 9.0, 11.0};

  // Perform least squares fitting
  StaticVector<double, 2> beta = leastSquare(x, y);

  // beta[0] = intercept (a), beta[1] = slope (b)
  // Expected: a ≈ 3.0, b ≈ 2.0
  double a = beta[0];
  double b = beta[1];

  std::cout << "Fitted equation: y = " << b << "x + " << a << std::endl;
  std::cout << "Expected: y = 2x + 3" << std::endl;

  bool a_correct = approximatelyEqual(a, 3.0, 1e-6);
  bool b_correct = approximatelyEqual(b, 2.0, 1e-6);

  std::cout << "Intercept a = " << a << " (expected 3.0) "
            << (a_correct ? "✓" : "✗") << std::endl;
  std::cout << "Slope b = " << b << " (expected 2.0) "
            << (b_correct ? "✓" : "✗") << std::endl;

  // Verify predictions
  std::cout << "\nVerification at test points:" << std::endl;
  int correct_predictions = 0;
  for (std::ptrdiff_t i = 0; i < 5; ++i) {
    double y_pred = b * x[i] + a;
    double error = std::abs(y_pred - y[i]);
    bool correct = approximatelyEqual(y_pred, y[i], 1e-10);
    std::cout << "  x=" << x[i] << ": y_actual=" << y[i]
              << ", y_pred=" << y_pred << ", error=" << error
              << (correct ? " ✓" : " ✗") << std::endl;
    if (correct)
      correct_predictions++;
  }

  std::cout << "Correct predictions: " << correct_predictions << "/5"
            << (a_correct && b_correct && correct_predictions == 5 ? " ✓ PASS"
                                                                   : " ✗ FAIL")
            << std::endl
            << std::endl;
}

// ============================================================================
// Test Case 2: Data with Noise (y = x + 1 with small perturbations)
// ============================================================================
void test_leastSquare_noisy() {
  std::cout << "=== Test 2: Noisy Data (y ≈ x + 1) ===" << std::endl;

  // Data with small noise
  StaticVector<double, 6> x{0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  StaticVector<double, 6> y{1.05, 2.02, 2.98,
                            4.03, 5.01, 6.04}; // y ≈ x + 1 + noise

  StaticVector<double, 2> beta = leastSquare(x, y);
  double a = beta[0];
  double b = beta[1];

  std::cout << "Fitted equation: y = " << b << "x + " << a << std::endl;
  std::cout << "Expected (approximately): y = 1.0x + 1.0" << std::endl;

  // Check that fitted values are close to expected
  bool a_reasonable = std::abs(a - 1.0) < 0.1;
  bool b_reasonable = std::abs(b - 1.0) < 0.1;

  std::cout << "Intercept a = " << a << " (expected ~1.0) "
            << (a_reasonable ? "✓" : "✗") << std::endl;
  std::cout << "Slope b = " << b << " (expected ~1.0) "
            << (b_reasonable ? "✓" : "✗") << std::endl;

  // Calculate R² (coefficient of determination)
  double mean_y = 0.0;
  for (std::ptrdiff_t i = 0; i < 6; ++i)
    mean_y += y[i];
  mean_y /= 6.0;

  double ss_tot = 0.0, ss_res = 0.0;
  for (std::ptrdiff_t i = 0; i < 6; ++i) {
    double y_pred = b * x[i] + a;
    ss_tot += (y[i] - mean_y) * (y[i] - mean_y);
    ss_res += (y[i] - y_pred) * (y[i] - y_pred);
  }
  double r_squared = 1.0 - ss_res / ss_tot;

  std::cout << "R² = " << r_squared << " (expected > 0.99 for low noise)"
            << (r_squared > 0.99 ? " ✓" : " ✗") << std::endl
            << std::endl;
}

// ============================================================================
// Test Case 3: Horizontal Line (y = 5, constant)
// ============================================================================
void test_leastSquare_constant() {
  std::cout << "=== Test 3: Constant Function (y = 5) ===" << std::endl;

  StaticVector<double, 4> x{0.0, 1.0, 2.0, 3.0};
  StaticVector<double, 4> y{5.0, 5.0, 5.0, 5.0};

  StaticVector<double, 2> beta = leastSquare(x, y);
  double a = beta[0];
  double b = beta[1];

  std::cout << "Fitted equation: y = " << b << "x + " << a << std::endl;
  std::cout << "Expected: y = 0x + 5 (horizontal line)" << std::endl;

  bool a_correct = approximatelyEqual(a, 5.0, 1e-10);
  bool b_correct = approximatelyEqual(b, 0.0, 1e-10);

  std::cout << "Intercept a = " << a << " (expected 5.0) "
            << (a_correct ? "✓" : "✗") << std::endl;
  std::cout << "Slope b = " << b << " (expected 0.0) "
            << (b_correct ? "✓" : "✗") << std::endl
            << std::endl;
}

// ============================================================================
// Test Case 4: Steep Line (y = -3x + 10)
// ============================================================================
void test_leastSquare_steep() {
  std::cout << "=== Test 4: Steep Negative Line (y = -3x + 10) ==="
            << std::endl;

  StaticVector<double, 5> x{0.0, 1.0, 2.0, 3.0, 4.0};
  StaticVector<double, 5> y{10.0, 7.0, 4.0, 1.0, -2.0};

  StaticVector<double, 2> beta = leastSquare(x, y);
  double a = beta[0];
  double b = beta[1];

  std::cout << "Fitted equation: y = " << b << "x + " << a << std::endl;
  std::cout << "Expected: y = -3x + 10" << std::endl;

  bool a_correct = approximatelyEqual(a, 10.0, 1e-6);
  bool b_correct = approximatelyEqual(b, -3.0, 1e-6);

  std::cout << "Intercept a = " << a << " (expected 10.0) "
            << (a_correct ? "✓" : "✗") << std::endl;
  std::cout << "Slope b = " << b << " (expected -3.0) "
            << (b_correct ? "✓" : "✗") << std::endl
            << std::endl;
}

// ============================================================================
// Test Case 5: Plane Fit (z = 1 + 2x + 3y)
// ==========================================================================
void test_leastSquare_plane() {
  std::cout << "=== Test 5: Plane Fit (z = 1 + 2x + 3y) ===" << std::endl;

  StaticVector<double, 6> x{0.0, 1.0, 0.0, 1.0, 2.0, 2.0};
  StaticVector<double, 6> y{0.0, 0.0, 1.0, 1.0, 0.0, 1.0};
  StaticVector<double, 6> z{1.0, 3.0, 4.0, 6.0, 5.0, 8.0};

  StaticVector<double, 3> beta = leastSquare(x, y, z);
  double a = beta[0];
  double b = beta[1];
  double c = beta[2];

  std::cout << "Fitted plane: z = " << a << " + " << b << "x + " << c << "y"
            << std::endl;
  std::cout << "Expected: z = 1 + 2x + 3y" << std::endl;

  bool a_correct = approximatelyEqual(a, 1.0, 1e-6);
  bool b_correct = approximatelyEqual(b, 2.0, 1e-6);
  bool c_correct = approximatelyEqual(c, 3.0, 1e-6);

  std::cout << "Intercept a = " << a << " (expected 1.0) "
            << (a_correct ? "✓" : "✗") << std::endl;
  std::cout << "x coefficient b = " << b << " (expected 2.0) "
            << (b_correct ? "✓" : "✗") << std::endl;
  std::cout << "y coefficient c = " << c << " (expected 3.0) "
            << (c_correct ? "✓" : "✗") << std::endl;

  int correct_predictions = 0;
  for (std::ptrdiff_t i = 0; i < 6; ++i) {
    double z_pred = a + b * x[i] + c * y[i];
    bool correct = approximatelyEqual(z_pred, z[i], 1e-10);
    std::cout << "  x=" << x[i] << ", y=" << y[i] << ": z_actual=" << z[i]
              << ", z_pred=" << z_pred << (correct ? " ✓" : " ✗") << std::endl;
    if (correct)
      ++correct_predictions;
  }

  std::cout << "Correct predictions: " << correct_predictions << "/6"
            << (a_correct && b_correct && c_correct && correct_predictions == 6
                    ? " ✓ PASS"
                    : " ✗ FAIL")
            << std::endl
            << std::endl;
}

int main() {
  std::cout << "\n╔════════════════════════════════════════════════════════╗"
            << "\n║          Least Squares Regression Testing             ║"
            << "\n╚════════════════════════════════════════════════════════╝\n"
            << std::endl;

  test_leastSquare_linear();
  test_leastSquare_noisy();
  test_leastSquare_constant();
  test_leastSquare_steep();
  test_leastSquare_plane();

  std::cout << "╔════════════════════════════════════════════════════════╗"
            << "\n║           All Least Squares Tests Completed!           ║"
            << "\n╚════════════════════════════════════════════════════════╝\n"
            << std::endl;

  return 0;
}
