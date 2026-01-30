#ifndef PARENT_ELEMENT_SHAPE_FUNCTION
#define PARENT_ELEMENT_SHAPE_FUNCTION
#include "gaussQuadrature.h"
#include <functional>
// Range of 1D parent element Sr: [-1,1] -> mapping to [x1,x2] physical element
// First order polynomial shape function
auto N1_r = [](auto psy) { return (1 + psy) / 2; };
auto N2_r = [](auto psy) { return (1 - psy) / 2; };

// Discretizing x = phi(psy) = x1 N1(psy) + x2 N2(psy) = J*psy + M = psy *
// (x2-x1)/2 + (x2+x1)/2

// Jacobian: J = dx/d(psy) = (x2-x1)/2
template <typename T> T Jacobian(T x1, T x2) { return (x2 - x1) / 2; }
template <typename T> T midPoint(T x1, T x2) { return (x2 + x1) / 2; }

// Mapping deriviation: dN/dx = dN/d(psy) * d(psy)/dx = dN/d(psy) * 1/J
template <typename T> T dN1(T x1, T x2) { return 0.5 / Jacobian(x1, x2); }
template <typename T> T dN2(T x1, T x2) { return -0.5 / Jacobian(x1, x2); }

// Function Phi(psy) used to map
template <typename T> T physicCoor(T psy, T x1, T x2) {
  return N1_r(psy) * x1 + N2_r(psy) * x2;
}
template <typename T> T parentCoor(T x, T x1, T x2) {
  return (x - midPoint(x1, x2)) / Jacobian(x1, x2);
}

// Interpolation: int^(x2)_(x1) f(x)dx = int^1_(-1) f(phi(psy)) J d(psy)
double integrationGauss1D_ref(double x1, double x2,
                              std::function<double(double)> f_psy, int n = 2) {
  using namespace gaussQuadrature;

  const double *xi_ptr = nullptr;
  const double *w_ptr = nullptr;
  valueXiW(n, xi_ptr, w_ptr);
  double I = 0.0;
  double J = Jacobian(x1, x2);
  for (int i = 0; i < n; ++i) {
    I += J * w_ptr[i] * f_psy(xi_ptr[i]);
  }
  return I;
}

#endif

#include <iomanip>
#include <iostream>
int main() {
  double x = 4.0, x1 = 2.0, x2 = 5.0, psy = 0.0;
  std::cout << "=== Mapping Tests ===" << '\n';
  std::cout << "Ref coordinate of x=4: psy = " << parentCoor(x, x1, x2) << '\n';
  std::cout << "Physical coordinate of psy=0: x = " << physicCoor(psy, x1, x2)
            << '\n';
  std::cout << "Jacobian: " << Jacobian(x1, x2) << '\n';
  std::cout << "Derivatives of N1: dN1 = " << dN1(x1, x2) << '\n';
  std::cout << "Derivatives of N2: dN2 = " << dN2(x1, x2) << '\n';

  std::cout << "\n=== Integration Tests ===" << '\n';
  // Test 1: Integrate N1 over element [x1,x2]
  // Analytical: int N1 dpsy = int (1+psy)/2 dpsy from -1 to 1 = [psy/2 +
  // psy^2/4] = 1 With Jacobian J: Result = J * 1 = 1.5
  auto f1 = [](double psy) { return N1_r(psy); };
  double I1 = integrationGauss1D_ref(x1, x2, f1, 2);
  std::cout << std::fixed << std::setprecision(10);
  std::cout << "Integral of N1: " << I1 << " (expected: 1.5)" << '\n';

  // Test 2: Integrate N1*N2
  // Analytical: N1*N2 = (1+psy)/2 * (1-psy)/2 = (1-psy^2)/4
  // int (1-psy^2)/4 dpsy from -1 to 1 = 1/4 * [psy - psy^3/3]|_{-1}^{1}
  // = 1/4 * [(1 - 1/3) - (-1 + 1/3)] = 1/4 * [2/3 + 2/3] = 1/4 * 4/3 = 1/3
  // With Jacobian J: Result = J * 1/3 = 1.5/3 = 0.5
  auto f2 = [](double psy) { return N1_r(psy) * N2_r(psy); };
  double I2 = integrationGauss1D_ref(x1, x2, f2, 2);
  std::cout << "Integral of N1*N2: " << I2 << " (expected: 0.5)" << '\n';

  // Test 3: Integrate constant 1
  // Analytical: int 1 dpsy from -1 to 1 = 2
  // With Jacobian J: Result = J * 2 = 1.5 * 2 = 3 (which is element length
  // x2-x1)
  auto f3 = [](double psy) { return 1.0; };
  double I3 = integrationGauss1D_ref(x1, x2, f3, 2);
  std::cout << "Integral of 1: " << I3 << " (expected: 3.0 = element length)"
            << '\n';

  std::cout << "\n=== Comparison: Mapping in Integrand vs Parent Element ==="
            << '\n';

  // Example: Integrate x^2 over [x1=2, x2=5]
  // Analytical result: int x^2 dx from 2 to 5 = [x^3/3] = 125/3 - 8/3 = 117/3 =
  // 39

  // Method 1: MAPPING IN INTEGRAND (traditional approach)
  // Define f(x) on physical domain, map x=phi(psy) inside integrand
  auto f_physical = [=](double x) {
    return x * x; // f(x) = x^2 defined on physical [2,5]
  };
  // Use standard Gauss integration - it does mapping automatically
  double I_map = integrationGauss1D(x1, x2, f_physical, 2);
  std::cout << "Method 1 - Mapping in integrand:" << '\n';
  std::cout << "  f(x) = x^2 defined on [x1,x2]" << '\n';
  std::cout
      << "  integrationGauss1D maps x = phi(psy) = J*psy + M automatically"
      << '\n';
  std::cout << "  Result: " << I_map << " (expected: 39)" << '\n';

  // Method 2: PARENT ELEMENT (reference domain approach)
  // Define f_psy(psy) directly on reference domain [-1,1]
  // Must explicitly write f(phi(psy)) ourselves
  auto f_ref = [=](double psy) {
    double x = physicCoor(psy, x1, x2); // x = phi(psy) = J*psy + M
    return x * x;                       // f_psy(psy) = (phi(psy))^2
  };
  double I_ref = integrationGauss1D_ref(x1, x2, f_ref, 2);
  std::cout << "\nMethod 2 - Parent element:" << '\n';
  std::cout << "  f_psy(psy) = (phi(psy))^2 defined on [-1,1]" << '\n';
  std::cout << "  We explicitly compute x = phi(psy) inside f_psy" << '\n';
  std::cout << "  integrationGauss1D_ref only multiplies Jacobian" << '\n';
  std::cout << "  Result: " << I_ref << " (expected: 39)" << '\n';

  // When are they DIFFERENT?
  // For shape functions N_r(psy), parent element is SIMPLER!
  std::cout << "\n=== Why use parent element for FEM? ===" << '\n';
  std::cout << "When integrand is N1_r(psy) * N2_r(psy):" << '\n';

  // Parent element: Direct and clean
  auto integrand_parent = [](double psy) {
    return N1_r(psy) * N2_r(psy); // Already on [-1,1]!
  };
  double I_shape_parent = integrationGauss1D_ref(x1, x2, integrand_parent, 2);
  std::cout << "  Parent: f_psy = N1_r(psy)*N2_r(psy)  [SIMPLE!]" << '\n';
  std::cout << "  Result: " << I_shape_parent << '\n';

  // Mapping approach: Need inverse mapping!
  auto integrand_map = [=](double x) {
    double psy = parentCoor(x, x1, x2); // Need to map x -> psy!
    return N1_r(psy) * N2_r(psy);
  };
  double I_shape_map = integrationGauss1D(x1, x2, integrand_map, 2);
  std::cout
      << "  Mapping: f(x) = N1_r(parentCoor(x))*N2_r(parentCoor(x))  [COMPLEX!]"
      << '\n';
  std::cout << "  Result: " << I_shape_map << '\n';

  std::cout << "\nConclusion: For shape functions defined on [-1,1]," << '\n';
  std::cout << "parent element approach is MORE NATURAL and EFFICIENT!" << '\n';

  std::cout << "\n\n=== DETAILED STEP-BY-STEP EXAMPLE ===" << '\n';
  std::cout << "Problem: Integrate N1*N2 over element [x1=2, x2=5] using 2-point Gauss" << '\n';
  std::cout << "----------------------------------------------------------------------" << '\n';
  
  double J_example = Jacobian(x1, x2);
  double M_example = midPoint(x1, x2);
  std::cout << "\nSetup:" << '\n';
  std::cout << "  Element: [x1=" << x1 << ", x2=" << x2 << "]" << '\n';
  std::cout << "  Jacobian: J = (x2-x1)/2 = " << J_example << '\n';
  std::cout << "  Midpoint: M = (x2+x1)/2 = " << M_example << '\n';
  std::cout << "  Mapping: x = phi(psy) = J*psy + M = " << J_example << "*psy + " << M_example << '\n';
  std::cout << "  Inverse: psy = (x - M)/J = (x - " << M_example << ")/" << J_example << '\n';
  
  // 2-point Gauss quadrature points and weights
  std::cout << "\n2-point Gauss quadrature on [-1,1]:" << '\n';
  std::cout << "  psy1 = -1/sqrt(3) = " << -1.0/sqrt(3.0) << ",  w1 = 1.0" << '\n';
  std::cout << "  psy2 = +1/sqrt(3) = " << 1.0/sqrt(3.0) << ",  w2 = 1.0" << '\n';
  
  std::cout << "\n\n--- METHOD 1: integrationGauss1D (Mapping in Integrand) ---" << '\n';
  std::cout << "Formula: I = sum w_i * f(x_i) * J  where x_i = J*psy_i + M" << '\n';
  std::cout << "         (function automatically maps psy to x inside)" << '\n';
  
  double psy1 = -1.0/sqrt(3.0);
  double psy2 = 1.0/sqrt(3.0);
  double w1 = 1.0;
  double w2 = 1.0;
  
  double x_i1 = J_example * psy1 + M_example;
  double x_i2 = J_example * psy2 + M_example;
  
  std::cout << "\nStep 1: Map Gauss points from [-1,1] to [" << x1 << "," << x2 << "]" << '\n';
  std::cout << "  x1 = J*psy1 + M = " << J_example << "*(" << psy1 << ") + " << M_example << " = " << x_i1 << '\n';
  std::cout << "  x2 = J*psy2 + M = " << J_example << "*(" << psy2 << ") + " << M_example << " = " << x_i2 << '\n';
  
  std::cout << "\nStep 2: At each x_i, need to evaluate N1(x)*N2(x)" << '\n';
  std::cout << "  But N1, N2 are defined on [-1,1]! Need inverse map:" << '\n';
  
  double psy_back1 = parentCoor(x_i1, x1, x2);
  double psy_back2 = parentCoor(x_i2, x1, x2);
  std::cout << "  At x1=" << x_i1 << ": psy = (x-M)/J = (" << x_i1 << "-" << M_example << ")/" << J_example << " = " << psy_back1 << '\n';
  std::cout << "  At x2=" << x_i2 << ": psy = (x-M)/J = (" << x_i2 << "-" << M_example << ")/" << J_example << " = " << psy_back2 << '\n';
  
  double N1_at_psy1 = N1_r(psy_back1);
  double N2_at_psy1 = N2_r(psy_back1);
  double f_at_x1 = N1_at_psy1 * N2_at_psy1;
  
  double N1_at_psy2 = N1_r(psy_back2);
  double N2_at_psy2 = N2_r(psy_back2);
  double f_at_x2 = N1_at_psy2 * N2_at_psy2;
  
  std::cout << "\nStep 3: Evaluate N1(psy)*N2(psy) at each point" << '\n';
  std::cout << "  At psy=" << psy_back1 << ":" << '\n';
  std::cout << "    N1 = (1+psy)/2 = (1+" << psy_back1 << ")/2 = " << N1_at_psy1 << '\n';
  std::cout << "    N2 = (1-psy)/2 = (1-" << psy_back1 << ")/2 = " << N2_at_psy1 << '\n';
  std::cout << "    f = N1*N2 = " << f_at_x1 << '\n';
  
  std::cout << "  At psy=" << psy_back2 << ":" << '\n';
  std::cout << "    N1 = (1+psy)/2 = (1+" << psy_back2 << ")/2 = " << N1_at_psy2 << '\n';
  std::cout << "    N2 = (1-psy)/2 = (1-" << psy_back2 << ")/2 = " << N2_at_psy2 << '\n';
  std::cout << "    f = N1*N2 = " << f_at_x2 << '\n';
  
  std::cout << "\nStep 4: Sum with weights and Jacobian" << '\n';
  std::cout << "  I = J * (w1*f1 + w2*f2)" << '\n';
  std::cout << "    = " << J_example << " * (" << w1 << "*" << f_at_x1 << " + " << w2 << "*" << f_at_x2 << ")" << '\n';
  std::cout << "    = " << J_example << " * " << (w1*f_at_x1 + w2*f_at_x2) << '\n';
  std::cout << "    = " << J_example * (w1*f_at_x1 + w2*f_at_x2) << '\n';
  
  std::cout << "\n\n--- METHOD 2: integrationGauss1D_ref (Parent Element) ---" << '\n';
  std::cout << "Formula: I = sum w_i * f_psy(psy_i) * J" << '\n';
  std::cout << "         (function directly uses psy, no mapping needed!)" << '\n';
  
  std::cout << "\nStep 1: Use Gauss points DIRECTLY on [-1,1]" << '\n';
  std::cout << "  psy1 = " << psy1 << '\n';
  std::cout << "  psy2 = " << psy2 << '\n';
  std::cout << "  No mapping needed!" << '\n';
  
  std::cout << "\nStep 2: Evaluate N1(psy)*N2(psy) DIRECTLY at psy_i" << '\n';
  double f_psy1 = N1_r(psy1) * N2_r(psy1);
  double f_psy2 = N1_r(psy2) * N2_r(psy2);
  
  std::cout << "  At psy=" << psy1 << ":" << '\n';
  std::cout << "    N1 = (1+psy)/2 = (1+" << psy1 << ")/2 = " << N1_r(psy1) << '\n';
  std::cout << "    N2 = (1-psy)/2 = (1-" << psy1 << ")/2 = " << N2_r(psy1) << '\n';
  std::cout << "    f_psy = N1*N2 = " << f_psy1 << '\n';
  
  std::cout << "  At psy=" << psy2 << ":" << '\n';
  std::cout << "    N1 = (1+psy)/2 = (1+" << psy2 << ")/2 = " << N1_r(psy2) << '\n';
  std::cout << "    N2 = (1-psy)/2 = (1-" << psy2 << ")/2 = " << N2_r(psy2) << '\n';
  std::cout << "    f_psy = N1*N2 = " << f_psy2 << '\n';
  
  std::cout << "\nStep 3: Sum with weights and Jacobian" << '\n';
  std::cout << "  I = J * (w1*f_psy1 + w2*f_psy2)" << '\n';
  std::cout << "    = " << J_example << " * (" << w1 << "*" << f_psy1 << " + " << w2 << "*" << f_psy2 << ")" << '\n';
  std::cout << "    = " << J_example << " * " << (w1*f_psy1 + w2*f_psy2) << '\n';
  std::cout << "    = " << J_example * (w1*f_psy1 + w2*f_psy2) << '\n';
  
  std::cout << "\n\n=== COMPARISON ===" << '\n';
  std::cout << "Method 1: Map psy->x, then compute f(x), then map x->psy inside f(x)  [COMPLEX]" << '\n';
  std::cout << "Method 2: Use psy directly to compute f(psy)  [SIMPLE]" << '\n';
  std::cout << "\nBoth give same result: " << J_example * (w1*f_psy1 + w2*f_psy2) << '\n';
  std::cout << "But Method 2 is CLEANER and FASTER for FEM!" << '\n';
}