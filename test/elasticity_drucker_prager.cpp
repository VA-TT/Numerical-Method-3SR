#include "../library/Elasticity.h"

int main() {
  const double E = 210e9;
  const double nu = 0.3;

  const auto D = elasticityMatrix(E, nu, "planeStrain");

  Vector<double> stress_n{0.0, 0.0, 0.0};
  Vector<double> strain_n{0.0, 0.0, 0.0};
  Vector<double> dstrain{1e-6, 0.0, 0.0};

  // Example DP parameters
  const double alpha = 0.1;
  const double k = 1.0;

  auto [stress_np1, strain_np1, dl] = updateStressStrainDruckerPrager(
      stress_n, strain_n, dstrain, D, alpha, k, "planeStrain", nu);

  (void)stress_np1;
  (void)strain_np1;
  (void)dl;
  return 0;
}
