#include "../library/Matrix.h"

int main() {
  // B = [ [1,2,4],
  //       [2,1,3],
  //       [4,3,1] ]
  Matrix<double, 3, 3> B{1, 2, 4, 2, 1, 3, 4, 3, 1};

  std::cout << "B =\n" << B << "\n";

  // Eigenvalues (useful sanity-check for sqrt / inverse via spectral power)
  const auto [eigvals, V] = B.eigen(true);
  (void)V;
  std::cout << "eig(B) = " << eigvals << "\n\n";

  const auto B2 = B * B;
  std::cout << "B^2 = B * B =\n" << B2 << "\n";

  // Optional: validate against spectral power for symmetric matrix.
  const auto B2_sym = B.symPower(2.0);
  std::cout << "B^2 from symPower(2) =\n" << B2_sym << "\n";

  std::cout << "Difference (B*B - symPower) =\n" << (B2 - B2_sym) << "\n";

  // Inverse
  const auto I = Matrix<double, 3, 3>::identity();
  const auto Binv = B.inverse();
  std::cout << "B^{-1} from inverse() =\n" << Binv << "\n";
  std::cout << "B * B^{-1} - I =\n" << (B * Binv - I) << "\n";
  std::cout << "B^{-1} * B - I =\n" << (Binv * B - I) << "\n";

  const auto Binv_sym = B.symPower(-1.0);
  std::cout << "B^{-1} from symPower(-1) =\n" << Binv_sym << "\n";
  std::cout << "Difference (inverse() - symPower(-1)) =\n"
            << (Binv - Binv_sym) << "\n";
  std::cout << "B * symPower(-1) - I =\n" << (B * Binv_sym - I) << "\n";

  // Square root
  bool hasNegativeEigenvalue = false;
  for (Index i = 0; i < eigvals.size(); ++i) {
    if (eigvals[i] < 0.0) {
      hasNegativeEigenvalue = true;
      break;
    }
  }
  if (hasNegativeEigenvalue) {
    std::cout << "sqrt(B): not defined over real numbers (B is not PSD)\n\n";
  } else {
    const auto Bhalf = B.symPower(0.5);
    std::cout << "sqrt(B) from symPower(0.5) =\n" << Bhalf << "\n";
    std::cout << "sqrt(B)*sqrt(B) - B =\n" << (Bhalf * Bhalf - B) << "\n";
  }

  // Fourth power
  const auto B4 = B2 * B2;
  std::cout << "B^4 = (B^2)*(B^2) =\n" << B4 << "\n";

  const auto B4_sym = B.symPower(4.0);
  std::cout << "B^4 from symPower(4) =\n" << B4_sym << "\n";
  std::cout << "Difference (direct - symPower) =\n" << (B4 - B4_sym) << "\n";

  return 0;
}
