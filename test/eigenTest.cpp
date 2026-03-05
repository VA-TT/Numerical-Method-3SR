
#include "../library/Matrix.h"

#include <cmath>
#include <iomanip>
#include <iostream>

using Index = std::ptrdiff_t;

template <Index N>
double residualNorm2(const Matrix<double, N, N> &A, const Vector<double> &v,
                     double lambda) {
  const Matrix<double, N, 1> vc(v);
  const Matrix<double, N, 1> rc = (A * vc) - (lambda * vc);
  const Vector<double> r = static_cast<Vector<double>>(rc);
  double s = 0.0;
  for (Index i = 0; i < r.size(); ++i)
    s += r[i] * r[i];
  return std::sqrt(s);
}

template <Index N>
void printEigenReport(const Matrix<double, N, N> &A, const char *label) {
  std::cout << label << std::endl;
  std::cout << A << std::endl;

  std::cout << "Matrix is symmetric: " << (A.isSymmetric() ? "True" : "False")
            << "." << std::endl;
  std::cout << std::endl;

  std::cout << "Computing eigenvalues and eigenvectors using eigen()..."
            << std::endl;
  auto [eigvals, V] = A.eigen();

  std::cout << "Eigenvalues:" << std::endl;
  std::cout << std::setprecision(6) << eigvals << std::endl << std::endl;

  std::cout << "Eigenvector matrix V (columns):\n" << V << std::endl;

  std::cout << "Residual norms ||A v_i - lambda_i v_i||_2:" << std::endl;
  for (Index i = 0; i < N; ++i) {
    const Vector<double> vi = V.getColVector(i);
    std::cout << std::setprecision(6) << residualNorm2(A, vi, eigvals[i])
              << " ";
  }
  std::cout << std::endl << std::endl;

  if (A.isSymmetric()) {
    std::cout << "V^T * V (should be close to identity):\n"
              << (V.transpose() * V) << std::endl;
  } else {
    std::cout << "For non-symmetric A: diagonalization check uses V^{-1} A V"
                 " (not V^T A V)."
              << std::endl;
    try {
      const Matrix<double, N, N> Vinv = V.inverse();
      std::cout << "V^{-1} * A * V (should be close to diagonal if A is "
                   "diagonalizable with real eigenpairs):\n"
                << (Vinv * A * V) << std::endl;
    } catch (const std::exception &e) {
      std::cout << "Could not compute V.inverse(): " << e.what() << std::endl;
    }
  }
  std::cout << std::endl;
}

int main() {
  std::cout << "**********************************************" << std::endl;
  std::cout << "Testing eigenvalue and eigenvector code." << std::endl;
  std::cout << "Using ONLY Matrix::eigen()" << std::endl;
  std::cout << "**********************************************" << std::endl;
  std::cout << std::endl;

  // Non-symmetric example (same as your QR test)
  {
    Matrix<double, 3, 3> A{0.5, 0.75, 0.5, 1.0, 0.5, 0.75, 0.25, 0.25, 0.25};
    printEigenReport(A, "Testing with a simple non-symmetric 3x3 matrix:");
  }

  // Symmetric example (good case: real eigenpairs + orthogonal eigenvectors)
  {
    Matrix<double, 3, 3> A{6.0, 5.5, -1.0, 5.5, 1.0, -2.0, -1.0, -2.0, -3.0};
    printEigenReport(A, "Testing with a simple symmetric 3x3 matrix:");
  }

  // Example that likely has complex eigenvalues (real-only method limitation)
  {
    Matrix<double, 3, 3> A{4.0, -6.0, 8.0, 7.0, 9.0, -5.0, 9.0, -6.0, -4.0};
    printEigenReport(
        A, "Testing with an example that should have complex eigenvalues:");
  }
  {

    Matrix<double, 3, 3> A{1.0, 2.0, 4.0, 1.0, 1.0, 5.0, 1.0, 5.0, 1.0};
    printEigenReport(
        A, "Testing with an example that should have complex eigenvalues:");
  }

  return 0;
}