#include "../library/coordinates.h"
#include <cassert>
#include <cmath>
#include <iostream>

namespace {

template <typename T, Index n>
bool almostEqualVec(const StaticVector<T, n> &a, const StaticVector<T, n> &b,
                    T tol = static_cast<T>(1e-10)) {
  for (Index i = 0; i < n; ++i) {
    if (std::abs(a[i] - b[i]) > tol) {
      return false;
    }
  }
  return true;
}

template <typename T, Index R, Index C>
bool almostEqualMat(const Matrix<T, R, C> &A, const Matrix<T, R, C> &B,
                    T tol = static_cast<T>(1e-10)) {
  for (Index i = 0; i < R; ++i) {
    for (Index j = 0; j < C; ++j) {
      if (std::abs(A(i, j) - B(i, j)) > tol) {
        return false;
      }
    }
  }
  return true;
}

} // namespace

int main() {
  using T = double;

  // Old basis B = Cartesian basis -> A_B = I
  const StaticVector<StaticVector<T, 2>, 2> basisB{
      StaticVector<T, 2>{1.0, 0.0},
      StaticVector<T, 2>{0.0, 1.0},
  };

  // New basis C = 90 degree CCW rotation in old coordinates:
  // c1 = (0,1), c2 = (-1,0)
  const StaticVector<StaticVector<T, 2>, 2> basisC{
      StaticVector<T, 2>{0.0, 1.0},
      StaticVector<T, 2>{-1.0, 0.0},
  };

  const changeCoor<T, 2> B{basisB};
  const changeCoor<T, 2> C{basisC};

  // 1) Verify change-basis matrix formula P_{C<-B} = C^{-1} * B.
  const Matrix<T, 2, 2> P = changeCoor<T, 2>::changeBasisMatrix(B, C);
  const Matrix<T, 2, 2> expectedP{{0.0, 1.0, -1.0, 0.0}};
  assert(almostEqualMat(P, expectedP));

  // 2) Convert coordinates B -> C and back C -> B.
  const StaticVector<T, 2> vB{2.0, 3.0};
  const StaticVector<T, 2> vC = B.newBasisCoor(vB, C);
  const StaticVector<T, 2> expectedVC{3.0, -2.0};
  assert(almostEqualVec(vC, expectedVC));

  const StaticVector<T, 2> vBRecovered = B.oldBasisCoor(vC, C);
  assert(almostEqualVec(vBRecovered, vB));

  // 3) Linear map matrix change: [T]_C = P [T]_B P^{-1}
  // Let [T]_B = diag(2, 3).
  const Matrix<T, 2, 2> TB{{2.0, 0.0, 0.0, 3.0}};
  const Matrix<T, 2, 2> TC = B.newBasisOperator(TB, C);
  const Matrix<T, 2, 2> expectedTC = P * TB * P.inverse();
  assert(almostEqualMat(TC, expectedTC));

  const Matrix<T, 2, 2> TBRecovered = B.oldBasisOperator(TC, C);
  assert(almostEqualMat(TBRecovered, TB));

  // Consistency check on a vector: TC*vC should equal coordinates in C of T(v)
  const StaticVector<T, 2> TvB = TB * vB;
  const StaticVector<T, 2> TvCFromCoords = B.newBasisCoor(TvB, C);
  const StaticVector<T, 2> TvCFromTC = TC * vC;
  assert(almostEqualVec(TvCFromTC, TvCFromCoords));

  std::cout << "change_basis_test passed\n";
  return 0;
}
