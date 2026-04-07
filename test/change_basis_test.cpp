#include "../library/Coordinates.h"
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
  {
    using T = double;

    // B = {e1, e2}
    const StaticVector<StaticVector<T, 2>, 2> basisB{
        StaticVector<T, 2>{1.0, 0.0},
        StaticVector<T, 2>{0.0, 1.0},
    };

    // C = {v1, v2}, v1=(3,1), v2=(2,1)
    const StaticVector<StaticVector<T, 2>, 2> basisC{
        StaticVector<T, 2>{3.0, 1.0},
        StaticVector<T, 2>{2.0, 1.0},
    };

    const CoordSystem<T, 2> B{basisB};
    const CoordSystem<T, 2> C{basisC};

    // T(x,y) = (x+y, y)  => [T]_B
    // T(e1)=(1,0), T(e2)=(1,1) => columns of [T]_B
    const Matrix<T, 2, 2> TB{{1.0, 1.0, 0.0, 1.0}};

    // [T]_C = P_{C<-B} [T]_B P_{B<-C}
    const Matrix<T, 2, 2> TC = B.newBasisOperator(TB, C);

    // Expected
    const Matrix<T, 2, 2> expectedTB{{1.0, 1.0, 0.0, 1.0}};
    const Matrix<T, 2, 2> expectedTC{{2.0, 1.0, -1.0, 0.0}};

    assert(TB == expectedTB);
    assert(TC == expectedTC);

    std::cout << "[T]_B =\n" << TB << '\n';
    std::cout << "[T]_C =\n" << TC << '\n';

    std::cout << "Problem test passed: [T]_B and [T]_C are correct\n";
  }

  {
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

    const CoordSystem<T, 2> B{basisB};
    const CoordSystem<T, 2> C{basisC};

    assert(B.isOrthogonalBasis());
    assert(B.isOrthonormalBasis());
    assert(C.isOrthogonalBasis());
    assert(C.isOrthonormalBasis());

    const StaticVector<StaticVector<T, 2>, 2> skewBasisB{
        StaticVector<T, 2>{1.0, 0.0},
        StaticVector<T, 2>{1.0, 1.0},
    };
    const CoordSystem<T, 2> SkewB{skewBasisB};
    assert(!SkewB.isOrthogonalBasis());
    assert(!SkewB.isOrthonormalBasis());

    // 1) Verify change-basis matrix formula P_{C<-B} = C^{-1} * B.
    const Matrix<T, 2, 2> P = transitionMatrix(B, C);
    const Matrix<T, 2, 2> expectedP{{0.0, 1.0, -1.0, 0.0}};
    assert(almostEqualMat(P, expectedP));

    // 2) Convert coordinates B -> C and back C -> B.
    const StaticVector<T, 2> vB{2.0, 3.0};
    const StaticVector<T, 2> vC = B.newBasisCoord(vB, C);
    const StaticVector<T, 2> expectedVC{3.0, -2.0};
    assert(almostEqualVec(vC, expectedVC));

    const StaticVector<T, 2> vBRecovered = B.oldBasisCoord(vC, C);
    assert(almostEqualVec(vBRecovered, vB));

    // 3) Linear map matrix change: [T]_C = P [T]_B P^{-1}
    // Let [T]_B = diag(2, 3).
    const Matrix<T, 2, 2> TB{{2.0, 0.0, 0.0, 3.0}};
    const Matrix<T, 2, 2> TC = B.newBasisOperator(TB, C);
    const Matrix<T, 2, 2> expectedTC = P * TB * P.inverse();
    assert(almostEqualMat(TC, expectedTC));

    const Matrix<T, 2, 2> TBRecovered = B.oldBasisOperator(TC, C);
    assert(almostEqualMat(TBRecovered, TB));

    // Consistency check on a vector: TC*vC should equal coordinates in C of
    // T(v)
    const StaticVector<T, 2> TvB = TB * vB;
    const StaticVector<T, 2> TvCFromCoords = B.newBasisCoord(TvB, C);
    const StaticVector<T, 2> TvCFromTC = TC * vC;
    assert(almostEqualVec(TvCFromTC, TvCFromCoords));

    std::cout << "2D change_basis_test passed\n";
  }

  {
    using T = double;

    // Old basis B = Cartesian basis -> A_B = I
    const StaticVector<StaticVector<T, 3>, 3> basisB{
        StaticVector<T, 3>{1.0, 0.0, 0.0},
        StaticVector<T, 3>{0.0, 1.0, 0.0},
        StaticVector<T, 3>{0.0, 0.0, 1.0},
    };

    // New basis C = 90 degree CCW rotation around z-axis in old coordinates:
    // c1 = (0, 1, 0), c2 = (-1, 0, 0), c3 = (0, 0, 1)
    const StaticVector<StaticVector<T, 3>, 3> basisC{
        StaticVector<T, 3>{0.0, 1.0, 0.0},
        StaticVector<T, 3>{-1.0, 0.0, 0.0},
        StaticVector<T, 3>{0.0, 0.0, 1.0},
    };

    const CoordSystem<T, 3> B{basisB};
    const CoordSystem<T, 3> C{basisC};

    assert(B.isOrthogonalBasis());
    assert(B.isOrthonormalBasis());
    assert(C.isOrthogonalBasis());
    assert(C.isOrthonormalBasis());

    const StaticVector<StaticVector<T, 3>, 3> skewBasisB{
        StaticVector<T, 3>{1.0, 0.0, 0.0},
        StaticVector<T, 3>{1.0, 1.0, 0.0},
        StaticVector<T, 3>{0.0, 0.0, 1.0},
    };
    const CoordSystem<T, 3> SkewB{skewBasisB};
    assert(!SkewB.isOrthogonalBasis());
    assert(!SkewB.isOrthonormalBasis());

    // 1) Verify change-basis matrix formula P_{C<-B} = C^{-1} * B.
    const Matrix<T, 3, 3> P = transitionMatrix(B, C);
    const Matrix<T, 3, 3> expectedP{
        {0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0}};
    assert(almostEqualMat(P, expectedP));

    // 2) Convert coordinates B -> C and back C -> B.
    const StaticVector<T, 3> vB{2.0, 3.0, 5.0};
    const StaticVector<T, 3> vC = B.newBasisCoord(vB, C);
    const StaticVector<T, 3> expectedVC{3.0, -2.0, 5.0};
    assert(almostEqualVec(vC, expectedVC));

    const StaticVector<T, 3> vBRecovered = B.oldBasisCoord(vC, C);
    assert(almostEqualVec(vBRecovered, vB));

    // 3) Linear map matrix change: [T]_C = P [T]_B P^{-1}
    // Let [T]_B = diag(2, 3, 4).
    const Matrix<T, 3, 3> TB{{2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0}};
    const Matrix<T, 3, 3> TC = B.newBasisOperator(TB, C);
    const Matrix<T, 3, 3> expectedTC = P * TB * P.inverse();
    assert(almostEqualMat(TC, expectedTC));

    const Matrix<T, 3, 3> TBRecovered = B.oldBasisOperator(TC, C);
    assert(almostEqualMat(TBRecovered, TB));

    // 4) Consistency check on a vector: TC*vC should equal coordinates in C of
    // T(v)
    const StaticVector<T, 3> TvB = TB * vB;
    const StaticVector<T, 3> TvCFromCoords = B.newBasisCoord(TvB, C);
    const StaticVector<T, 3> TvCFromTC = TC * vC;
    assert(almostEqualVec(TvCFromTC, TvCFromCoords));

    // Additional test: test with a non-diagonal matrix
    const Matrix<T, 3, 3> TB2{{1.0, 2.0, 3.0, 0.5, 1.5, 2.5, 0.1, 0.2, 0.3}};
    const Matrix<T, 3, 3> TC2 = B.newBasisOperator(TB2, C);
    const StaticVector<T, 3> TvB2 = TB2 * vB;
    const StaticVector<T, 3> TvCFromCoords2 = B.newBasisCoord(TvB2, C);
    const StaticVector<T, 3> TvCFromTC2 = TC2 * vC;
    assert(almostEqualVec(TvCFromTC2, TvCFromCoords2));

    std::cout << "All 3D change_basis tests passed!\n";
  }
  return 0;
}