#include "../library/Matrix.h" //Approximative Comparsion
#include "../library/Vector.h"
#include "iostream"
#include <iomanip>     //tab
#include <stdexcept>   //throw exception
#include <type_traits> // precision

////////////////////////////////////////////////////////////////////////////
///////////////////////////   TEST FUNCTION   //////////////////////////////
////////////////////////////////////////////////////////////////////////////

void testMatrixConstructors() {
  std::cout << "=== Testing Matrix Constructors ===" << std::endl;

  // Default constructor
  Matrix<double, 2, 3> m1;
  std::cout << "Default constructor (2x3):\n" << m1 << std::endl;

  // Initializer list constructor
  Matrix<int, 2, 3> m2{1, 2, 3, 4, 5, 6};
  std::cout << "Initializer list constructor {1,2,3,4,5,6}:\n"
            << m2 << std::endl;

  // Copy constructor test
  Matrix<int, 2, 3> m3 = m2;
  std::cout << "Copy constructor:\n" << m3 << std::endl;

  std::cout << "Constructor tests passed!\n" << std::endl;
}

void testMatrixAccessors() {
  std::cout << "=== Testing Matrix Accessors ===" << std::endl;

  Matrix<int, 3, 3> m{1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::cout << "Original matrix:\n" << m << std::endl;

  // Test operator() access
  std::cout << "Element (0,0): " << m(0, 0) << std::endl;
  std::cout << "Element (1,2): " << m(1, 2) << std::endl;
  std::cout << "Element (2,1): " << m(2, 1) << std::endl;

  // Test operator[] access
  std::cout << "Element [0]: " << m[0] << std::endl;
  std::cout << "Element [4]: " << m[4] << std::endl;
  std::cout << "Element [8]: " << m[8] << std::endl;

  // Test modification
  m(1, 1) = 99;
  m[0] = 88;
  std::cout << "After modification m(1,1)=99, m[0]=88:\n" << m << std::endl;

  // Test getters
  std::cout << "Rows: " << m.getRows() << ", Cols: " << m.getCols()
            << ", Length: " << m.length() << std::endl;

  std::cout << "Accessor tests passed!\n" << std::endl;
}

void testMatrixStaticMethods() {
  std::cout << "=== Testing Static Matrix Methods ===" << std::endl;

  // Zero matrix
  auto zero = Matrix<double, 3, 3>::zero();
  std::cout << "Zero matrix (3x3):\n" << zero << std::endl;

  // Ones matrix
  auto ones = Matrix<int, 2, 4>::ones();
  std::cout << "Ones matrix (2x4):\n" << ones << std::endl;

  // Identity matrix
  auto identity = Matrix<double, 4, 4>::identity();
  std::cout << "Identity matrix (4x4):\n" << identity << std::endl;

  // Test reset methods
  Matrix<int, 2, 2> m{1, 2, 3, 4};
  std::cout << "Before reset:\n" << m << std::endl;

  m.resetZero();
  std::cout << "After resetZero():\n" << m << std::endl;

  m.resetOnes();
  std::cout << "After resetOnes():\n" << m << std::endl;

  m.resetIdentity();
  std::cout << "After resetIdentity():\n" << m << std::endl;

  std::cout << "Static methods tests passed!\n" << std::endl;
}

void testMatrixArithmetic() {
  std::cout << "=== Testing Matrix Arithmetic ===" << std::endl;

  Matrix<double, 2, 3> A{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  Matrix<double, 2, 3> B{2.0, 1.0, 4.0, 3.0, 6.0, 5.0};

  std::cout << "Matrix A:\n" << A << std::endl;
  std::cout << "Matrix B:\n" << B << std::endl;

  // Addition
  auto C = A + B;
  std::cout << "A + B:\n" << C << std::endl;

  // Subtraction
  auto D = A - B;
  std::cout << "A - B:\n" << D << std::endl;

  // Scalar multiplication
  auto E = 2.5 * A;
  std::cout << "2.5 * A:\n" << E << std::endl;

  auto F = A * 3.0;
  std::cout << "A * 3.0:\n" << F << std::endl;

  // Unary minus
  auto G = -A;
  std::cout << "-A:\n" << G << std::endl;

  // Equality
  auto A_copy = A;
  std::cout << "A == A_copy: " << (A == A_copy) << std::endl;
  std::cout << "A == B: " << (A == B) << std::endl;
  std::cout << "A != B: " << (A != B) << std::endl;

  std::cout << "Arithmetic tests passed!\n" << std::endl;
}

void testMatrixMultiplication() {
  std::cout << "=== Testing Matrix Multiplication ===" << std::endl;

  Matrix<double, 2, 3> A{1, 2, 3, 4, 5, 6};
  Matrix<double, 3, 2> B{1, 2, 3, 4, 5, 6};

  std::cout << "Matrix A (2x3):\n" << A << std::endl;
  std::cout << "Matrix B (3x2):\n" << B << std::endl;

  auto C = A * B;
  std::cout << "A * B (2x2):\n" << C << std::endl;

  auto D = B * A;
  std::cout << "B * A (3x3):\n" << D << std::endl;

  // Square matrix multiplication
  Matrix<int, 3, 3> E{1, 2, 3, 0, 1, 4, 5, 6, 0};
  Matrix<int, 3, 3> F{1, 0, 0, 0, 1, 0, 0, 0, 1}; // Identity

  std::cout << "Matrix E:\n" << E << std::endl;
  std::cout << "Matrix F (Identity):\n" << F << std::endl;

  auto G = E * F;
  std::cout << "E * Identity:\n" << G << std::endl;

  std::cout << "Matrix multiplication tests passed!\n" << std::endl;
}

void testMatrixProperties() {
  std::cout << "=== Testing Matrix Properties ===" << std::endl;

  // Diagonal matrix
  Matrix<int, 3, 3> diagonal{5, 0, 0, 0, 3, 0, 0, 0, 7};
  std::cout << "Diagonal matrix:\n" << diagonal << std::endl;
  std::cout << "isDiagonal(): " << diagonal.isDiagonal() << std::endl;
  std::cout << "isUpperTriangular(): " << diagonal.isUpperTriangular()
            << std::endl;
  std::cout << "isLowerTriangular(): " << diagonal.isLowerTriangular()
            << std::endl;
  std::cout << "isSquare(): " << diagonal.isSquare() << std::endl;

  // Upper triangular
  Matrix<int, 3, 3> upper{1, 2, 3, 0, 4, 5, 0, 0, 6};
  std::cout << "\nUpper triangular matrix:\n" << upper << std::endl;
  std::cout << "isUpperTriangular(): " << upper.isUpperTriangular()
            << std::endl;
  std::cout << "isLowerTriangular(): " << upper.isLowerTriangular()
            << std::endl;

  // Symmetric matrix
  Matrix<double, 3, 3> symmetric{1, 2, 3, 2, 4, 5, 3, 5, 6};
  std::cout << "\nSymmetric matrix:\n" << symmetric << std::endl;
  std::cout << "isSymmetric(): " << symmetric.isSymmetric() << std::endl;

  // Non-square matrix
  Matrix<int, 2, 3> nonsquare{1, 2, 3, 4, 5, 6};
  std::cout << "\nNon-square matrix (2x3):\n" << nonsquare << std::endl;
  std::cout << "isSquare(): " << nonsquare.isSquare() << std::endl;

  std::cout << "Matrix properties tests passed!\n" << std::endl;
}

void testDeterminantAndTrace() {
  std::cout << "=== Testing Determinant and Trace ===" << std::endl;

  // 1x1 matrix
  Matrix<double, 1, 1> m1{5.0};
  std::cout << "1x1 matrix: " << m1(0, 0) << std::endl;
  std::cout << "det(1x1): " << det(m1) << std::endl;
  std::cout << "trace(1x1): " << trace(m1) << std::endl;

  // 2x2 matrix
  Matrix<double, 2, 2> m2{1, 2, 3, 4};
  std::cout << "\n2x2 matrix:\n" << m2 << std::endl;
  std::cout << "det(2x2): " << det(m2) << std::endl;
  std::cout << "trace(2x2): " << trace(m2) << std::endl;

  // 3x3 matrix
  Matrix<double, 3, 3> m3{1, 2, 3, 0, 1, 4, 5, 6, 0};
  std::cout << "\n3x3 matrix:\n" << m3 << std::endl;
  std::cout << "det(3x3): " << det(m3) << std::endl;
  std::cout << "trace(3x3): " << trace(m3) << std::endl;

  // Identity matrix
  auto identity = Matrix<double, 3, 3>::identity();
  std::cout << "\n3x3 Identity:\n" << identity << std::endl;
  std::cout << "det(Identity): " << det(identity) << std::endl;
  std::cout << "trace(Identity): " << trace(identity) << std::endl;

  std::cout << "Determinant and trace tests passed!\n" << std::endl;
}

void testMatrixInverse() {
  std::cout << "=== Testing Matrix Inverse ===" << std::endl;

  try {
    // 2x2 invertible matrix
    Matrix<double, 2, 2> m2{1, 2, 3, 5};
    std::cout << "Original 2x2 matrix:\n" << m2 << std::endl;
    std::cout << "det(m2): " << det(m2) << std::endl;

    auto inv2 = m2.inverse();
    std::cout << "Inverse:\n" << inv2 << std::endl;

    auto product2 = m2 * inv2;
    std::cout << "Original * Inverse:\n" << product2 << std::endl;

  } catch (const std::exception &e) {
    std::cout << "Caught exception: " << e.what() << std::endl;
  }

  // Test singular matrix (should throw)
  try {
    Matrix<double, 2, 2> singular{1, 2, 2, 4};
    std::cout << "\nTrying to invert singular matrix:\n"
              << singular << std::endl;
    auto inv_singular = singular.inverse();
  } catch (const std::exception &e) {
    std::cout << "✓ Caught expected exception: " << e.what() << std::endl;
  }

  std::cout << "Matrix inverse tests passed!\n" << std::endl;
}

void testRowColumnOperations() {
  std::cout << "=== Testing Row/Column Operations ===" << std::endl;

  Matrix<int, 3, 4> m{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
  std::cout << "Original matrix:\n" << m << std::endl;

  // Test row access
  auto row1 = m.row(1);
  std::cout << "Row 1: ";
  for (const auto &val : row1) {
    std::cout << val << " ";
  }
  std::cout << std::endl;

  // Test column access
  auto col2 = m.col(2);
  std::cout << "Column 2: ";
  for (const auto &val : col2) {
    std::cout << val.get() << " ";
  }
  std::cout << std::endl;

  // Test row swapping
  Matrix<int, 3, 3> swap_test{1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::cout << "Before swapping rows 0 and 2:\n" << swap_test << std::endl;
  swap_test.swapRows(0, 2);
  std::cout << "After swapping rows 0 and 2:\n" << swap_test << std::endl;

  // Test transpose
  Matrix<double, 2, 3> transpose_test{1, 2, 3, 4, 5, 6};
  std::cout << "Original (2x3):\n" << transpose_test << std::endl;
  auto transposed = transpose_test.transpose();
  std::cout << "Transposed (3x2):\n" << transposed << std::endl;

  std::cout << "Row/Column operations tests passed!\n" << std::endl;
}

void testConcatenationAndSplit() {
  std::cout << "=== Testing Concatenation and Split ===" << std::endl;

  Matrix<int, 2, 2> A{1, 2, 3, 4};
  Matrix<int, 2, 3> B{5, 6, 7, 8, 9, 10};

  std::cout << "Matrix A (2x2):\n" << A << std::endl;
  std::cout << "Matrix B (2x3):\n" << B << std::endl;

  auto concat = concatenateMatrixHorizontal(A, B);
  std::cout << "Concatenated matrix A|B (2x5):\n" << concat << std::endl;

  // Test splitByColumn
  Matrix<int, 2, 2> split_A;
  Matrix<int, 2, 3> split_B;
  concat.splitByColumn(2, split_A, split_B);

  std::cout << "Split back - first part:\n" << split_A << std::endl;
  std::cout << "Split back - second part:\n" << split_B << std::endl;

  std::cout << "Concatenation and split tests passed!\n" << std::endl;
}

void testLinearSystem() {
  std::cout << "=== Testing Linear System Solver ===" << std::endl;

  try {
    // Solve Ax = b
    Matrix<double, 3, 3> A{2, 1, 1, 1, 0, 1, 0, 3, 1};
    Matrix<double, 3, 1> b{4, 2, 6};

    std::cout << "Coefficient matrix A:\n" << A << std::endl;
    std::cout << "Right-hand side b:\n" << b << std::endl;

    auto x = solveLinearSystem(A, b);
    std::cout << "Solution x:\n" << x << std::endl;

    // Verify: Ax should equal b
    auto verification = A * x;
    std::cout << "Verification A*x:\n" << verification << std::endl;

  } catch (const std::exception &e) {
    std::cout << "Caught exception: " << e.what() << std::endl;
  }

  std::cout << "Linear system tests passed!\n" << std::endl;
}

void testVectorIntegration() {
  std::cout << "=== Testing Matrix-Vector Integration ===" << std::endl;

  Matrix<double, 3, 3> A{2, 1, 1, 1, 0, 1, 0, 3, 1};
  DynamicVector<double> b{4, 2, 6};

  std::cout << "Matrix A:\n" << A << std::endl;
  std::cout << "Vector b: " << b << std::endl;

  auto x = solveLinearSystem(A, b);
  std::cout << "Solution x:\n" << x << std::endl;

  // Convert Vector back to Matrix for multiplication
  Matrix<double, 3, 1> xMatrix{x};
  auto verification = A * xMatrix;
  std::cout << "Verification A*x:\n" << verification << std::endl;
}

void testTensorProduct() {
  std::cout << "\n=== Testing Tensor Product (Direct Method) ===\n";

  // ✅ FIX: Use flat initializer list
  Matrix<double, 2, 2> A{1, 2, 3, 4}; // Row-major: [1,2; 3,4]
  Matrix<double, 2, 2> B{5, 6, 7, 8}; // Row-major: [5,6; 7,8]

  auto result = tensorProduct(A, B);

  std::cout << "A (2×2):\n" << A << "\n";
  std::cout << "B (2×2):\n" << B << "\n";
  std::cout << "A ⊗ B (4×4):\n" << result << "\n";

  std::cout << "Tensor product test passed!\n";
}

void testTensorProductBlocks() {
  std::cout << "\n=== Testing Tensor Product (Block Method) ===\n";

  Matrix<double, 2, 2> A = {1, 2, 3, 4};
  Matrix<double, 2, 2> B = {5, 6, 7, 8};

  auto result = tensorProduct(A, B);

  std::cout << "A:\n" << A << "\n";
  std::cout << "B:\n" << B << "\n";
  std::cout << "A ⊗ B (4×4):\n" << result << "\n";

  // Expected:
  // | 1*B  2*B | = | 5  6 | 10 12 |
  // | 3*B  4*B |   | 7  8 | 14 16 |
  //                |15 18 | 20 24 |
  //                |21 24 | 28 32 |
}

void testStaticFirstMatrixFeatures() {
  std::cout << "\n=== Testing Static-first Matrix Features ===\n";

  Matrix<double, 2, 2> A{4.0, 1.0, 1.0, 3.0};

  // 1) getColVector now returns StaticVector<T, nRows>
  const auto c0 = A.getColVector(0);
  static_assert(
      std::is_same_v<std::decay_t<decltype(c0)>, StaticVector<double, 2>>,
      "getColVector should return StaticVector<double, 2>");
  if (!(std::abs(c0[0] - 4.0) < 1e-12 && std::abs(c0[1] - 1.0) < 1e-12)) {
    throw std::runtime_error("getColVector static return mismatch");
  }

  // 2) Matrix * StaticVector
  const StaticVector<double, 2> x{2.0, -1.0};
  const auto Ax = A * x;
  static_assert(
      std::is_same_v<std::decay_t<decltype(Ax)>, StaticVector<double, 2>>,
      "Matrix * StaticVector should return StaticVector");
  if (!(std::abs(Ax[0] - 7.0) < 1e-12 && std::abs(Ax[1] - (-1.0)) < 1e-12)) {
    throw std::runtime_error("Matrix * StaticVector result mismatch");
  }

  // 3) solveLinearSystem with StaticVector RHS
  const StaticVector<double, 2> b{1.0, 2.0};
  const auto sol = solveLinearSystem(A, b);
  static_assert(
      std::is_same_v<std::decay_t<decltype(sol)>, StaticVector<double, 2>>,
      "solveLinearSystem(static RHS) should return StaticVector");

  const auto check = A * sol;
  for (Index i = 0; i < 2; ++i) {
    if (std::abs(check[i] - b[i]) > 1e-8)
      throw std::runtime_error("solveLinearSystem(static) verification failed");
  }

  // 4) eigen() now returns StaticVector eigenvalues
  const auto [eigvals, V] = A.eigen();
  static_assert(
      std::is_same_v<std::decay_t<decltype(eigvals)>, StaticVector<double, 2>>,
      "eigen() should return StaticVector eigenvalues");

  // Validate eigenpair residual for first vector.
  const auto v0 = V.getColVector(0);
  const auto Av0 = A * v0;
  const auto lv0 = eigvals[0] * v0;
  const auto r0 = Av0 - lv0;
  if (magnitude(r0) > 1e-6) {
    throw std::runtime_error("eigen() residual too large for first eigenpair");
  }

  std::cout << "Static-first Matrix feature tests passed!\n";
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////   MAIN   //////////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main() {
  std::cout << "=== MATRIX CLASS COMPREHENSIVE TEST SUITE ===" << std::endl;
  std::cout << "=============================================" << std::endl
            << std::endl;

  // Run all test cases
  testMatrixConstructors();
  testMatrixAccessors();
  testMatrixStaticMethods();
  testMatrixArithmetic();
  testMatrixMultiplication();
  testMatrixProperties();
  testDeterminantAndTrace();
  testMatrixInverse();
  testRowColumnOperations();
  testConcatenationAndSplit();
  testLinearSystem();
  testVectorIntegration();
  testTensorProduct();
  testTensorProductBlocks();
  testStaticFirstMatrixFeatures();

  std::cout << "\n=== ORIGINAL TESTS (from previous main) ===" << std::endl;

  // Your original tests
  Matrix<double, 3, 3> A{-3, 2, -1, 6, -6, 7, 3, -4, 4};
  Matrix<double, 3, 1> B{-1, -7, -6};

  std::cout << "Original problem: Solve Ax = B" << std::endl;
  std::cout << "A:\n" << A << std::endl;
  std::cout << "B:\n" << B << std::endl;

  try {
    Matrix<double, 3, 1> X{solveLinearSystem(A, B)};
    std::cout << "Solution X:\n" << X << std::endl;

    // Verify solution
    auto verification = A * X;
    std::cout << "Verification A*X:\n" << verification << std::endl;

    std::cout << "A inverse:\n" << A.inverse() << std::endl;

  } catch (const std::exception &e) {
    std::cout << "Error: " << e.what() << std::endl;
  }

  // Test various determinants
  std::cout << "\nDeterminant tests:" << std::endl;
  const Matrix<double, 1, 1> H1{6};
  const Matrix<double, 2, 2> H2{6, 3, 2, 4};
  const Matrix<double, 3, 3> D{2.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 3.0, 1.0};

  std::cout << "det(H1): " << det(H1) << std::endl;
  std::cout << "det(H2): " << det(H2) << std::endl;
  std::cout << "det(D): " << det(D) << std::endl;

  std::cout << "\n=== ALL TESTS COMPLETED SUCCESSFULLY! ===" << std::endl;

  std::cout << "=== VECTOR CLASS TEST SUITE ===" << std::endl;
  std::cout << std::endl;

  // Test 1: Constructors
  std::cout << "1. Testing Constructors:" << std::endl;
  DynamicVector<int> v1;                // Default constructor
  DynamicVector<int> v2(5);             // Size constructor (5 zeros)
  DynamicVector<int> v3{1, 2, 3, 4, 5}; // Initializer list

  std::cout << "v1 (default): " << v1 << " (size: " << v1.size() << ")"
            << std::endl;
  std::cout << "v2(5): " << v2 << " (size: " << v2.size() << ")" << std::endl;
  std::cout << "v3{1,2,3,4,5}: " << v3 << " (size: " << v3.size() << ")"
            << std::endl;
  std::cout << std::endl;

  // Test 2: Indexing
  std::cout << "2. Testing Indexing:" << std::endl;
  DynamicVector<double> v4{10.5, 20.3, 30.7, 40.1};
  std::cout << "v4: " << v4 << std::endl;
  std::cout << "v4[0] = " << v4[0] << std::endl;
  std::cout << "v4[2] = " << v4[2] << std::endl;

  // Modify element
  v4[1] = 99.9;
  std::cout << "After v4[1] = 99.9: " << v4 << std::endl;

  // Test signed indexing
  Index idx = 3;
  std::cout << "Using Index idx=3: v4[idx] = " << v4[idx] << std::endl;
  std::cout << std::endl;

  // Test 3: Vector Operations
  std::cout << "3. Testing Vector Operations:" << std::endl;
  DynamicVector<double> va{1.0, 2.0, 3.0};
  DynamicVector<double> vb{4.0, 5.0, 6.0};

  std::cout << "va = " << va << std::endl;
  std::cout << "vb = " << vb << std::endl;

  auto vc = va + vb;
  std::cout << "va + vb = " << vc << std::endl;

  auto vd = vb - va;
  std::cout << "vb - va = " << vd << std::endl;

  auto ve = 2.5 * va;
  std::cout << "2.5 * va = " << ve << std::endl;

  auto vf = va * 3.0;
  std::cout << "va * 3.0 = " << vf << std::endl;

  auto vg = -va;
  std::cout << "-va = " << vg << std::endl;
  std::cout << std::endl;

  // Test 4: Vector Math
  std::cout << "4. Testing Vector Math:" << std::endl;
  DynamicVector<double> u1{3.0, 4.0, 0.0};
  DynamicVector<double> u2{1.0, 0.0, 0.0};

  std::cout << "u1 = " << u1 << std::endl;
  std::cout << "u2 = " << u2 << std::endl;

  double dot = dotProduct(u1, u2);
  std::cout << "Dot product u1·u2 = " << dot << std::endl;

  auto cross = crossProduct(u1, u2);
  std::cout << "Cross product u1×u2 = " << cross << std::endl;

  double mag1 = magnitude(u1);
  double mag2 = magnitude(u2);
  std::cout << "Magnitude |u1| = " << mag1 << std::endl;
  std::cout << "Magnitude |u2| = " << mag2 << std::endl;

  auto unit1 = normalize(u1);
  auto unit2 = normalize(u2);
  std::cout << "Unit vector of u1 = " << unit1 << std::endl;
  std::cout << "Unit vector of u2 = " << unit2 << std::endl;
  std::cout << "Magnitude of unit u1 = " << magnitude(unit1) << std::endl;
  std::cout << std::endl;

  // Test 4a: Testing CrossProduct2 (Levi-Civita Method)
  std::cout << "4a. Testing CrossProduct2 (Levi-Civita Method):" << std::endl;
  auto cross2 = crossProduct2(u1, u2);
  std::cout << "Cross product2 u1 x u2 = " << cross2 << std::endl;
  std::cout << " crossProduct == crossProduct2 : "
            << (cross == cross2 ? "✓ SAME" : "✗ DIFFERENT") << std::endl;
  std::cout << std::endl;

  // Test 5: Resize and Push Back
  std::cout << "5. Testing Resize and Push Back:" << std::endl;
  DynamicVector<int> vr{10, 20, 30};
  std::cout << "Original vr: " << vr << " (size: " << vr.size() << ")"
            << std::endl;

  vr.resize(6);
  std::cout << "After resize(6): " << vr << " (size: " << vr.size() << ")"
            << std::endl;

  vr.resize(8, 77);
  std::cout << "After resize(8, 77): " << vr << " (size: " << vr.size() << ")"
            << std::endl;

  vr.push_back(100);
  std::cout << "After push_back(100): " << vr << " (size: " << vr.size() << ")"
            << std::endl;
  std::cout << std::endl;

  // Test 6: Edge Cases and Error Handling
  std::cout << "6. Testing Edge Cases:" << std::endl;

  try {
    DynamicVector<int> v_small{1, 2};
    DynamicVector<int> v_big{1, 2, 3, 4};
    std::cout << "Trying to add vectors of different sizes..." << std::endl;
    auto result = v_small + v_big;
  } catch (const std::exception &e) {
    std::cout << "✓ Caught expected exception: " << e.what() << std::endl;
  }

  try {
    DynamicVector<double> zero_vec{0.0, 0.0, 0.0};
    std::cout << "Trying to normalize zero vector..." << std::endl;
    auto unit = normalize(zero_vec);
  } catch (const std::exception &e) {
    std::cout << "✓ Caught expected exception: " << e.what() << std::endl;
  }

  try {
    DynamicVector<int> v2d_1{1, 2};
    DynamicVector<int> v2d_2{3, 4};
    std::cout << "Trying cross product on 2D vectors..." << std::endl;
    auto cross_2d = crossProduct(v2d_1, v2d_2);
  } catch (const std::exception &e) {
    std::cout << "✓ Caught expected exception: " << e.what() << std::endl;
  }

  std::cout << std::endl;
  std::cout << "=== ALL TESTS COMPLETED SUCCESSFULLY! ===" << std::endl;

  // Test 7: Projection
  std::cout << "7. Testing Vector Projection:" << std::endl;

  DynamicVector<double> proj_a{3.0, 4.0, 0.0}; // Vector a
  DynamicVector<double> proj_b{1.0, 0.0, 0.0}; // Unit vector along x-axis
  DynamicVector<double> proj_c{2.0, 2.0, 0.0}; // 45-degree vector

  std::cout << "Vector proj_a = " << proj_a << std::endl;
  std::cout << "Vector proj_b = " << proj_b << std::endl;
  std::cout << "Vector proj_c = " << proj_c << std::endl;

  // Vector projection
  auto vector_proj_ab = proj_a.projection(proj_b);
  auto vector_proj_ac = proj_a.projection(proj_c);
  std::cout << "Projection of proj_a onto proj_b: " << vector_proj_ab
            << std::endl;
  std::cout << "Projection of proj_a onto proj_c: " << vector_proj_ac
            << std::endl;

  std::cout << std::endl;

  // Test 8: Angle Calculation
  std::cout << "8. Testing Angle Calculation:" << std::endl;

  DynamicVector<double> angle_v1{1.0, 0.0, 0.0}; // Unit vector along x-axis
  DynamicVector<double> angle_v2{0.0, 1.0, 0.0}; // Unit vector along y-axis
  DynamicVector<double> angle_v3{1.0, 1.0, 0.0}; // 45-degree vector

  std::cout << "angle_v1 = " << angle_v1 << std::endl;
  std::cout << "angle_v2 = " << angle_v2 << std::endl;
  std::cout << "angle_v3 = " << angle_v3 << std::endl;

  // Test 90 degree angle
  double angle12_rad = angleRad(angle_v1, angle_v2);
  double angle12_deg = angleDegree(angle_v1, angle_v2);
  std::cout << "Angle between angle_v1 and angle_v2: " << angle12_rad
            << " rad = " << angle12_deg << " deg" << std::endl;

  // Test 45 degree angle
  double angle13_rad = angleRad(angle_v1, angle_v3);
  double angle13_deg = angleDegree(angle_v1, angle_v3);
  std::cout << "Angle between angle_v1 and angle_v3: " << angle13_rad
            << " rad = " << angle13_deg << " deg" << std::endl;

  std::cout << std::endl;

  // Test 9: Perpendicular and Parallel
  std::cout << "9. Testing Perpendicular and Parallel:" << std::endl;

  DynamicVector<double> perp_x{1.0, 0.0, 0.0}; // Unit vector along x
  DynamicVector<double> perp_y{0.0, 1.0, 0.0}; // Unit vector along y
  DynamicVector<double> para_x{2.0, 0.0, 0.0}; // Parallel to perp_x

  std::cout << "perp_x = " << perp_x << std::endl;
  std::cout << "perp_y = " << perp_y << std::endl;
  std::cout << "para_x = " << para_x << std::endl;

  // Test perpendicular (should be true)
  bool is_perp = isPerpendicular(perp_x, perp_y);
  std::cout << "isPerpendicular(perp_x, perp_y): " << is_perp << std::endl;

  // Test parallel (should be true)
  bool is_para = isParallel(perp_x, para_x);
  std::cout << "isParallel(perp_x, para_x): " << is_para << std::endl;

  // Show verification values
  std::cout << "Dot product perp_x·perp_y = " << dotProduct(perp_x, perp_y)
            << std::endl;
  std::cout << "Cross product magnitude |perp_x×para_x| = "
            << magnitude(crossProduct(perp_x, para_x)) << std::endl;

  std::cout << std::endl;

  Matrix<double, 2, 2> W{1, -2, -1, 3};
  Matrix<double, 2, 2> X{1, 1, 0, 1};
  std::cout << W * X * W.inverse();

  return 0;
}