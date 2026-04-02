#ifndef MY_MATRIX_CLASS_H
#define MY_MATRIX_CLASS_H

// Consider return m_elements[(i - 1) * nCols + (j - 1)]; in order to accessing
// the matrices with index starting from 1 in mathematic
// LU Decomposition (to be implemented) to be more optimized
// Resize function?
// Jacobian, Hessian Matrix

#include "Vector.h"       //My vector class
#include "comparison.h"   //Approximative Comparsion
#include "random.h"       //random
#include "signFunction.h" //sign()
#include <algorithm>      //max, min, swap, sort
#include <array>
#include <cassert>          //assert
#include <cmath>            //power
#include <functional>       // std::reference_wrapper
#include <initializer_list> //initiate list
#include <iomanip>          //tab
#include <iostream>
#include <span>        //Cheap view of array
#include <stdexcept>   //throw exception
#include <type_traits> // precision
#include <vector>

using Index = std::ptrdiff_t; // typedef

template <typename T, Index nRows, Index nCols> class Matrix {
private:
  // Check nRows, nCols > 0
  static_assert(nRows > 0 && nCols > 0 &&
                "Number of rows and columns must be greater than 0!");
  std::array<T, nRows * nCols> m_elements{};

public:
  // Constructors, Destructors
  Matrix(std::initializer_list<T> list) {
    assert(static_cast<Index>(list.size()) == nRows * nCols);
    std::copy(list.begin(), list.end(), m_elements.begin());
  }
  Matrix() = default;
  Matrix(const Matrix &) = default;
  Matrix(Matrix &&) = default;
  Matrix &operator=(const Matrix &) = default;
  Matrix &operator=(Matrix &&) = default;
  ~Matrix() = default;

  // 0 Matrix, 1 Matrix and Identity Matrix and functions to reset, static to
  // make them create one copy only
  static Matrix zero() { return Matrix{}; }
  void resetZero() { (*this) = Matrix::zero(); }

  static Matrix ones() {
    Matrix result{};
    std::fill(result.m_elements.begin(), result.m_elements.end(), 1);
    return result;
  }
  void resetOnes() { (*this) = Matrix::ones(); }

  static Matrix identity() {
    Matrix result{};
    for (Index i = 0; i < std::min(nRows, nCols); ++i)
      result(i, i) = T{1};
    return result;
  }
  void resetIdentity() { (*this) = Matrix::identity(); }

  // Accessing the elements in the array with one parameter (i)
  T &operator[](Index i) {
    assert(i >= 0 && i < this->length() && "Index out of bounds!");
    return m_elements[static_cast<std::size_t>(i)];
  }
  const T &operator[](Index i) const {
    assert(i >= 0 && i < this->length() && "Index out of bounds!");
    return m_elements[static_cast<std::size_t>(i)];
  }

  // Accessing the elements in the matrices with two parameters (i,j)
  T &operator()(Index i, Index j) {
    assert(i < nRows && j < nCols && i >= 0 && j >= 0);
    return m_elements[static_cast<std::size_t>(i * nCols + j)];
  }

  const T &operator()(Index i, Index j) const {
    assert(i < nRows && j < nCols && i >= 0 && j >= 0);
    return m_elements[static_cast<std::size_t>(i * nCols + j)];
  }

  // Getters to get number of rows & columns + total elements numbers
  constexpr Index getCols() const { return nCols; }
  constexpr Index getRows() const { return nRows; }
  constexpr Index length() const { return nRows * nCols; }

  // Getters for engineering tensor
  T &xx() {
    static_assert(nRows >= 1 && nCols >= 1, "xx() requires size >= 1x1");
    return (*this)(0, 0);
  }
  const T &xx() const {
    static_assert(nRows >= 1 && nCols >= 1, "xx() requires size >= 1x1");
    return (*this)(0, 0);
  }

  T &xy() {
    static_assert(nRows >= 1 && nCols >= 2, "xy() requires size >= 1x2");
    return (*this)(0, 1);
  }
  const T &xy() const {
    static_assert(nRows >= 1 && nCols >= 2, "xy() requires size >= 1x2");
    return (*this)(0, 1);
  }

  T &yx() {
    static_assert(nRows >= 2 && nCols >= 1, "yx() requires size >= 2x1");
    return (*this)(1, 0);
  }
  const T &yx() const {
    static_assert(nRows >= 2 && nCols >= 1, "yx() requires size >= 2x1");
    return (*this)(1, 0);
  }

  // Backward-compatible alias (typo kept): prefer yx().
  T &Yx() {
    static_assert(nRows >= 2 && nCols >= 1, "Yx() requires size >= 2x1");
    return yx();
  }
  const T &Yx() const {
    static_assert(nRows >= 2 && nCols >= 1, "Yx() requires size >= 2x1");
    return yx();
  }

  T &yy() {
    static_assert(nRows >= 2 && nCols >= 2, "yy() requires size >= 2x2");
    return (*this)(1, 1);
  }
  const T &yy() const {
    static_assert(nRows >= 2 && nCols >= 2, "yy() requires size >= 2x2");
    return (*this)(1, 1);
  }

  // 3D components (only available for matrices that are at least 3x3)
  T &xz() {
    static_assert(nRows >= 1 && nCols >= 3, "xz() requires size >= 1x3");
    return (*this)(0, 2);
  }
  const T &xz() const {
    static_assert(nRows >= 1 && nCols >= 3, "xz() requires size >= 1x3");
    return (*this)(0, 2);
  }

  T &yz() {
    static_assert(nRows >= 2 && nCols >= 3, "yz() requires size >= 2x3");
    return (*this)(1, 2);
  }
  const T &yz() const {
    static_assert(nRows >= 2 && nCols >= 3, "yz() requires size >= 2x3");
    return (*this)(1, 2);
  }

  T &zx() {
    static_assert(nRows >= 3 && nCols >= 1, "zx() requires size >= 3x1");
    return (*this)(2, 0);
  }
  const T &zx() const {
    static_assert(nRows >= 3 && nCols >= 1, "zx() requires size >= 3x1");
    return (*this)(2, 0);
  }

  T &zy() {
    static_assert(nRows >= 3 && nCols >= 2, "zy() requires size >= 3x2");
    return (*this)(2, 1);
  }
  const T &zy() const {
    static_assert(nRows >= 3 && nCols >= 2, "zy() requires size >= 3x2");
    return (*this)(2, 1);
  }

  T &zz() {
    static_assert(nRows >= 3 && nCols >= 3, "zz() requires size >= 3x3");
    return (*this)(2, 2);
  }
  const T &zz() const {
    static_assert(nRows >= 3 && nCols >= 3, "zz() requires size >= 3x3");
    return (*this)(2, 2);
  }

  // Reference to row(i)
  std::span<T> row(Index i) {
    assert(i < nRows && i >= 0);
    return std::span<T>(&m_elements[i * nCols], nCols);
  }

  std::span<const T> row(Index i) const {
    assert(i < nRows && i >= 0);
    return std::span<const T>(&m_elements[i * nCols], nCols);
  }

  // Reference to col(j): couldn't use std::span (data is not continuous in
  // m_element) and std::array here (as with array reference_wrapper default
  // constructor will fail)
  std::vector<std::reference_wrapper<T>> col(Index j) {
    assert(j < nCols && j >= 0);
    std::vector<std::reference_wrapper<T>> col_refs;
    col_refs.reserve(nRows);
    for (Index i = 0; i < nRows; ++i) {
      col_refs.push_back(std::ref((*this)(i, j)));
    }
    return col_refs;
  }

  // outputting matrix
  friend std::ostream &operator<<(std::ostream &out, const Matrix &matrix) {
    std::ios oldState(nullptr);
    oldState.copyfmt(out); // save stream state

    // Print floating point entries in scientific notation with 3 digits
    if constexpr (std::is_floating_point_v<T>) {
      out << std::fixed << std::setprecision(4);
      out << std::scientific;
    }

    constexpr int tab = 15;
    for (Index i = 0; i < nRows; ++i) {
      out << "|";
      for (Index j = 0; j < nCols; ++j) {
        if (approximatelyEqualAbsRel(matrix(i, j), T{0}))
          out << std::setw(tab) << T{0};
        else
          out << std::setw(tab) << matrix(i, j);
      }
      out << " |" << '\n';
    }

    out.copyfmt(oldState);
    return out;
  }

  // Matrix algebra operators
  friend Matrix operator*(T k, const Matrix &m) {
    Matrix<T, nRows, nCols> result{m};
    for (auto &e : result.m_elements) {
      e *= k;
    }
    return result;
  }
  friend Matrix operator*(const Matrix &m, T k) { return k * m; }

  friend Matrix operator+(const Matrix &m1, const Matrix &m2) {
    assert(m1.getCols() == m2.getCols() && m1.getRows() == m2.getRows() &&
           "Unable to perform matrix addition/substraction.");
    Matrix<T, nRows, nCols> result{m1};
    for (Index i{0}; i < result.length(); ++i)
      result[i] += m2[i];
    return result;
  }

  friend Matrix operator-(const Matrix &m1, const Matrix &m2) {
    return m1 + (-1) * m2;
  }

  Matrix operator-() const { return (-1) * (*this); }

  friend bool operator==(const Matrix &m1,
                         const Matrix &m2) // can be moved outside of the class
  {
    if (m1.getCols() != m2.getCols() || m1.getRows() != m2.getRows()) {
      return false;
    }
    for (Index i{0}; i < m1.length(); i++) {
      if (!approximatelyEqualAbsRel(m1[i], m2[i])) {
        return false;
      }
    }
    return true;
  }

  Matrix &operator+=(const Matrix &other) {
    assert(this->getCols() == other.getCols() &&
           this->getRows() == other.getRows() &&
           "Unable to perform matrix addition.");

    for (Index i = 0; i < this->length(); ++i) {
      (*this)[i] += other[i];
    }
    return *this;
  }

  Matrix &operator-=(const Matrix &other) {
    assert(this->getCols() == other.getCols() &&
           this->getRows() == other.getRows() &&
           "Unable to perform matrix subtraction.");

    for (Index i = 0; i < this->length(); ++i) {
      (*this)[i] -= other[i];
    }
    return *this;
  }

  Matrix &operator*=(T scalar) {
    for (auto &element : m_elements) {
      element *= scalar;
    }
    return *this;
  }

  template <Index R2, Index C2>
  Matrix<T, nRows, nCols> &operator*=(const Matrix<T, R2, C2> &m2) {
    static_assert(nRows == nCols, "Left matrix must be square.");
    static_assert(R2 == C2, "Right matrix must be square.");
    static_assert(nRows == R2, "Matrices must have the same size.");
    Matrix<T, nRows, C2> result = (*this) * m2;
    for (Index i = 0; i < nRows; ++i)
      for (Index j = 0; j < nCols; ++j)
        (*this)(i, j) = result(i, j);
    return *this;
  }

  friend bool operator!=(const Matrix &m1, const Matrix &m2) {
    return !(m1 == m2);
  }

  // Separate a matrix to two matrices by a column (position at the beginning
  // of the second matrix)
  template <Index leftCols, Index rightCols>
  void splitByColumn(Index colPos, Matrix<T, nRows, leftCols> &A,
                     Matrix<T, nRows, rightCols> &B) const {
    assert(colPos >= 0);
    assert(leftCols == colPos);
    assert(rightCols == nCols - colPos);

    for (Index i = 0; i < nRows; ++i)
      for (Index j = 0; j < nCols; ++j) {
        if (j < colPos)
          A(i, j) = (*this)(i, j);
        else
          B(i, j - colPos) = (*this)(i, j);
      }
  }

  // // Pointer to a whole row(i) or column(j)
  // std::array<T *, nCols> rowPointer(Index i)
  // {
  //   assert(i < nRows);
  //   std::array<T *, nCols> row_pointer{};
  //   for (Index j = 0; j < nCols; ++j)
  //   {
  //     row_pointer[j] = &(*this)(i, j);
  //   }
  //   return row_pointer;
  // }

  // std::array<T *, nRows> colPointer(Index j)
  // {
  //   assert(j < nCols);
  //   std::array<T *, nRows> col_pointer{};
  //   for (Index i = 0; i < nRows; ++i)
  //   {
  //     col_pointer[i] = &(*this)(i, j);
  //   }
  //   return col_pointer;
  // }

  // Find sub matrix to calculate minor of det(A)
  // At first I used 2 indices to loop through the Matrix, credit :
  // https://www.youtube.com/watch?v=YVk0nYrwBb0&t=1211s
  Matrix<T, (nRows - 1), (nCols - 1)> subMatrix(Index row, Index col) const {
    if (nRows <= 1 || nCols <= 1)
      throw std::invalid_argument(
          "Cannot create submatrix of size 0x0 from a 1x1 matrix.");
    Matrix<T, (nRows - 1), (nCols - 1)> sub{};
    Index k{0};
    for (Index i{0}; i < nRows; ++i) {
      for (Index j{0}; j < nCols; ++j) {
        if (i != row && j != col) {
          sub[k] = (*this)(i, j);
          ++k;
        }
      }
    }
    return sub;
  }

  // Swap 2 rows
  void swapRows(Index i1, Index i2) {
    for (Index j{0}; j < nCols; ++j) {
      std::swap((*this)(i1, j), (*this)(i2, j));
    }
  }

  // Find the row with the max element at a given column (Noted below pivot)
  Index indexRowMax(Index col, Index startRow = 0) const {
    assert(col < nCols && startRow < nRows && col >= 0 && "Index problem!");
    Index maxRow = startRow;
    T maxValue = std::abs((*this)(startRow, col));
    for (Index i = startRow + 1; i < nRows; ++i) {
      T value = std::abs((*this)(i, col));
      if (value > maxValue) {
        maxValue = value;
        maxRow = i;
      }
    }
    return maxRow;
  }

  // Transpose matrix
  Matrix<T, nCols, nRows> transpose() const {
    Matrix<T, nCols, nRows> result{};
    for (Index j = 0; j < nCols; ++j) {
      for (Index i = 0; i < nRows; ++i) {
        result(j, i) = (*this)(i, j);
      }
    }
    return result;
  }

  // Inverse Matrix by reduced row echelon form
  Matrix<T, nRows, nCols> inverse() const {
    if (!isSquare())
      throw std::invalid_argument("Inverse only defined for square matrices.");
    Matrix<T, nRows, nCols> inverseMatrix{};
    Matrix<T, nRows, nCols> testIdentity{};
    Matrix<T, nRows, (nCols + nCols)> augmentedMatrix{
        concatenateMatrixHorizontal(*this,
                                    Matrix<T, nRows, nCols>::identity())};
    // check det first
    if (approximatelyEqualAbsRel(det(*this), T{0})) {
      throw std::invalid_argument("Determinant of Matrix is 0, hence matrix is "
                                  "singular and cannot be inverted.");
    }
    for (Index pivotIndex{0}; pivotIndex < nRows; ++pivotIndex) {
      // std::cout << "Pivot index " << pivotIndex << std::endl;
      Index maxIndex{augmentedMatrix.indexRowMax(pivotIndex, pivotIndex)};
      if (maxIndex != pivotIndex) {
        augmentedMatrix.swapRows(pivotIndex, maxIndex);
        // std::cout << "Swap rows " << pivotIndex << " and " << maxIndex <<
        // std::endl; std::cout << augmentedMatrix << std::endl;
      }
      T pivot{augmentedMatrix(pivotIndex, pivotIndex)};
      if (approximatelyEqualAbsRel(pivot, T{})) {
        throw std::runtime_error("Matrix is singular and cannot be inverted.");
      }
      auto pivotSpan = augmentedMatrix.row(pivotIndex);
      // Normalize pivot in place
      if (!approximatelyEqualAbsRel(pivot, T{1})) {
        for (auto &val : pivotSpan)
          val /= pivot;
        // std::cout << "Normalize row " << pivotIndex << " by pivot = " <<
        // pivot << "\n" << augmentedMatrix << '\n';
      }

      // Eliminate other rows
      for (Index i{0}; i < nRows; ++i) {
        // Skip the pivot row
        if (i == pivotIndex)
          continue;

        T factor = augmentedMatrix(
            i, pivotIndex); // no need to /pivot as pivot = 1.0 already
        // Skip value under or above the pilot that are 0.0
        if (approximatelyEqualAbsRel(factor, T{}))
          continue;

        // std::cout << "Eliminate row " << i << " using pivot row " <<
        // pivotIndex
        //           << " (factor = " << factor << ")\n";
        auto targetSpan = augmentedMatrix.row(i);
        for (Index j = 0; j < std::ssize(pivotSpan); ++j) {
          targetSpan[j] -= factor * pivotSpan[j];
        }
        // std::cout << augmentedMatrix << '\n';
      }
    }
    augmentedMatrix.splitByColumn(nCols, testIdentity, inverseMatrix);
    if (testIdentity == Matrix::identity())
      return inverseMatrix;
    else
      throw std::runtime_error(
          "Matrix inversion failed: result is not an identity matrix.");
  }

  // Bool function to get the characteristics of matrix
  bool isDiagonal() const {
    for (Index i = 0; i < nRows; ++i) {
      for (Index j = 0; j < nCols; ++j) {
        if (i != j) {
          if (!approximatelyEqualAbsRel((*this)(i, j), T{0})) {
            return false;
          }
        }
      }
    }
    return true;
  }

  bool isUpperTriangular() const {
    for (Index i = 0; i < nRows; ++i) {
      for (Index j = 0; j < nCols; ++j) {
        if (i > j) {
          if (!approximatelyEqualAbsRel((*this)(i, j), T{0})) {
            return false;
          }
        }
      }
    }
    return true;
  }

  bool isLowerTriangular() const {
    for (Index i = 0; i < nRows; ++i) {
      for (Index j = 0; j < nCols; ++j) {
        if (i < j) {
          if (!approximatelyEqualAbsRel((*this)(i, j), T{0})) {
            return false;
          }
        }
      }
    }
    return true;
  }

  constexpr bool isSquare() const { return nRows == nCols; }

  bool isSymmetric() const { return ((*this) == this->transpose()); }

  bool isSkewSymmetric() const { return (-(*this) == (this->transpose())); }

  bool isOrthogonal() const {
    if (!this->isSquare())
      return false;
    return ((*this) * (this->transpose()) == Matrix::identity());
  }

  void reflect() // turn the upper trinagular matrix into symmetric matrix
  {
    // Ensure square matrix
    if (!this->isSquare())
      throw std::invalid_argument("reflect() requires a square matrix");

    // If matrix is upper triangular, copy upper -> lower to make symmetric
    if (this->isUpperTriangular()) {
      for (Index i = 0; i < nRows; ++i) {
        for (Index j = 0; j < i; ++j) {
          // set lower (i,j) = upper (j,i)
          (*this)(i, j) = (*this)(j, i);
        }
      }
      return;
    }
    // If matrix is lower triangular, copy lower -> upper
    if (this->isLowerTriangular()) {
      for (Index i = 0; i < nRows; ++i) {
        for (Index j = 0; j < i; ++j) {
          // set upper (j,i) = lower (i,j)
          (*this)(j, i) = (*this)(i, j);
        }
      }
      return;
    }
    return;
  }

  // Convert a Matrix type to a Vector type by using Constructor
  Matrix(const DynamicVector<T> &vec) {
    static_assert(nCols == 1, "Can only construct column matrix from Vector");
    if (vec.size() != nRows)
      throw std::invalid_argument("Vector size must match matrix rows");
    for (Index i = 0; i < nRows; ++i)
      (*this)(i, 0) = vec[i];
  }

  Matrix(const StaticVector<T, nRows> &vec) {
    static_assert(nCols == 1,
                  "Can only construct column matrix from Vector (static)");
    for (Index i = 0; i < nRows; ++i)
      (*this)(i, 0) = vec[i];
  }

  // Type conversion from Matrix to Vector
  operator DynamicVector<T>() const {
    static_assert(nCols == 1, "Can only convert column matrix to Vector");
    DynamicVector<T> result;
    result.resize(nRows);
    for (Index i = 0; i < nRows; ++i) {
      result[i] = (*this)(i, 0);
    }
    return result;
  }

  operator StaticVector<T, nRows>() const {
    static_assert(nCols == 1, "Can only convert column matrix to Vector");
    StaticVector<T, nRows> result;
    for (Index i = 0; i < nRows; ++i) {
      result[i] = (*this)(i, 0);
    }
    return result;
  }

  StaticVector<T, nRows> getColVector(Index j) const {
    StaticVector<T, nRows> v{};
    for (Index i = 0; i < nRows; ++i)
      v[i] = (*this)(i, j);
    return v;
  }

  StaticVector<Matrix<T, nRows, nCols>, 2> symmetricDecomposition() const {
    return {0.5 * ((*this) + this->transpose()),
            0.5 * ((*this) - this->transpose())};
  }

  // StaticVector<T, nRows> getColVector(Index j) const {
  //   StaticVector<T, nRows> v{};
  //   for (Index i = 0; i < nRows; ++i)
  //     v[i] = (*this)(i, j);
  //   return v;
  // }

  // How about static Vector
  // Real Eigenvalues and their Eigen Vectors for symmetric matrices:: AX = kX
  std::pair<StaticVector<T, nRows>, Matrix>
  eigen(Index maxIterations = 1000) const {
    const bool symmetric = this->isSymmetric();

    // Phase 1: QR iteration to estimate eigenvalues.
    Index iteration{0};
    Matrix A_qr = *this;
    Matrix V = symmetric ? Matrix::identity() : Matrix::zero();

    while (iteration < maxIterations && !A_qr.isUpperTriangular()) {
      auto [Q, R] = QRdecomposition(A_qr);
      A_qr = R * Q;
      if (symmetric) {
        V = V * Q; // for symmetric: columns converge to eigenvectors
      }
      ++iteration;
    }

    if (iteration == maxIterations) {
      std::cout << "Max iteration reached without achieving upper triangular "
                   "form for A!"
                << '\n';
    }

    StaticVector<T, nRows> eigenValues{};
    for (Index i = 0; i < nRows; ++i)
      eigenValues[i] = A_qr(i, i);

    // Phase 2 (non-symmetric): compute right-eigenvectors using shifted inverse
    // iteration for each (real) eigenvalue estimate.
    if (!symmetric) {
      const Matrix A0 = *this;
      const Matrix I = Matrix::identity();
      const T tol = static_cast<T>(1e-10);
      const T epsBase = static_cast<T>(1e-8);

      for (Index j = 0; j < nRows; ++j) {
        const T mu0 = eigenValues[j];

        StaticVector<T, nRows> x{};
        for (Index i = 0; i < nRows; ++i)
          x[i] = static_cast<T>(Random::get(-10, 10));

        // Avoid zero initial vector
        if (approximatelyEqualAbsRel(magnitude(x), T{0}))
          x[0] = T{1};
        x = x / magnitude(x);

        StaticVector<T, nRows> xPrev = StaticVector<T, nRows>::zero();

        for (Index invIter = 0; invIter < maxIterations; ++invIter) {
          if (magnitude(x - xPrev) < tol)
            break;
          xPrev = x;

          StaticVector<T, nRows> y{};
          bool solved = false;

          for (int attempt = 0; attempt < 4 && !solved; ++attempt) {
            const T mu = mu0 + epsBase * static_cast<T>(attempt);
            const Matrix shifted = A0 - mu * I;
            try {
              y = solveLinearSystem(shifted, x);
              solved = true;
            } catch (const std::exception &) {
              // try a slightly perturbed shift
            }
          }
          if (!solved)
            break;
          const T yNorm = magnitude(y);
          if (approximatelyEqualAbsRel(yNorm, T{0}))
            break;
          x = y / yNorm;
        }
        for (Index i = 0; i < nRows; ++i)
          V(i, j) = x[i];
      }
    }

    // Normalize eigenvectors (each column of V) to unit length.
    for (Index j = 0; j < nRows; ++j) {
      T norm2{};
      for (Index i = 0; i < nRows; ++i)
        norm2 += V(i, j) * V(i, j);
      const T norm = std::sqrt(norm2);
      if (approximatelyEqualAbsRel(norm, T{0}))
        continue;
      for (Index i = 0; i < nRows; ++i)
        V(i, j) /= norm;
    }

    // V^T is the tensor which transform current xyz coordinate system to the
    // the system defined by normalized eigenvectors of the current tensor
    // getColVector to retrieve the eigenvectors[i] corresponding to
    // eigenvalues[i]
    return {eigenValues, V};
  }

  // turn matrix from (x,y,z) basis to eigenvectors basis (diagonalization)
  Matrix toEigenBasis() const {
    auto [_, V] = this->eigen();
    const Matrix A{*this};

    // Symmetric case: eigenvectors are (approximately) orthonormal
    if (this->isSymmetric()) {
      return (V.transpose() * A * V);
    }
    const Matrix Vinv = V.inverse();
    return (Vinv * A * V);
  }

  Matrix symPower(T pow) const {
    assert(this->isSymmetric() && "Only implemented for symmetric matrix!");
    Matrix result{};
    auto [eigvals, V] = this->eigen();
    for (Index i{0}; i < eigvals.size(); ++i) {
      const StaticVector<T, nRows> v = V.getColVector(i);
      result += std::pow(eigvals[i], pow) * tensorProduct<nRows, nRows>(v, v);
    }
    return result;
  }

  T firstInvariant() const { return trace(*this); }
  T secondInvariant() const {
    return T{0.5} * (trace(*this) * trace(*this) - trace((*this) * (*this)));
  }
  T thirdInvariant() const { return det(*this); }

  StaticVector<T, 3> invariants() const {
    return {this->firstInvariant(), this->secondInvariant(),
            this->thirdInvariant()};
  }
};

/////////////////////////// END OF MATRIX CLASS//////////////////////////
// Householder reflection:   A  = Q * R
template <typename T, Index nRows, Index nCols>
std::pair<Matrix<T, nRows, nCols>, Matrix<T, nRows, nCols>>
QRdecomposition(const Matrix<T, nRows, nCols> &A) {
  assert(A.isSquare() && "Only implemented for square matrices");
  Matrix<T, nRows, nCols> R{A};
  Matrix<T, nRows, nCols> Q{Matrix<T, nRows, nCols>::identity()};

  for (Index j{0}; j < nCols - 1; ++j) // Except the last column
  {
    const Index subSize{nCols - j};
    DynamicVector<T> a1(subSize);
    DynamicVector<T> b1(subSize);
    b1[0] = T{1.0}; // basis vector
    for (Index i = j; i < nRows; ++i) {
      a1[i - j] = R(i, j);
    }
    DynamicVector<T> n =
        normalize(a1 - (-sgn(a1[0])) * magnitude(a1) * b1); // n = û

    // Build full Householder matrix P = I - 2 n n^T (embedded in nRows x
    // nCols)
    Matrix<T, nRows, nCols> P = Matrix<T, nRows, nCols>::identity();
    for (Index row = j; row < nRows; ++row) {
      for (Index col = j; col < nCols; ++col) {
        P(row, col) -= T{2.0} * n[row - j] * n[col - j];
      }
    }
    R = P * R;
    Q = Q * P;
  }
  return {Q, R};
}

// Concatenate 2 matrices horizontally (into augemented matrix)
//  (A) + (B) = ( A B )
template <typename T, Index R1, Index C1, Index R2, Index C2>
Matrix<T, R1, (C1 + C2)>
concatenateMatrixHorizontal(const Matrix<T, R1, C1> &A,
                            const Matrix<T, R2, C2> &B) {
  assert(R1 == R2 && "The number of rows of both matrices must match to create "
                     "an concatenated matrix.\n");
  Matrix<T, R1, (C1 + C2)> result{};
  for (Index i = 0; i < R1; ++i)
    for (Index j = 0; j < (C1 + C2); ++j) {
      if (j < C1) {
        assert(j < A.getCols() && j >= 0 &&
               "Error index exceeding matrix's size.\n");
        result(i, j) = A(i, j);
      } else {
        assert((j - C1) < B.getCols() && j >= 0 &&
               "Error index exceeding matrix's size.\n");
        result(i, j) = B(i, j - C1);
      }
    }
  return result;
}

// Concatenate 2 matrices vertically
//   |A| & |B| -> | A |
//                | B |
template <typename T, Index R1, Index C1, Index R2, Index C2>
Matrix<T, (R1 + R2), C1> concatenateMatrixVertical(const Matrix<T, R1, C1> &A,
                                                   const Matrix<T, R2, C2> &B) {
  assert(C1 == C2 && "The number of columns of both matrices must match for "
                     "vertical concatenation.");
  Matrix<T, (R1 + R2), C1> result{};
  for (Index i = 0; i < (R1 + R2); ++i)
    for (Index j = 0; j < C1; ++j) {
      if (i < R1) {
        assert(i < A.getRows() && i >= 0 &&
               "Error index exceeding matrix's size.\n");
        result(i, j) = A(i, j);
      } else {
        assert((i - R1) < B.getRows() && i >= 0 &&
               "Error index exceeding matrix's size.\n");
        result(i, j) = B(i - R1, j);
      }
    }
  return result;
}

// Solve Ax = B by taking x = A.inverse()*B
template <typename T, Index R1, Index C1>
Matrix<T, R1, 1> solveLinearSystem(const Matrix<T, R1, C1> &A,
                                   const Matrix<T, R1, 1> &B) {
  Matrix<T, R1, C1> inverseA{A.inverse()};
  return inverseA * B;
}

// Solve Ax = B by taking x = A.inverse()*B
template <typename T, Index R1, Index C1>
DynamicVector<T> solveLinearSystem(const Matrix<T, R1, C1> &A,
                                   const DynamicVector<T> &B) {
  Matrix<T, R1, 1> Bmatrix{B};
  Matrix<T, R1, C1> inverseA{A.inverse()};
  Matrix<T, R1, 1> result{inverseA * Bmatrix};
  return DynamicVector<T>(result);
}

template <typename T, Index R1, Index C1>
StaticVector<T, R1> solveLinearSystem(const Matrix<T, R1, C1> &A,
                                      const StaticVector<T, R1> &B) {
  Matrix<T, R1, 1> Bmatrix{};
  for (Index i = 0; i < R1; ++i)
    Bmatrix(i, 0) = B[i];

  Matrix<T, R1, C1> inverseA{A.inverse()};
  Matrix<T, R1, 1> result{inverseA * Bmatrix};
  return static_cast<StaticVector<T, R1>>(result);
}

// Matrix multiplication
template <typename T, Index R1, Index C1, Index R2, Index C2>
Matrix<T, R1, C2> operator*(const Matrix<T, R1, C1> &m1,
                            const Matrix<T, R2, C2> &m2) {
  assert(C1 == R2 && "Not suitable for matrix multiplication.");
  Matrix<T, R1, C2> result{};
  for (Index i = 0; i < R1; ++i)
    for (Index j = 0; j < C2; ++j)
      for (Index k = 0; k < C1; ++k)
        result(i, j) += m1(i, k) * m2(k, j);
  return result;
}

// Matrix-vector multiplication (m * v)
template <typename T, Index R, Index C>
DynamicVector<T> operator*(const Matrix<T, R, C> &m,
                           const DynamicVector<T> &v) {
  assert(v.size() == C && "Vector size must match matrix column count.");
  const Matrix<T, C, 1> col{v};
  const Matrix<T, R, 1> result = m * col;
  return static_cast<DynamicVector<T>>(
      result); // uses Matrix::operator DynamicVector<T>()
}

template <typename T, Index R, Index C>
StaticVector<T, R> operator*(const Matrix<T, R, C> &m,
                             const StaticVector<T, C> &v) {
  const Matrix<T, C, 1> col{v};
  const Matrix<T, R, 1> result = m * col;
  return static_cast<StaticVector<T, R>>(result);
}

// Orthogonal bool function
template <typename T, Index R1, Index C1, Index R2, Index C2>
bool arePairOrthogonal(const Matrix<T, R1, C1> &m1,
                       const Matrix<T, R2, C2> &m2) {
  if ((m1 != Matrix<T, R1, C1>::zero()) && (m2 != Matrix<T, R2, C2>::zero())) {
    return ((m1.transpose() * m2) == Matrix<T, R1, C2>::identity());
  }
  return false; // Both matrices must not be 0
}

// find the most 0 row to calculate det(A) faster
template <typename T, Index nRows, Index nCols>
Index mostZeroRow(const Matrix<T, nRows, nCols> &m) {
  std::array<int, nRows> exceed0{};
  Index index{0};
  for (Index i{0}; i < nRows; i++) {
    for (const auto &element : m.row(i)) {
      if (approximatelyEqualAbsRel(element, T{}))
        ++exceed0[i];
    }
  }
  for (Index i{1}; i < nRows; i++) {
    if (exceed0[i] > exceed0[index]) {
      index = i;
    }
  }
  return index;
}

// Calculate trace(A)
template <typename T, Index nRows, Index nCols>
T trace(const Matrix<T, nRows, nCols> &m) {
  static_assert(nRows == nCols, "Trace applies only to squared matrix");
  T result{};
  for (Index i = 0; i < nRows; ++i)
    result += m(i, i);
  return result;
}

// Calculate det(A) using Laplace expansion
// Though this algorithm is not efficient with huge matrices. Considering
// implementing LU decomposition insted... template <typename T> T det(const
// Matrix<T, 1, 1> &m)
// {
//   return m(0, 0);
// }

// // Base case: 2x2
// template <typename T>
// T det(const Matrix<T, 2, 2> &m)
// {
//   return m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0);
// }

template <typename T, Index nRows, Index nCols>
T det(const Matrix<T, nRows, nCols> &m) {
  static_assert(nRows == nCols,
                "Determinant can only be computed for square matrices");

  if constexpr (nRows == 1)
    return m(0, 0);
  else if constexpr (nRows == 2)
    return m(0, 0) * m(1, 1) - m(0, 1) * m(1, 0);
  else if constexpr (nRows == 3) // not necessary, but optimized for 3x3 matrix
                                 // (which is the most often case in my need)
    return m(0, 0) * (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)) -
           m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0)) +
           m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
  else {
    T detSum{};
    Index rowNum{mostZeroRow(m)};
    for (Index j = 0; j < nCols; ++j) {
      T a = m(rowNum, j);
      if (approximatelyEqualAbsRel(a, T{}))
        continue;
      int sign = ((rowNum + j) % 2 == 0) ? 1 : -1; // (-1)^(row+col)
      Matrix<T, nRows - 1, nCols - 1> sub = m.subMatrix(rowNum, j);
      detSum += static_cast<T>(sign) * a * det(sub);
    }
    return detSum;
  }
}

// Tensor Product function with explicit template parameters for vectors
template <typename T, Index R, Index C>
Matrix<T, R, C> tensorProduct(const DynamicVector<T> &v1,
                              const DynamicVector<T> &v2) {
  // Runtime validation
  if (v1.size() != R || v2.size() != C) {
    throw std::invalid_argument("Vector sizes must match template parameters");
  }

  Matrix<T, R, C> result{};
  for (Index i = 0; i < R; ++i) {
    for (Index j = 0; j < C; ++j) {
      result(i, j) = v1[i] * v2[j];
    }
  }
  return result;
}

// Backward-compatible overload for calls like tensorProduct<2, 2>(v1, v2)
template <Index R, Index C, typename T>
Matrix<T, R, C> tensorProduct(const DynamicVector<T> &v1,
                              const DynamicVector<T> &v2) {
  return tensorProduct<T, R, C>(v1, v2);
}

template <typename T, Index R, Index C>
Matrix<T, R, C> tensorProduct(const StaticVector<T, R> &v1,
                              const StaticVector<T, C> &v2) {
  Matrix<T, R, C> result{};
  for (Index i = 0; i < R; ++i) {
    for (Index j = 0; j < C; ++j) {
      result(i, j) = v1[i] * v2[j];
    }
  }
  return result;
}

template <Index R, Index C, typename T>
Matrix<T, R, C> tensorProduct(const StaticVector<T, R> &v1,
                              const StaticVector<T, C> &v2) {
  return tensorProduct<T, R, C>(v1, v2);
}

template <typename T, Index R1, Index C1, Index R2, Index C2>
Matrix<T, R1 * R2, C1 * C2> tensorProduct(const Matrix<T, R1, C1> &A,
                                          const Matrix<T, R2, C2> &B) {
  Matrix<T, R1 * R2, C1 * C2> result{};

  // Each element A(i,j) multiplies the entire matrix B
  for (Index i = 0; i < R1; ++i) {
    for (Index j = 0; j < C1; ++j) {
      // Fill block at position (i*R2, j*C2)
      for (Index k = 0; k < R2; ++k) {
        for (Index l = 0; l < C2; ++l) {
          result(i * R2 + k, j * C2 + l) = A(i, j) * B(k, l);
        }
      }
    }
  }
  return result;
}

template <typename T, Index R1, Index C1, Index R2, Index C2>
T tensorContraction(const Matrix<T, R1, C1> &A, const Matrix<T, R2, C2> &B) {
  return trace(A.transpose() * B);
}

template <typename T, Index R, Index C> T tensorNorm(const Matrix<T, R, C> &A) {
  return constexpr_sqrt(tensorContraction(A, A));
}

// Hadamard product (element-wise multiplication) for matrices
template <typename T, Index nRows, Index nCols>
Matrix<T, nRows, nCols> hadamardProduct(const Matrix<T, nRows, nCols> &A,
                                        const Matrix<T, nRows, nCols> &B) {
  Matrix<T, nRows, nCols> result{};
  for (Index i = 0; i < nRows * nCols; ++i) {
    result[i] = A[i] * B[i];
  }
  return result;
}

#endif