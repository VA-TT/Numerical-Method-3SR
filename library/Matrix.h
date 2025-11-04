#ifndef MY_MATRIX_CLASS
#define MY_MATRIX_CLASS

// Consider return m_elements[(i - 1) * nCols + (j - 1)]; in order to accessing
// the matrices with index starting from 1 in mathematic
// LU Decomposition (to be implemented) to be more optimized
// Resize function?
// EigenValue, EigenVector? Dyadic? Power Matrix?

#include "Vector.h"     //My vector class
#include "comparison.h" //Approximative Comparsion
#include <algorithm>    //max, min, swap, sort
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
    assert(static_cast<int>(list.size()) == nRows * nCols);
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
      result(i, i) = 1.0;
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
      out << std::scientific << std::setprecision(1);
    }

    constexpr int tab = 15;
    for (Index i = 0; i < nRows; ++i) {
      out << "|";
      for (Index j = 0; j < nCols; ++j) {
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

  friend bool operator!=(const Matrix &m1, const Matrix &m2) {
    return !(m1 == m2);
  }

  // Separate a matrix to two matrices by a column (position at the beginning of
  // the second matrix)
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
    if (approximatelyEqualAbsRel(det(*this), 0.0)) {
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
      if (approximatelyEqualAbsRel(pivot, 0.0)) {
        throw std::runtime_error("Matrix is singular and cannot be inverted.");
      }
      auto pivotSpan = augmentedMatrix.row(pivotIndex);
      // Normalize pivot in place
      if (!approximatelyEqualAbsRel(pivot, 1.0)) {
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
        if (approximatelyEqualAbsRel(
                factor,
                0.0)) // Skip value under or above the pilot that are 0.0
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
          if ((*this)(i, j) != 0) {
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
          if ((*this)(i, j) != 0) {
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
          if ((*this)(i, j) != 0) {
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

  // Convert a Matrix type to a Vector type
  Matrix(const Vector<T> &vec) {
    static_assert(nCols == 1, "Can only construct column matrix from Vector");
    if (vec.size() != nRows)
      throw std::invalid_argument("Vector size must match matrix rows");
    for (Index i = 0; i < nRows; ++i)
      (*this)(i, 0) = vec[i];
  }

  // Type conversion from Matrix to Vector
  operator Vector<T>() const {
    static_assert(nCols == 1, "Can only convert column matrix to Vector");
    Vector<T> result;
    result.resize(nRows);
    for (Index i = 0; i < nRows; ++i) {
      result[i] = (*this)(i, 0);
    }
    return result;
  }
};

/////////////////////////// END OF MATRIX CLASS//////////////////////////

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
//   |A|  = | A |
// + |B|    | B |
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
  // try {
  //   Matrix<T, R1, C1> inverseA{A.inverse()};
  //   return inverseA * B;
  // } catch (const std::exception &) {
  //   throw std::invalid_argument("A is not invertible.");
  // }
  Matrix<T, R1, C1> inverseA{A.inverse()};
  return inverseA * B;
}

// Solve Ax = B by taking x = A.inverse()*B
template <typename T, Index R1, Index C1>
Vector<T> solveLinearSystem(const Matrix<T, R1, C1> &A, const Vector<T> &B) {
  // if (B.size() != R1) {
  //   throw std::invalid_argument("Vector B size must match matrix A rows");
  // }
  // try {
  //   Matrix<T, R1, 1> Bmatrix{B};
  //   Matrix<T, R1, C1> inverseA{A.inverse()};
  //   Matrix<T, R1, 1> result{inverseA * Bmatrix};
  //   return Vector<T>(result);
  // } catch (const std::exception &) {
  //   throw std::invalid_argument("A is not invertible.");
  // }
  Matrix<T, R1, 1> Bmatrix{B};
  Matrix<T, R1, C1> inverseA{A.inverse()};
  Matrix<T, R1, 1> result{inverseA * Bmatrix};
  return Vector<T>(result);
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
      if (approximatelyEqualAbsRel(element, 0.0))
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
  else {
    T detSum{};
    Index rowNum{mostZeroRow(m)};
    for (Index j = 0; j < nCols; ++j) {
      T a = m(rowNum, j);
      if (approximatelyEqualAbsRel(a, 0.0))
        continue;
      int sign = ((rowNum + j) % 2 == 0) ? 1 : -1; // (-1)^(row+col)
      Matrix<T, nRows - 1, nCols - 1> sub = m.subMatrix(rowNum, j);
      detSum += static_cast<T>(sign) * a * det(sub);
    }
    return detSum;
  }
}

// Tensor Product function with explicit template parameters for vectors
template <Index R, Index C, typename T>
Matrix<T, R, C> tensorProduct(const Vector<T> &v1, const Vector<T> &v2) {
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

#endif