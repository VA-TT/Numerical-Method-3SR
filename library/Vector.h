#ifndef MY_VECTOR_CLASS_H
#define MY_VECTOR_CLASS_H

#include "Matrix.h"
#include "comparison.h" //Approximative Comparsion
#include "kroneckerDelta_LeviCivita.h"
#include "physicConstants.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using Index = std::ptrdiff_t; // typedef

template <typename T, Index n>
class StaticVector; // class' declaration, definition below

///////////////////////////////////////////////////
/////////////// DYNAMIC CLASS VECTOR //////////////
///////////////////////////////////////////////////

// Consider integrated to Class (dynamic+static) into 1 by using is_static flag
// constexpr Index dynamic = -1;
// template <typename T, Index n = dynamic> class DynamicVector {}
// static constexpr bool is_static = (n != Dynamic);
// if constexpr (DynamicVector<T, n>::is_static) {
//   // vector static
// } else {
//   // vector dynamic
// }
template <typename T> class DynamicVector {
private:
  std::vector<T> m_elements{};

public:
  // Default Constructors
  DynamicVector() = default;
  DynamicVector(const DynamicVector &) = default;
  DynamicVector(DynamicVector &&) = default;
  DynamicVector &operator=(const DynamicVector &) = default;
  DynamicVector &operator=(DynamicVector &&) = default;
  ~DynamicVector() = default;

  // Constructor with length n of elements which are 0
  explicit DynamicVector(std::size_t n) : m_elements(n, T{}) {}

  // // Constructor with length n of elements which are value
  DynamicVector(std::size_t n, const T &value) : m_elements(n, value) {}

  // Constructor with initializer_list
  DynamicVector(std::initializer_list<T> list) : m_elements(list) {}

  template <Index n> explicit DynamicVector(const StaticVector<T, n> &v) {
    m_elements.resize(static_cast<std::size_t>(n));
    for (Index i = 0; i < n; ++i)
      m_elements[static_cast<std::size_t>(i)] = v[i];
  }

  // DynamicVector 0
  static DynamicVector<T> zero(Index n) { return DynamicVector<T>(n); }

  // Reset all value to 0 default
  void resetZero() {
    for (auto &e : m_elements) {
      if constexpr (requires(T &x) { x.resetZero(); }) {
        e.resetZero();
      } else {
        e = T{};
      }
    }
  }

  // Signed indexing
  auto &operator[](Index i) {
    assert(i >= 0 && "Negative index not allowed");
    assert(static_cast<std::size_t>(i) < m_elements.size() &&
           "Index out of bounds");
    return m_elements.data()[static_cast<std::size_t>(i)];
  }

  const auto &operator[](Index i) const {
    assert(i >= 0 && "Negative index not allowed");
    assert(static_cast<std::size_t>(i) < m_elements.size() &&
           "Index out of bounds");
    return m_elements.data()[static_cast<std::size_t>(i)];
  }

  // Size functions
  // constexpr std::size_t size() const noexcept { return m_elements.size(); }
  constexpr Index size() const noexcept {
    return static_cast<Index>(m_elements.size());
  }

  // Convenience coordinate accessors for 2D/3D Vectors
  T &x() {
    assert(size() >= 1 && "x() requires size() >= 1");
    return (*this)[0];
  }
  const T &x() const {
    assert(size() >= 1 && "x() requires size() >= 1");
    return (*this)[0];
  }

  T &y() {
    assert(size() >= 2 && "y() requires size() >= 2");
    return (*this)[1];
  }
  const T &y() const {
    assert(size() >= 2 && "y() requires size() >= 2");
    return (*this)[1];
  }

  T &z() {
    assert(size() >= 3 && "z() requires size() >= 3");
    return (*this)[2];
  }
  const T &z() const {
    assert(size() >= 3 && "z() requires size() >= 3");
    return (*this)[2];
  }

  T &at(Index i) {
    if (i < 0)
      throw std::out_of_range("Negative index not allowed");
    return m_elements.at(static_cast<std::size_t>(i));
  }

  const T &at(Index i) const {
    if (i < 0)
      throw std::out_of_range("Negative index not allowed");
    return m_elements.at(static_cast<std::size_t>(i));
  }

  // Iterators
  auto begin() { return m_elements.begin(); }
  auto end() { return m_elements.end(); }
  auto begin() const { return m_elements.begin(); }
  auto end() const { return m_elements.end(); }

  // Resize + Reserve vector
  void resize(Index n) {
    assert(n >= 0 && "Negative size not allowed");
    m_elements.resize(static_cast<std::size_t>(n));
  }
  void resize(Index n, const T &value) {
    assert(n >= 0 && "Negative size not allowed");
    m_elements.resize(static_cast<std::size_t>(n), value);
  }
  void reserve(Index n) {
    assert(n >= 0 && "Negative reserve not allowed");
    m_elements.reserve(static_cast<std::size_t>(n));
  }

  // Add element
  void push_back(const T &value) { m_elements.push_back(value); }

  // Access last element
  T &back() {
    assert(!m_elements.empty() && "Cannot call back() on empty vector");
    return m_elements.back();
  }

  const T &back() const {
    assert(!m_elements.empty() && "Cannot call back() on empty vector");
    return m_elements.back();
  }

  // Print
  friend std::ostream &operator<<(std::ostream &os, const DynamicVector<T> &v) {
    os << "[";
    for (Index i = 0; i < v.size(); ++i) {
      if constexpr (requires(const T &x) {
                      approximatelyEqualAbsRel(x, T{0});
                    }) {
        if (approximatelyEqualAbsRel(v[i], T{0}))
          os << T{0};
        else
          os << v[i];
      } else {
        os << v[i];
      }
      if (i < v.size() - 1)
        os << ", ";
    }
    os << "]";
    return os;
  }

  //+= operator
  DynamicVector &operator+=(const DynamicVector<T> &other) {
    if (this->size() != other.size())
      throw std::invalid_argument(
          "Vectors must have the same dimension for operator+=.");
    for (Index i{0}; i < this->size(); ++i)
      (*this)[i] += other[i];
    return *this;
  }

  DynamicVector &operator-=(const DynamicVector<T> &other) {
    return *this += -other;
  }

  // Projection on another vector
  DynamicVector projection(const DynamicVector &other) const {
    DynamicVector unit{normalize(other)};
    return unit * dotProduct(*this, unit);
  }

  // Sort in-place. Default: descending order.
  void sort(bool descending = true) {
    if (descending) {
      std::sort(this->begin(), this->end(),
                [](const T &a, const T &b) { return a > b; });
    } else {
      std::sort(this->begin(), this->end(),
                [](const T &a, const T &b) { return a < b; });
    }
  }

  // Return a sorted copy. Default: descending order.
  DynamicVector sorted(bool descending = true) const {
    DynamicVector copy = *this;
    copy.sort(descending);
    return copy;
  }

  template <Index n> DynamicVector &operator=(const StaticVector<T, n> &v) {
    m_elements.resize(static_cast<std::size_t>(n));
    for (Index i = 0; i < n; ++i)
      m_elements[i] = v[i];
    return *this;
  }
};

////////////////// END OF DYNAMIC VECTOR CLASS /////////////////////////////

// DynamicVector operators
template <typename T>
bool operator==(const DynamicVector<T> &v1, const DynamicVector<T> &v2) {
  if (v1.size() != v2.size())
    return false;
  for (Index i = 0; i < v1.size(); ++i) {
    if (!approximatelyEqualAbsRel(v1[i], v2[i]))
      return false;
  }
  return true;
}

template <typename T>
bool operator!=(const DynamicVector<T> &v1, const DynamicVector<T> &v2) {
  return !(v1 == v2);
}

template <typename T>
DynamicVector<T> operator+(const DynamicVector<T> &v1,
                           const DynamicVector<T> &v2) {
  if (v1.size() != v2.size())
    throw std::invalid_argument("Vectors must have the same dimension.");

  DynamicVector<T> result(v1.size());
  for (Index i = 0; i < v1.size(); ++i)
    result[i] = v1[i] + v2[i];
  return result;
}

// Scalar multiplication (scalar * vector)
template <typename T>
DynamicVector<T> operator*(const T &k, const DynamicVector<T> &v) {
  DynamicVector<T> result(v.size());
  for (Index i = 0; i < v.size(); ++i)
    result[i] = k * v[i];
  return result;
}

// Scalar multiplication (vector * scalar)
template <typename T>
DynamicVector<T> operator*(const DynamicVector<T> &v, const T &k) {
  return k * v;
}

// Element-wise multiplication
template <typename T>
DynamicVector<T> operator*(const DynamicVector<T> &v1,
                           const DynamicVector<T> &v2) {
  assert(v1.size() == v2.size() && "Can't perfom element-wised mulplication "
                                   "on 2 different-sized vectors!");
  DynamicVector<T> result(v1.size());
  for (Index i = 0; i < v1.size(); ++i)
    result[i] = v1[i] * v2[i];
  return result;
}

// Unary minus
template <typename T> DynamicVector<T> operator-(const DynamicVector<T> &v) {
  return T{-1} * v;
}

// Scalar multiplication (vector / scalar)
template <typename T>
DynamicVector<T> operator/(const DynamicVector<T> &v, const T &k) {
  assert(!approximatelyEqualAbsRel(k, T{0}) &&
         "Can't divide by a number approximate to 0.");
  return v * (T{1} / k);
}

// Scalar divide (scalar/vector)
template <typename T>
DynamicVector<T> operator/(const T &k, const DynamicVector<T> &v) {
  DynamicVector<T> result(v.size());
  for (Index i = 0; i < v.size(); ++i) {
    assert(!approximatelyEqualAbsRel(v[i], T{}) && "Division by zero!");
    result[i] = k / v[i];
  }
  return result;
}

// Element-wised divide
template <typename T>
DynamicVector<T> operator/(const DynamicVector<T> &v1,
                           const DynamicVector<T> &v2) {
  assert(v1.size() == v2.size() &&
         "Can't perfom element-wised division on 2 different-sized vectors!");
  return v1 * (T{1} / v2);
}

// DynamicVector subtraction
template <typename T>
DynamicVector<T> operator-(const DynamicVector<T> &v1,
                           const DynamicVector<T> &v2) {
  return v1 + (-v2);
}

// Dot product
template <typename T>
T dotProduct(const DynamicVector<T> &v1, const DynamicVector<T> &v2) {
  if (v1.size() != v2.size())
    throw std::invalid_argument("Vectors must have the same dimension.");

  T result{};
  for (Index i = 0; i < v1.size(); ++i)
    result += v1[i] * v2[i];
  return result;
}

// Cross product (3D)
template <typename T>
DynamicVector<T> crossProduct(const DynamicVector<T> &v1,
                              const DynamicVector<T> &v2) {
  if (v1.size() != 3 || v2.size() != 3)
    throw std::invalid_argument(
        "Cross product is only defined for 3D vectors.");
  // Method 1
  DynamicVector<T> result(3);
  result[0] = v1[1] * v2[2] - v1[2] * v2[1];
  result[1] = v1[2] * v2[0] - v1[0] * v2[2];
  result[2] = v1[0] * v2[1] - v1[1] * v2[0];
  return result;
}

// Cross product (3D) - Method 2 using Levi-Civita
template <typename T>
DynamicVector<T> crossProduct2(const DynamicVector<T> &v1,
                               const DynamicVector<T> &v2) {
  if (v1.size() != 3 || v2.size() != 3)
    throw std::invalid_argument(
        "Cross product is only defined for 3D vectors.");

  // Method 2: Using Levi-Civita symbol
  DynamicVector<T> result(3);
  for (Index i = 0; i < 3; ++i) {
    result[i] = T{0}; // Initialize to zero
    for (Index j = 0; j < 3; ++j) {
      for (Index k = 0; k < 3; ++k) {
        result[i] += leviCivita(i, j, k) * v1[j] * v2[k];
      }
    }
  }
  return result;
}

// DynamicVector magnitude/norm
template <typename T> T magnitude(const DynamicVector<T> &v) {
  return std::sqrt(dotProduct(v, v));
}

// Normalize to Unit vector
template <typename T> DynamicVector<T> normalize(const DynamicVector<T> &v) {
  T mag = magnitude(v);
  if (approximatelyEqualAbsRel(mag, T{0}))
    throw std::invalid_argument("Cannot normalize zero vector.");
  return v * (T{1} / mag);
}

// Angle functions
template <typename T>
T angleRad(const DynamicVector<T> &v1, const DynamicVector<T> &v2) {
  return std::acos(dotProduct(v1, v2) / (magnitude(v1) * magnitude(v2)));
}

template <typename T>
T angleDegree(const DynamicVector<T> &v1, const DynamicVector<T> &v2) {
  return angleRad(v1, v2) * T{180} / T{constants::pi};
}

template <typename T>
bool isPerpendicular(const DynamicVector<T> &v1, const DynamicVector<T> &v2) {
  return approximatelyEqualAbsRel(dotProduct(v1, v2), T{});
}

template <typename T>
bool isParallel(const DynamicVector<T> &v1, const DynamicVector<T> &v2) {
  return approximatelyEqualAbsRel(magnitude(crossProduct(v1, v2)), T{});
}

// flatten: 2D -> 1D (row-major)
template <typename T>
DynamicVector<T> flatten(const DynamicVector<DynamicVector<T>> &twoDVect,
                         Index nRows, Index nCols) {
  if (static_cast<Index>(twoDVect.size()) != nRows)
    throw std::invalid_argument("flatten: nRows mismatch");
  for (Index r = 0; r < nRows; ++r)
    if (static_cast<Index>(twoDVect[r].size()) != nCols)
      throw std::invalid_argument("flatten: nCols mismatch");

  DynamicVector<T> flat(nRows * nCols);
  for (Index i = 0; i < nRows; ++i)
    for (Index j = 0; j < nCols; ++j)
      flat[i * nCols + j] = twoDVect[i][j];
  return flat;
}

// unflatten: 1D -> 2D (row-major)
template <typename T>
DynamicVector<DynamicVector<T>> unflatten(const DynamicVector<T> &oneDVect,
                                          Index nRows, Index nCols) {
  if (static_cast<Index>(oneDVect.size()) != nRows * nCols)
    throw std::invalid_argument("unflatten: size mismatch");
  DynamicVector<DynamicVector<T>> twoDVect(nRows);
  for (Index i = 0; i < nRows; ++i) {
    twoDVect[i].resize(static_cast<std::size_t>(nCols));
    for (Index j = 0; j < nCols; ++j)
      twoDVect[i][j] = oneDVect[i * nCols + j];
  }
  return twoDVect;
}

///////////////////////////////////////////////////
/////////////// STATIC CLASS VECTOR ///////////////
///////////////////////////////////////////////////
template <typename T, Index n> class StaticVector {
private:
  static_assert(n > 0 && "Number of rows and columns must be greater than 0!");
  std::array<T, n> m_elements{};

public:
  // Default Constructors
  StaticVector() = default;
  StaticVector(const StaticVector &) = default;
  StaticVector(StaticVector &&) = default;
  StaticVector &operator=(const StaticVector &) = default;
  StaticVector &operator=(StaticVector &&) = default;
  ~StaticVector() = default;

  // Constructor with initializer_list
  StaticVector(std::initializer_list<T> list) {
    assert(static_cast<Index>(list.size()) == n);
    std::copy(list.begin(), list.end(), m_elements.begin());
  }

  // Explicit conversion(static_cast/initialization) from DynamicVector.
  explicit StaticVector(const DynamicVector<T> &v) {
    if (v.size() != n)
      throw std::invalid_argument(
          "StaticVector size must match DynamicVector size.");
    for (Index i = 0; i < n; ++i)
      m_elements[static_cast<std::size_t>(i)] = v[i];
  }

  explicit operator DynamicVector<T>() const {
    DynamicVector<T> out(static_cast<std::size_t>(n));
    for (Index i = 0; i < n; ++i)
      out[i] = m_elements[static_cast<std::size_t>(i)];
    return out;
  }

  // // Constructor with length n of elements which are value
  explicit StaticVector(const T &value) {
    for (auto &e : m_elements)
      e = value;
  }

  // StaticVector 0
  static StaticVector<T, n> zero() { return StaticVector<T, n>{}; }

  // Reset all vallue to 0 default
  void resetZero() {
    for (auto &e : m_elements) {
      if constexpr (requires(T &x) { x.resetZero(); }) {
        e.resetZero();
      } else {
        e = T{};
      }
    }
  }

  // Signed indexing
  auto &operator[](Index i) {
    assert(i >= 0 && "Negative index not allowed");
    assert(static_cast<std::size_t>(i) < n && "Index out of bounds");
    return m_elements[static_cast<std::size_t>(i)];
  }

  const auto &operator[](Index i) const {
    assert(i >= 0 && "Negative index not allowed");
    assert(static_cast<std::size_t>(i) < n && "Index out of bounds");
    return m_elements[static_cast<std::size_t>(i)];
  }

  // Size functions
  constexpr Index size() const noexcept { return static_cast<Index>(n); }

  // Convenience coordinate accessors for 2D/3D StaticVectors
  T &x() {
    static_assert(n >= 1 && "x() requires size() >= 1");
    return (*this)[0];
  }
  const T &x() const {
    static_assert(n >= 1 && "x() requires size() >= 1");
    return (*this)[0];
  }

  T &y() {
    static_assert(n >= 2 && "y() requires size() >= 2");
    return (*this)[1];
  }
  const T &y() const {
    static_assert(n >= 2 && "y() requires size() >= 2");
    return (*this)[1];
  }

  T &z() {
    static_assert(n >= 3 && "z() requires size() >= 3");
    return (*this)[2];
  }
  const T &z() const {
    static_assert(n >= 3 && "z() requires size() >= 3");
    return (*this)[2];
  }

  // Iterators
  auto begin() { return m_elements.begin(); }
  auto end() { return m_elements.end(); }
  auto begin() const { return m_elements.begin(); }
  auto end() const { return m_elements.end(); }

  // Print
  friend std::ostream &operator<<(std::ostream &os,
                                  const StaticVector<T, n> &v) {
    os << "[";
    for (Index i = 0; i < n; ++i) {
      if (approximatelyEqualAbsRel(v[i], T{0}))
        os << T{0};
      else
        os << v[i];
      if (i < n - 1)
        os << ", ";
    }
    os << "]";
    return os;
  }

  //+= operator
  StaticVector &operator+=(const StaticVector<T, n> &other) {
    for (Index i{0}; i < n; ++i)
      (*this)[i] += other[i];
    return *this;
  }

  StaticVector &operator-=(const StaticVector<T, n> &other) {
    return *this += -other;
  }

  // Projection on another StaticVector
  StaticVector projection(const StaticVector &other) const {
    StaticVector<T, n> unit = normalize(other);
    return unit * dotProduct(*this, unit);
  }

  // Sort in-place. Default: descending order.
  void sort(bool descending = true) {
    if (descending) {
      std::sort(this->begin(), this->end(),
                [](const T &a, const T &b) { return a > b; });
    } else {
      std::sort(this->begin(), this->end(),
                [](const T &a, const T &b) { return a < b; });
    }
  }

  // Return a sorted copy. Default: descending order.
  StaticVector sorted(bool descending = true) const {
    StaticVector copy = *this;
    copy.sort(descending);
    return copy;
  }

  StaticVector &operator=(const DynamicVector<T> &v) {
    if (v.size() != n)
      throw std::invalid_argument(
          "StaticVector size must match DynamicVector size.");
    for (Index i = 0; i < n; ++i)
      m_elements[static_cast<std::size_t>(i)] = v[i];
    return *this;
  }

  Matrix<T, 1, n> transpose() const {
    Matrix<T, 1, n> result{};
    for (Index i = 0; i < n; ++i)
      result(0, i) = (*this)[i];
    return result;
  }
};

/////////////////////// END OF STATICVECTOR CLASS/////////////////////
// StaticVector operators
template <typename T, Index n1, Index n2>
bool operator==(const StaticVector<T, n1> &v1, const StaticVector<T, n2> &v2) {
  if (n1 != n2)
    return false;
  for (Index i = 0; i < n1; ++i) {
    if (!approximatelyEqualAbsRel(v1[i], v2[i]))
      return false;
  }
  return true;
}

template <typename T, Index n>
bool operator!=(const StaticVector<T, n> &v1, const StaticVector<T, n> &v2) {
  return !(v1 == v2);
}

template <typename T, Index n>
StaticVector<T, n> operator+(const StaticVector<T, n> &v1,
                             const StaticVector<T, n> &v2) {
  StaticVector<T, n> result{};
  for (Index i = 0; i < n; ++i)
    result[i] = v1[i] + v2[i];
  return result;
}

// Scalar multiplication (scalar * StaticVector)
template <typename T, Index n>
StaticVector<T, n> operator*(const T &k, const StaticVector<T, n> &v) {
  StaticVector<T, n> result{};
  for (Index i = 0; i < n; ++i)
    result[i] = k * v[i];
  return result;
}

// Scalar multiplication (StaticVector * scalar)
template <typename T, Index n>
StaticVector<T, n> operator*(const StaticVector<T, n> &v, const T &k) {
  return k * v;
}

// Element-wise multiplication
template <typename T, Index n>
StaticVector<T, n> operator*(const StaticVector<T, n> &v1,
                             const StaticVector<T, n> &v2) {
  StaticVector<T, n> result{};
  for (Index i = 0; i < n; ++i)
    result[i] = v1[i] * v2[i];
  return result;
}

// Unary minus
template <typename T, Index n>
StaticVector<T, n> operator-(const StaticVector<T, n> &v) {
  return T{-1} * v;
}

// Scalar multiplication (StaticVector / scalar)
template <typename T, Index n>
StaticVector<T, n> operator/(const StaticVector<T, n> &v, const T &k) {
  assert(!approximatelyEqualAbsRel(k, T{0}) &&
         "Can't divide by a number approximate to 0.");
  return v * (T{1} / k);
}

// Scalar divide (scalar/StaticVector)
template <typename T, Index n>
StaticVector<T, n> operator/(const T &k, const StaticVector<T, n> &v) {
  StaticVector<T, n> result{};
  for (Index i = 0; i < n; ++i) {
    assert(!approximatelyEqualAbsRel(v[i], T{}) && "Division by zero!");
    result[i] = k / v[i];
  }
  return result;
}

// Element-wised divide
template <typename T, Index n>
StaticVector<T, n> operator/(const StaticVector<T, n> &v1,
                             const StaticVector<T, n> &v2) {
  return v1 * (T{1} / v2);
}

// StaticVector subtraction
template <typename T, Index n>
StaticVector<T, n> operator-(const StaticVector<T, n> &v1,
                             const StaticVector<T, n> &v2) {
  return v1 + (-v2);
}

// Dot product
template <typename T, Index n>
T dotProduct(const StaticVector<T, n> &v1, const StaticVector<T, n> &v2) {
  T result{};
  for (Index i = 0; i < n; ++i)
    result += v1[i] * v2[i];
  return result;
}

// Quadratic form: v^T A v
template <typename T, Index n>
T quadraticForm(const StaticVector<T, n> &v, const Matrix<T, n, n> &A) {
  const StaticVector<T, 1> q = (v.transpose() * A) * v;
  return q[0];
}

template <typename T, Index n>
T quadraticForm(const DynamicVector<T> &v, const Matrix<T, n, n> &A) {
  if (v.size() != n)
    throw std::invalid_argument(
        "quadraticForm: vector size must match matrix dimensions.");

  T result{};
  for (Index i = 0; i < n; ++i) {
    T rowDot{};
    for (Index j = 0; j < n; ++j)
      rowDot += A(i, j) * v[j];
    result += v[i] * rowDot;
  }
  return result;
}

// Cross product (3D)
template <typename T, Index n>
StaticVector<T, n> crossProduct(const StaticVector<T, n> &v1,
                                const StaticVector<T, n> &v2) {
  static_assert(n == 3, "Cross product only for 3D StaticVectors");
  // Method 1
  StaticVector<T, 3> result{};
  result[0] = v1[1] * v2[2] - v1[2] * v2[1];
  result[1] = v1[2] * v2[0] - v1[0] * v2[2];
  result[2] = v1[0] * v2[1] - v1[1] * v2[0];
  return result;
}

// Cross product (3D) - Method 2 using Levi-Civita
template <typename T, Index n>
StaticVector<T, n> crossProduct2(const StaticVector<T, n> &v1,
                                 const StaticVector<T, n> &v2) {
  static_assert(n == 3, "Cross product only for 3D StaticVectors");
  // Method 2: Using Levi-Civita symbol
  StaticVector<T, 3> result{};
  for (Index i = 0; i < 3; ++i) {
    result[i] = T{0}; // Initialize to zero
    for (Index j = 0; j < 3; ++j) {
      for (Index k = 0; k < 3; ++k) {
        result[i] += leviCivita(i, j, k) * v1[j] * v2[k];
      }
    }
  }
  return result;
}

// StaticVector magnitude/norm
template <typename T, Index n> T magnitude(const StaticVector<T, n> &v) {
  return std::sqrt(dotProduct(v, v));
}

// Normalize to Unit StaticVector
template <typename T, Index n>
StaticVector<T, n> normalize(const StaticVector<T, n> &v) {
  T mag = magnitude(v);
  if (approximatelyEqualAbsRel(mag, T{0}))
    throw std::invalid_argument("Cannot normalize zero StaticVector.");
  return v * (T{1} / mag);
}

// Angle functions
template <typename T, Index n>
T angleRad(const StaticVector<T, n> &v1, const StaticVector<T, n> &v2) {
  return std::acos(dotProduct(v1, v2) / (magnitude(v1) * magnitude(v2)));
}

template <typename T, Index n>
T angleDegree(const StaticVector<T, n> &v1, const StaticVector<T, n> &v2) {
  return angleRad(v1, v2) * T{180} / T{constants::pi};
}

template <typename T, Index n>
bool isPerpendicular(const StaticVector<T, n> &v1,
                     const StaticVector<T, n> &v2) {
  return approximatelyEqualAbsRel(dotProduct(v1, v2), T{});
}

template <typename T, Index n>
bool isParallel(const StaticVector<T, n> &v1, const StaticVector<T, n> &v2) {
  return approximatelyEqualAbsRel(magnitude(crossProduct(v1, v2)), T{});
}

#endif