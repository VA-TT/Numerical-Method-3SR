#ifndef MY_VECTOR_CLASS
#define MY_VECTOR_CLASS

#include "comparison.h" //Approximative Comparsion
#include "kroneckerDelta_LeviCivita.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using Index = std::ptrdiff_t; // typedef

template <typename T> class Vector {
private:
  std::vector<T> m_elements{};

public:
  // Default Constructors
  Vector() = default;
  Vector(const Vector &) = default;
  Vector(Vector &&) = default;
  Vector &operator=(const Vector &) = default;
  Vector &operator=(Vector &&) = default;
  ~Vector() = default;

  // Constructor with length n of elements which are 0
  explicit Vector(std::size_t n) : m_elements(n, T{}) {}

  // // Constructor with length n of elements which are value
  // Vector(std::size_t n, const T& value) : m_elements(n, value) {}

  // Constructor với initializer_list
  Vector(std::initializer_list<T> list) : m_elements(list) {}

  // // Constructor 3D (x, y, z)
  // Vector(const T& x, const T& y, const T& z) : m_elements{x, y, z} {}

  // Signed indexing với assert
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

  // Resize vector
  void resize(std::size_t n) { m_elements.resize(n); }
  void resize(std::size_t n, const T &value) { m_elements.resize(n, value); }

  // Add element
  void push_back(const T &value) { m_elements.push_back(value); }

  // Print
  friend std::ostream &operator<<(std::ostream &os, const Vector<T> &v) {
    os << "[";
    for (Index i = 0; i < v.size(); ++i) {
      os << v[i];
      if (i < v.size() - 1)
        os << ", ";
    }
    os << "]";
    return os;
  }

  //+= operator
  Vector &operator+=(const Vector<T> &other) {
    if (this->size() != other.size())
      throw std::invalid_argument(
          "Vectors must have the same dimension for operator+=.");
    for (Index i{0}; i < this->size(); ++i)
      (*this)[i] += other[i];
    return *this;
  }
  // Projection on another vector
  Vector<T> projection(const Vector &other) const {
    return normalize(other) * dotProduct(*this, normalize(other));
  }
};

/////////////////////////////////// END OF VECTOR CLASS
//////////////////////////////////////////
// Vector operators
template <typename T>
bool operator==(const Vector<T> &v1, const Vector<T> &v2) {
  if (v1.size() != v2.size())
    return false;
  for (Index i = 0; i < v1.size(); ++i) {
    if (!approximatelyEqualAbsRel(v1[i], v2[i]))
      return false;
  }
  return true;
}

template <typename T>
bool operator!=(const Vector<T> &v1, const Vector<T> &v2) {
  return !(v1 == v2);
}

template <typename T>
Vector<T> operator+(const Vector<T> &v1, const Vector<T> &v2) {
  if (v1.size() != v2.size())
    throw std::invalid_argument("Vectors must have the same dimension.");

  Vector<T> result(v1.size());
  for (Index i = 0; i < v1.size(); ++i)
    result[i] = v1[i] + v2[i];
  return result;
}

// Scalar multiplication (scalar * vector)
template <typename T> Vector<T> operator*(const T &k, const Vector<T> &v) {
  Vector<T> result(v.size());
  for (Index i = 0; i < v.size(); ++i)
    result[i] = k * v[i];
  return result;
}

// Scalar multiplication (vector * scalar)
template <typename T> Vector<T> operator*(const Vector<T> &v, const T &k) {
  return k * v;
}

// Unary minus
template <typename T> Vector<T> operator-(const Vector<T> &v) {
  return T{-1} * v;
}

// Scalar multiplication (vector / scalar)
template <typename T> Vector<T> operator/(const Vector<T> &v, const T &k) {
  assert(!approximatelyEqualAbsRel(k, 0.0) &&
         "Can't divide by a number approximate to 0.");
  return v * (1 / k);
}

// Vector subtraction
template <typename T>
Vector<T> operator-(const Vector<T> &v1, const Vector<T> &v2) {
  return v1 + (-v2);
}

// Dot product
template <typename T> T dotProduct(const Vector<T> &v1, const Vector<T> &v2) {
  if (v1.size() != v2.size())
    throw std::invalid_argument("Vectors must have the same dimension.");

  T result{};
  for (Index i = 0; i < v1.size(); ++i)
    result += v1[i] * v2[i];
  return result;
}

// Cross product (3D)
template <typename T>
Vector<T> crossProduct(const Vector<T> &v1, const Vector<T> &v2) {
  if (v1.size() != 3 || v2.size() != 3)
    throw std::invalid_argument(
        "Cross product is only defined for 3D vectors.");
  // Method 1
  Vector<T> result(3);
  result[0] = v1[1] * v2[2] - v1[2] * v2[1];
  result[1] = v1[2] * v2[0] - v1[0] * v2[2];
  result[2] = v1[0] * v2[1] - v1[1] * v2[0];
  return result;
}

// Cross product (3D) - Method 2 using Levi-Civita
template <typename T>
Vector<T> crossProduct2(const Vector<T> &v1, const Vector<T> &v2) {
  if (v1.size() != 3 || v2.size() != 3)
    throw std::invalid_argument(
        "Cross product is only defined for 3D vectors.");

  // Method 2: Using Levi-Civita symbol
  Vector<T> result(3);
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

// Vector magnitude/norm
template <typename T> T magnitude(const Vector<T> &v) {
  return std::sqrt(dotProduct(v, v));
}

// Normalize to Unit vector
template <typename T> Vector<T> normalize(const Vector<T> &v) {
  T mag = magnitude(v);
  if (approximatelyEqualAbsRel(mag, T{0.0}))
    throw std::invalid_argument("Cannot normalize zero vector.");
  return v * (T{1.0} / mag);
}

// Angle functions
template <typename T> T angleRad(const Vector<T> &v1, const Vector<T> &v2) {
  return std::acos(dotProduct(v1, v2) / (magnitude(v1) * magnitude(v2)));
}

template <typename T> T angleDegree(const Vector<T> &v1, const Vector<T> &v2) {
  const T PI = T{3.14159265358979323846};
  return angleRad(v1, v2) * T{180} / PI;
}

template <typename T>
bool isPerpendicular(const Vector<T> &v1, const Vector<T> &v2) {
  return approximatelyEqualAbsRel(dotProduct(v1, v2), T{});
}

template <typename T>
bool isParallel(const Vector<T> &v1, const Vector<T> &v2) {
  return approximatelyEqualAbsRel(magnitude(crossProduct(v1, v2)), T{});
}

#endif