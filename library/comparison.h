#ifndef COMPARISON_HPP
#define COMPARISON_HPP

#include <algorithm> // for std::max
#include <iostream>

// Our own constexpr implementation of std::abs (for use in C++14/17/20)
// In C++23, use std::abs
// constAbs() can be called like a normal function, but can handle different types of values (e.g. int, double, etc...)
template <typename T>
constexpr T constAbs(T x)
{
    return (x < 0 ? -x : x);
}

// Return true if the difference between a and b is within epsilon percent of the larger of a and b
constexpr bool approximatelyEqualRel(double a, double b, double relEpsilon = 1e-8)
{
    return (constAbs(a - b) <= (std::max(constAbs(a), constAbs(b)) * relEpsilon));
}

// Return true if the difference between a and b is less than or equal to absEpsilon, or within relEpsilon percent of the larger of a and b
constexpr bool approximatelyEqualAbsRel(double a, double b, double absEpsilon = 1e-12, double relEpsilon = 1e-8)
{
    // Check if the numbers are really close -- needed when comparing numbers near zero.
    if (constAbs(a - b) <= absEpsilon)
        return true;

    // Otherwise fall back to Knuth's algorithm
    return approximatelyEqualRel(a, b, relEpsilon);
}

namespace myConstants
{
    [[maybe_unused]] constexpr double relEps{1e-8};
    [[maybe_unused]] constexpr double absEps{1e-12};
}

#endif // MATRIX_COMPARISON_HPP

// Credit: https://www.learncpp.com/cpp-tutorial/relational-operators-and-floating-point-comparisons/

// Command example: constexpr bool same { approximatelyEqualAbsRel(a, 1.0, absEps, relEps) };