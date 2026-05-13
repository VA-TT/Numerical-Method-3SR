#ifndef CONSTEXPR_SIGN_FUNCTION
#define CONSTEXPR_SIGN_FUNCTION

#include "comparison.h"

template <typename T> constexpr T sgn(T x) {
  return approximatelyEqualAbsRel(x, T{0}) ? T{0} : ((x < T{0}) ? T{-1} : T{1});
}

#endif