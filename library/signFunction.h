#ifndef CONSTEXPR_SIGN_FUNCTION
#define CONSTEXPR_SIGN_FUNCTION

#include "comparison.h"

constexpr double sgn(double x) {
  return approximatelyEqualAbsRel(x, 0.0) ? 0.0 : ((x < 0) ? -1.0 : 1.0);
}

#endif