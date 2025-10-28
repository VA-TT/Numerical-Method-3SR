#ifndef MY_KRONECKER_DELTA
#define MY_KRONECKER_DELTA

constexpr inline int kroneckerDelta(int i, int j) { return (i == j) ? 1 : 0; }

constexpr inline int leviCivita(int i, int j, int k) {
  return (i - j) * (j - k) * (k - i) / 2;
}

#endif