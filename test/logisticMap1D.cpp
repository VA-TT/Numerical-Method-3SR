#include "../library/DualDifferentiation.h"
#include <iostream>

// x_n+1 = muy.x_n.(1 - x_n)
double logisticMap(Dual x) {
  double muy = 1.0; // 0 < muy < 4
  double x0 = 1.0;  // 0 < x0 < 1
  Dual f = muy * x * (1.0 - x);
  if (f == x)

}
int main() {}
