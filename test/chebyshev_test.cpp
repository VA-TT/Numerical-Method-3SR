#include <cmath>
#include <iomanip>
#include <iostream>

// include the interpolation utilities (guarded)
#include "../library/interpolate.h"

int main() {
  std::cout << std::fixed << std::setprecision(10);

  auto run_for_n = [&](int n) {
    int nodes = n + 1; // degree n => n+1 nodes
    Vector<double> xs(nodes);
    Vector<double> ys(nodes);
    for (int k = 0; k < nodes; ++k) {
      double x = std::cos(M_PI * k /
                          double(n)); // Chebyshev-Gauss-Lobatto nodes in [-1,1]
      xs[k] = x;
      ys[k] = 1.0 / (1.0 + 25.0 * x * x);
    }

    int pts = 2001;
    double max_err = 0.0;
    double sum_sq = 0.0;
    for (int i = 0; i < pts; ++i) {
      double x = -1.0 + 2.0 * i / double(pts - 1);
      double approx = interpolatePolynomial(xs, ys, x);
      double exact = 1.0 / (1.0 + 25.0 * x * x);
      double err = std::fabs(approx - exact);
      if (err > max_err)
        max_err = err;
      sum_sq += err * err;
    }
    double rms = std::sqrt(sum_sq / pts);
    std::cout << "n=" << n << " nodes=" << nodes << " max_err=" << max_err
              << " rms=" << rms << "\n";
  };

  int tests[] = {5, 10, 20, 40};
  for (int t : tests)
    run_for_n(t);

  return 0;
}
