#ifndef MY_NEWTON_RAPHSON_ROOT
#define MY_NEWTON_RAPHSON_ROOT

#include <cmath>

inline double myFunc(double x) {
  return (8 * std::pow(x, 3) - 4 * std::pow(x, 2) - 3 * x + 1);
}

inline double derivative(double x) { return (24 * std::pow(x, 2) - 8 * x - 3); }

inline double newtonRaphson(double x0, double epsilon, int max_iter = 1000) {
  double x_n = x0;
  const double tiny = 1e-16; // ngưỡng để tránh chia cho 0
  for (int i = 0; i < max_iter; ++i) {
    double f_n = myFunc(x_n);
    double f_prime = derivative(x_n);

    if (std::abs(f_prime) < tiny) {
      // derivative too small: stop or fallback
      break;
    }

    double x_next = x_n - f_n / f_prime;

    // dừng theo bước (x) hoặc theo giá trị hàm (tuỳ bạn)
    if (std::abs(x_next - x_n) <= epsilon)
      return x_next;
    // or: if (std::abs(f_n) <= epsilon) return x_n;

    x_n = x_next;
  }
  return x_n; // trả x gần đúng sau max_iter
}

#endif