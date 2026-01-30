// Range of 1D reference element Sr: [-1,1]
// First order polynomial shape function
auto N1_r = [](auto psy) { return (1 + psy) / 2; };
auto N2_r = [](auto psy) { return (1 - psy) / 2; };

// Discretizing x = phi(psy) = x1 N1(psy) + x2 N2(psy) = J*psy + M = psy *
// (x2-x1)/2 + (x2+x1)/2

// Jacobian: J = dx/d(psy) = (x2-x1)/2
template <typename T> T Jacobian1D(T x1, T x2) { return (x2 - x1) / 2; }
template <typename T> T midPoint(T x1, T x2) { return (x2 + x1) / 2; }

// Mapping deriviation: dN/dx = dN/d(psy) * d(psy)/dx = dN/d(psy) * 1/J
template <typename T> T dN1(T x1, T x2) { return 0.5 / Jacobian1D(x1, x2); }
template <typename T> T dN2(T x1, T x2) { return -0.5 / Jacobian1D(x1, x2); }

// Function Phi(psy) used to map
template <typename T> T physicCoor(T psy, T x1, T x2) {
  return N1_r(psy) * x1 + N2_r(psy) * x2;
}
template <typename T> T parentCoor(T x, T x1, T x2) {
  return (x - midPoint(x1, x2)) / Jacobian1D(x1, x2);
}
// Interpolation: int^(x2)_(x1) f(x)dx = int^1_(-1) f(phi(psy)) J d(psy)

// #include <iostream>
// int main() {
//   double x = 4.0, x1 = 2.0, x2 = 5.0, psy = 0.0;
//   std::cout << "Ref coordinate of x=4: psy = " << parentCoor(x, x1, x2) <<
//   '\n'; std::cout << "Physical coordinate of psy=0: x = " << physicCoor(psy,
//   x1, x2)
//             << '\n';
//   std::cout << "Jacobian: " << Jacobian1D(x1, x2) << '\n';
//   std::cout << "Derivatives of N1: dN1 = " << dN1(x1, x2) << '\n';
//   std::cout << "Derivatives of N2: dN2 = " << dN2(x1, x2) << '\n';
// }