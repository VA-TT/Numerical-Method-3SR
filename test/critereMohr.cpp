#include "../library/Matrix.h"
#include "../library/Vector.h"
#include "../library/physicConstants.h"
#include <cmath>
#include <iostream>

int main() {
  double alpha = 40.0;
  double nx = std::sin(alpha / 180.0 * constants::pi);
  double ny = std::cos(alpha / 180.0 * constants::pi);
  Matrix<double, 2, 1> n{nx, ny};
  Matrix<double, 2, 1> t{-ny, nx};
  Matrix<double, 2, 2> sigma = {1.5026016179884286e+03, 1.1995347639438139e+03,
                                1.1995347639438139e+03, 2.1350627645031122e+03};

  const double sigma_n = (n.transpose() * sigma * n)(0, 0);
  const double tau_nt = (t.transpose() * sigma * n)(0, 0);

  double phi = 43.0;
  double c = 1.0;

  double phi_rad = phi / 180.0 * constants::pi;

  const double tau_max = std::tan(phi_rad) * sigma_n + c;

  std::cout << "sigma_n = " << sigma_n << "\n";
  std::cout << "tau_nt  = " << tau_nt << "\n";
  std::cout << "tau_max = " << tau_max << "\n";
  std::cout << "Erreur (tau_nt - tau_max) : " << (tau_nt - tau_max)
            << std::endl;

  auto [principalSigmas, eigenVectors] = sigma.eigen();
  (void)eigenVectors;
  principalSigmas.sort(true);
  std::cout << "Contraintes principales (sigma1 >= sigma2): " << principalSigmas
            << "\n";

  const double sigma_1 = principalSigmas[0];
  const double sigma_2 = principalSigmas[1];

  const double Kp =
      std::pow(std::tan(constants::pi / 4.0 + phi_rad / 2.0), 2.0);
  const double sigma_1_analytic =
      sigma_2 * Kp + 2.0 * c * std::cos(phi_rad) / (1.0 - std::sin(phi_rad));

  std::cout << "Kp = tan^2(pi/4 + phi/2) = " << Kp << "\n";
  std::cout << "sigma_1 (from eigen)      = " << sigma_1 << "\n";
  std::cout << "sigma_1_analytic (MC)     = " << sigma_1_analytic << "\n";
}