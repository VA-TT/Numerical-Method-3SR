#include "../library/Matrix.h"
#include "../library/Vector.h"
#include "../library/physicConstants.h"
#include <cmath>
#include <iostream>

int main() {

  // Rotation dans le sens horaire, la compression est positive
  double phi = 35.0;
  double c = 1000.0;

  double alpha = 55;

  double nx = std::sin(alpha / 180.0 * constants::pi);
  double ny = std::cos(alpha / 180.0 * constants::pi);
  Matrix<double, 2, 1> n{nx, ny};
  Matrix<double, 2, 1> t{-ny, nx};
  Matrix<double, 2, 2> sigmaCauchy = {
      1.0290172785590373e+03, -1.8256989997892878e+00, -1.8256989997892878e+00,
      -4.8654819189809857e+01};
  Matrix<double, 2, 2> sigma = -sigmaCauchy;

  const double sigma_n = (n.transpose() * sigma * n)(0, 0);
  const double tau_nt = (t.transpose() * sigma * n)(0, 0);

  double phi_rad = phi / 180.0 * constants::pi;

  const double tau_max = std::tan(phi_rad) * sigma_n + c;

  double beta = alpha + 90.0;
  double sigma_u = (sigma.xx() + sigma.yy()) / 2.0 +
                   (sigma.xx() - sigma.yy()) / 2.0 *
                       std::cos(2.0 * beta / 180.0 * constants::pi) -
                   sigma.xy() * std::sin(2.0 * beta / 180.0 * constants::pi);
  double tau_uv = (sigma.xx() - sigma.yy()) / 2.0 *
                      std::sin(2.0 * beta / 180.0 * constants::pi) +
                  sigma.xy() * std::cos(2.0 * beta / 180.0 * constants::pi);

  std::cout << "sigma_n = " << sigma_n << "\n";
  std::cout << "tau_nt  = " << tau_nt << "\n";
  //   std::cout << "sigma_u = " << sigma_u << "\n";
  //   std::cout << "tau_uv  = " << tau_uv << "\n";
  std::cout << "tau_max = " << tau_max << "\n";

  std::cout << "Critere MC: |tau_nt| - tau_max = " << (tau_nt - tau_max)
            << std::endl;
  std::cout << "Ecart relatif: (|tau_nt| - tau_max)/tau_max = "
            << (std::abs(tau_nt) - tau_max) / tau_max * 100.0 << "%"
            << std::endl;

  auto [principalSigmas, principalDirection] = sigma.eigen();
  principalSigmas.sort(true);
  std::cout << "Contraintes principales (sigma1 >= sigma2): " << principalSigmas
            << "\n";

  const double sigma_1 = principalSigmas[0];
  const double sigma_2 = principalSigmas[1];

  // Analytic principal stresses for 2D symmetric tensor
  // sigma_{1,2} = (sx+sy)/2 ± sqrt( ((sx-sy)/2)^2 + tau_xy^2 )
  //   const double sx = sigma.xx();
  //   const double sy = sigma.yy();
  //   const double txy = sigma.xy();
  //   const double s_avg = 0.5 * (sx + sy);
  //   const double radius = std::sqrt(std::pow(0.5 * (sx - sy), 2.0) + txy *
  //   txy); const double sigma_1_analytic2D = s_avg + radius; const double
  //   sigma_2_analytic2D = s_avg - radius;

  const double Kp =
      std::pow(std::tan(constants::pi / 4.0 + phi_rad / 2.0), 2.0);
  const double sigma_1_analytic =
      sigma_2 * Kp + 2.0 * c * std::cos(phi_rad) / (1.0 - std::sin(phi_rad));

  //   std::cout << "Kp = tan^2(pi/4 + phi/2) = " << Kp << "\n";
  std::cout << "sigma_1 principal      = " << sigma_1 << "\n";
  //   std::cout << "sigma_1 analytique (2D) = " << sigma_1_analytic2D << "\n";
  std::cout << "sigma_2 principal      = " << sigma_2 << "\n";
  //   std::cout << "sigma_2 analytique (2D) = " << sigma_2_analytic2D << "\n";
  std::cout << "sigma_1_analytique (MC)     = " << sigma_1_analytic << "\n";
  std::cout << "Ecart relatif: (sigma1_principal - sigma1_MC)/sigma1_MC = "
            << (sigma_1 - sigma_1_analytic) / sigma_1_analytic * 100.0 << "%"
            << std::endl;
  std::cout << "Angle entre direction principal et x: "
            << angleDegree(principalDirection.getColVector(0),
                           StaticVector<double, 2>{1.0, 0.0})
            << "\n";
}