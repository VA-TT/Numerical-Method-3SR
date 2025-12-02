#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>

namespace mercuryConstants
{
    double a{0.39};
    double e{0.206};
    double GM{4 * std::pow(std::numbers::pi, 2)};
};

struct planet
{
    double x{};
    double y{};
    double vx{};
    double vy{};
    double ax{};
    double ay{};
    double r{};
};

void integrateEulerCromer(planet &p, double GM, double epsilon, double tmax)
{
    // Instant t
    p.vx += epsilon * p.ax;
    p.vy += epsilon * p.ay;
    p.x += epsilon * p.vx;
    p.y += epsilon * p.vy;

    // Instant t + epsilon
    p.r = std::sqrt(std::pow(p.x, 2) + std::pow(p.y, 2));
    p.ax = -GM / (std::pow(p.r, 3)) * p.x;
    p.ay = -GM / (std::pow(p.r, 3)) * p.y;
}

void integrateVerlet(planet &p, double GM, double epsilon, double tmax)
{
    // Instant t + epsilon

    p.x += p.vx * epsilon + 0.5 * p.ax * std::pow(epsilon, 2);
    p.y += p.vy * epsilon + 0.5 * p.ay * std::pow(epsilon, 2);

    double R{std::sqrt(std::pow(p.x, 2) + std::pow(p.y, 2))};
    double Ax{-GM / (std::pow(p.r, 3)) * p.x};
    double Ay{-GM / (std::pow(p.r, 3)) * p.y};

    p.vx += 0.5 * (Ax + p.ax) * epsilon;
    p.vy += 0.5 * (Ay + p.ay) * epsilon;

    p.ax = Ax;
    p.ay = Ay;
}

int main()
{
    // Condition initiale (aphelion) at t0:
    planet mercury{};
    mercury.x = mercuryConstants::a * (1 + mercuryConstants::e);
    mercury.y = 0.0;
    mercury.vx = 0.0;
    mercury.vy = std::sqrt(mercuryConstants::GM / mercuryConstants::a * (1 - mercuryConstants::e) / (1 + mercuryConstants::e));
    mercury.r = std::sqrt(std::pow(mercury.x, 2) + std::pow(mercury.y, 2));
    mercury.ax = -mercuryConstants::GM / (std::pow(mercury.r, 3)) * mercury.x;
    mercury.ay = -mercuryConstants::GM / (std::pow(mercury.r, 3)) * mercury.y;

    double epsilon = 0.001;
    double tmax = 0.28;
    std::ofstream foutEuler("orbitEuler.dat");
    // std::ofstream foutVerlet("orbitVerlet.dat");
    foutEuler << std::fixed << std::setprecision(15);
    // foutVerlet << std::fixed << std::setprecision(15);
    for (double t = 0; t < tmax; t += epsilon)
    {
        foutEuler << mercury.x << " " << mercury.y << "\n";
        // foutVerlet << mercury.x << " " << mercury.y << "\n";
        integrateEulerCromer(mercury, mercuryConstants::GM, epsilon, tmax);
        // integrateVerlet(mercury, mercuryConstants::GM, epsilon, tmax);
    }
    foutEuler.close();
    return 0;
}