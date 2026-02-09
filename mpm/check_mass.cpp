#include "../library/MPM.cpp"
#include <iostream>

int main() {
    double E = 4.0 * 3.14159265358979323846 * 3.14159265358979323846;
    double rho = 1.0;
    double L = 1.0;
    double v0 = 0.1;
    double dt = 0.01;
    double duration = 10.0;
    double xloc = 1.0;
    
    MPM1D<double, 2, 1> beam(E, rho, L, v0, dt, duration, xloc);
    beam.setupMP();
    
    std::cout << "mass_p[0] = " << beam.getMPmass(0) << "\n";
    std::cout << "volume_p[0] = " << beam.getMPvolume(0) << "\n";
    std::cout << "position_p[0] = " << beam.getMPposition(0) << "\n";
    
    return 0;
}
