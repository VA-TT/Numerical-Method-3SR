#include "../library/DualDifferentiation.h"
#include "../library/Matrix.h" //Approximative Comparsion
#include "../library/Vector.h"
#include "../library/clock.h"
#include "../library/gaussQuadrature.h"
#include "../library/interpolate.h"
#include "../library/physicConstants.h"
#include <cassert> // for assert
#include <fstream> //working with files
#include <iomanip> //tab
#include <iostream>
#include <numbers>     // for std::numbers::pi
#include <stdexcept>   //throw exception
#include <type_traits> // precision
#include <vector>

double solution(double x);

namespace modelParameters {
// Problem's domain (a,b)
constexpr double a{0.0};
constexpr double b{1.0};
constexpr double L{b - a};
constexpr double H{1.0};   // 1D
constexpr double V{H * L}; // 1D

// Material property
constexpr double E{4 * constants::pi * constants::pi};
constexpr double rho{1.0}; // gonna change
constexpr double c{constexpr_sqrt(E / rho)};
double current_time{0.0};
constexpr double duration{10.0};
constexpr double dt{0.01};
constexpr double nSteps{duration / dt};
constexpr double w{1.0 / L * c};
constexpr double v0{0.1};
constexpr double x_loc{0.5};

constexpr Index nNodes{2};             // Numbers of nodes
constexpr Index nElements{nNodes - 1}; // Numbers of elements

// Material Points
constexpr Index nMPs{1}; // Numbers material points
double mass_p{rho * V / nMPs};
double x_p{0.5 * L}; // Location of MP
double V_p{V / nMPs};
// Initial loading condition
double v_p{v0};
double stress_p{0.0}, strain_p{0.0};
double mv_p{mass_p * v_p};                 // Momentum
double strain_rate_p{0.0}, dStrain_p{0.0}; // Momentum

// Initiate needed containers
DynamicVector<double> nodes{};
DynamicVector<double> positions{};      // Analytical positions
DynamicVector<double> velocities{};     // Analytical velocities
DynamicVector<double> positions_mpm{};  // MPM positions
DynamicVector<double> velocities_mpm{}; // MPM velocities
DynamicVector<double> times{};
// Start and end of elements, will be generated automatically in main
DynamicVector<Index> eleOrigin;
DynamicVector<Index> eleEnd;
DynamicVector<double> length(nElements);
DynamicVector<DynamicVector<double>> element(nElements);
DynamicVector<DynamicVector<std::function<double(double)>>> N(nNodes),
    B(nNodes);
DynamicVector<double> N_p(nNodes), B_p(nNodes);
DynamicVector<double> m_i(nNodes), v_i(nNodes), mv_i(nNodes);

// Nodal external forces
DynamicVector<double> f_ext_i(nNodes), f_int_i(nNodes), f_total_i(nNodes);

} // namespace modelParameters

// Analytic solution for uniform load: u(x) = (1/EA)*( f*(L*x - x^2/2) + F*x )
double velocity_solution(double x_loc, double t) {
  using namespace modelParameters;
  return v0 * std::cos(w * t);
}

double position_solution(double x_loc, double t) {
  using namespace modelParameters;
  return x_loc * std::exp(v0 / (L * w) * std::sin(w * t));
}

// 0. Mesh generated function
DynamicVector<double> generateMesh(double a, double b, Index n) {
  DynamicVector<double> nodes;
  for (Index i{0}; i < n; ++i) {
    nodes.push_back((b - a) / (n - 1) * i + a);
  }
  return nodes;
}

// 1. Generate Shape function on each nodes
void shapeFunction() {
  // Equally divied
  using namespace modelParameters;
  // Initialize element connectivity
  for (Index e{0}; e < nElements; ++e) {
    eleOrigin.push_back(e);
    eleEnd.push_back(e + 1);
  }
  int local_i = 0;
  int local_j = 1;
  for (Index e{0}; e < nElements; ++e) {
    int i = eleOrigin[e];
    int j = eleEnd[e];
    double x_i = nodes[i];
    double x_j = nodes[j];
    element[e] = {x_i, x_j};
    length[e] = constexpr_fabs(x_j - x_i);

    // MUST set parameter "x" here to type <<auto>> in order to accept Dual
    // class as input for derivative calculating
    auto shapefunction_i = [=](auto x) {
      return basisLagrange(local_i, element[e], x);
    };
    auto dShape_i = [=](double x) { return automaticDiff(shapefunction_i, x); };
    auto shapefunction_j = [=](auto x) {
      return basisLagrange(local_j, element[e], x);
    };
    auto dShape_j = [=](double x) { return automaticDiff(shapefunction_j, x); };
    for (Index kk{0}; kk < nNodes; ++kk) {
      if (kk == i) {
        N[kk].push_back(shapefunction_i);
        B[kk].push_back(dShape_i);
      } else if (kk == j) {
        N[kk].push_back(shapefunction_j);
        B[kk].push_back(dShape_j);
      } else {
        N[kk].push_back([](double x) { return 0.0; });
        B[kk].push_back([](double x) { return 0.0; });
      }
    }
  }
}

// 2. Map mass and momentum to nodes via shape function
// 3. Compute nodal forces
// 4. Applying Boundary condition
// 5. Node to particle

int main() {
  // Equally divied
  Timer t;
  using namespace modelParameters;

  double dt_crit = L / c;
  assert((dt_crit / 10.0) >= dt &&
         "Time step isn't satisfied CFL condition (too big)");
  nodes = generateMesh(a, b, nNodes);
  std::cout << nodes << '\n';

  // Generate shape functionss
  shapeFunction();

  //  Plot shape functions (only once)
  std::ofstream shapeFile("mpm1D_shapes.txt");
  shapeFile << std::fixed << std::setprecision(6);
  shapeFile << "#x \t N1 \t N2 \t N1_x \t N2_x\n";
  Index nPoints{100};
  for (Index i{0}; i <= nPoints; ++i) {
    double x = a + (b - a) * i / nPoints;
    double N1_val = N[0][0](x);
    double N2_val = N[1][0](x);
    double N1_x_val = B[0][0](x);
    double N2_x_val = B[1][0](x);
    shapeFile << x << " \t " << N1_val << " \t " << N2_val << " \t " << N1_x_val
              << " \t " << N2_x_val << "\n";
  }
  shapeFile.close();
  std::cout << "Shape functions saved to mpm1D_shapes.txt\n";

  // Time integration loop
  for (Index step{0}; step < nSteps; step++) {
    // Save MPM solution
    current_time = step * dt;
    times.push_back(current_time);
    positions_mpm.push_back(x_p);
    velocities_mpm.push_back(v_p);

    // Save analytical solution for comparison
    double x_p_analytical = position_solution(x_loc, current_time);
    double v_p_analytical = velocity_solution(x_loc, current_time);
    positions.push_back(x_p_analytical);
    velocities.push_back(v_p_analytical);

    // Recalculate shape functions at material point's current position N_p[i] =
    // N_i(x_p) (from MPM) (TODO: should check which element contains x_p)
    for (Index i{0}; i < nNodes; i++) {
      N_p[i] = N[i][0](x_p); // only 1 element (number 0)
      B_p[i] = B[i][0](x_p);
    }

    // Reset nodal masses and momentum to zero
    for (Index i{0}; i < nNodes; i++) {
      m_i[i] = mv_i[i] = 0.0;
    }

    // Accumulate mass and momentum from all material points
    // For each particle p, accumulate its contribution to all nodes i
    for (Index p{0}; p < nMPs; p++) {
      // Note: Currently only 1 particle at x_p_current
      // N_p[i] = shape function value of node i at particle p's position
      for (Index i{0}; i < nNodes; i++) {
        m_i[i] += N_p[i] * mass_p;        // m_i = Î£_p N_i(x_p) * mass_p
        mv_i[i] += N_p[i] * mass_p * v_p; // mv_i = Î£_p N_i(x_p) * mass_p * v_p
      }
    }
    // 3.Compute nodal force
    // f_i,ext = b_i + t_i (body force + nodal force)
    f_ext_i.resetZero();
    // f_i,int = V_p*B_i(x_p)*stress_p (body force + nodal force)
    f_int_i.resetZero();
    for (Index p{0}; p < nMPs; p++) {
      for (Index i{0}; i < nNodes; i++) {
        f_int_i[i] += -B_p[i] * V_p *
                      modelParameters::stress_p; // m_i = Î£_p N_i(x_p) * mass_p
      }
    }
    f_total_i = f_int_i + f_ext_i; // Unbalanced at this point

    // 4. Applying Boundary condition
    mv_i[0] = 0;      // v[0] = 0
    f_total_i[0] = 0; // a[0] = 0

    // Update nodal momentum
    mv_i += f_total_i * dt;

    // 5. Node to particle
    // Update particle velocity: v_p^{n+1} = v_p^n + Î”t * Î£_i N_i(x_p) * a_i
    // where a_i = f_i / m_i
    // Update particle position: x_p^{n+1} = x_p^n + Î”t * Î£_i N_i(x_p) *
    // v_i^{n+1} where v_i^{n+1} = mv_i^{n+1} / m_i
    for (Index i{0}; i < nNodes; i++) {
      v_p += N_p[i] * (f_total_i[i] / m_i[i]) * dt;
      x_p += N_p[i] * (mv_i[i] / m_i[i]) * dt;
    }

    // Check for NaN or out of bounds
    if (std::isnan(x_p) || std::isnan(v_p)) {
      std::cerr << "ERROR at step " << step << ": NaN detected!\n";
      std::cerr << "  x_p = " << x_p << ", v_p = " << v_p << "\n";
      std::cerr << "  m_i = " << m_i << "\n";
      std::cerr << "  N_p = " << N_p << "\n";
      break;
    }

    if (x_p < a || x_p > b) {
      std::cerr << "WARNING at step " << step << ": Particle out of domain!\n";
      std::cerr << "  x_p = " << x_p << " (domain: [" << a << ", " << b
                << "])\n";
    }

    // Update particle momentum
    mv_p = mass_p * v_p;
    // Update nodal velocity
    for (Index i{0}; i < nNodes; i++) {
      if (m_i[i] > 1e-12) {
        v_i[i] = mass_p * v_p * N_p[i] / m_i[i];
        v_i[i] = mv_i[i] / m_i[i]; // Causing losing in energy
      } else {
        v_i[i] = 0.0;
      }
    }
    // Apply BC
    v_i[0] = 0;

    // Compute strain rate and stress
    strain_rate_p = 0.0; // Reset strain rate
    for (Index i{0}; i < nNodes; i++) {
      strain_rate_p += B_p[i] * v_i[i]; // Mapping strain rate from the node
    }
    dStrain_p = strain_rate_p * dt;
    stress_p += E * dStrain_p; // Elastic
  }
  std::cout << "Mapped mass to nodes (last step): " << m_i << '\n';

  std::ofstream txtFile("mpm1D.txt");
  txtFile << std::fixed << std::setprecision(6);
  txtFile << "#Time \t x_analytical \t v_analytical \t x_MPM \t v_MPM \t "
             "error_x \t error_v\n";
  for (Index i = 0; i < nSteps; ++i) {
    double error_x = constexpr_fabs(positions[i] - positions_mpm[i]);
    double error_v = constexpr_fabs(velocities[i] - velocities_mpm[i]);
    txtFile << times[i] << " \t " << positions[i] << " \t " << velocities[i]
            << " \t " << positions_mpm[i] << " \t " << velocities_mpm[i]
            << " \t " << error_x << " \t " << error_v << "\n";
  }
  txtFile.close();

  std::cout << "\nComparison at t=0:\n";
  std::cout << "  Position: analytical=" << positions[0]
            << ", MPM=" << positions_mpm[0]
            << ", error=" << constexpr_fabs(positions[0] - positions_mpm[0])
            << "\n";
  std::cout << "  Velocity: analytical=" << velocities[0]
            << ", MPM=" << velocities_mpm[0]
            << ", error=" << constexpr_fabs(velocities[0] - velocities_mpm[0])
            << "\n";
  std::cout << "\nComparison at final time:\n";
  std::cout << "  Position: analytical=" << positions.back()
            << ", MPM=" << positions_mpm.back() << ", error="
            << constexpr_fabs(positions.back() - positions_mpm.back()) << "\n";
  std::cout << "  Velocity: analytical=" << velocities.back()
            << ", MPM=" << velocities_mpm.back() << ", error="
            << constexpr_fabs(velocities.back() - velocities_mpm.back())
            << "\n";

  std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";
  return 0;
}
