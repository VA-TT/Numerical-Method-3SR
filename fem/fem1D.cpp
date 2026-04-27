#include "../library/FEM.h"
#include "../library/clock.h"
#include <iostream>

int main() {
  std::cout << "=== FEM 1D PROBLEM ===\n\n";

  Timer timer;

  // Create FEM problem: EA=10, L=1, 6 nodes
  // f = 5 kN, F = 10 kN at x=1
  // Expected: R=-15 kN, u(x)=1.5x-0.25x^2, N(x)=15-5x
  constexpr Index nNodes = 6;
  double EA = 10.0;
  double L = 1.0;
  // syntax: FEM1D<type, nNodes> object(EA, L)
  FEM1D<double, nNodes> beam(EA, L);

  // Set distributed load: f(x) = 5 (constant)
  beam.setDistributedLoad(5.0);

  // Set analytical solution: u(x) = 1.5x - 0.25x^2
  beam.setAnalyticSolution([](double x) { return 1.5 * x - 0.25 * x * x; });

  // Assemble system
  beam.assembleKF();

  // Apply boundary conditions
  // IMPORTANT: Apply Neumann BC FIRST, then Dirichlet BC)
  Index first_node = 0;
  Index last_node = beam.getNumNodes() - 1;
  beam.applyNeumanCondition(last_node, 10.0);    // N(L) = F = 10
  beam.applyDirichletCondition(first_node, 0.0); // u(0) = 0

  // Solve
  beam.solveFEM();

  // Post-processing
  beam.calculateReaction();
  beam.calculateAxialForce();
  beam.compareAnalytic();
  beam.exportResult("fem1D_6nodes.txt");

  std::cout << "\n=== Time elapsed: " << timer.elapsed() << " seconds ===\n";
  return 0;
}
