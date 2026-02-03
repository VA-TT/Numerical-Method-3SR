#include "../library/FEM.h"
#include "../library/clock.h"
#include <iostream>

int main() {
  std::cout << "=== FEM 1D PROBLEM ===\n\n";

  Timer timer;

  // Create FEM problem: EA=10, L=1, 6 nodes
  // f = 5 kN, F = 10 kN at x=1
  // Expected: R=-15 kN, u(x)=1.5x-0.25x^2, N(x)=15-5x
  FEM1D<double, 6> beam(10.0, 1.0);

  // Set distributed load: f(x) = 5 (constant)
  beam.setDistributedLoad(5.0);

  // Set analytical solution: u(x) = 1.5x - 0.25x^2
  beam.setAnalyticSolution([](double x) { return 1.5 * x - 0.25 * x * x; });

  // Assemble system
  beam.assembleKF();

  // Apply boundary conditions
  Index first_node = 0;
  Index last_node = beam.getNumNodes() - 1;
  beam.applyDirichletCondition(first_node, 0.0); // u(0) = 0
  beam.applyNeumanCondition(last_node, 10.0);    // EA*u'(L) = 10

  // Solve
  beam.solveFEM();

  // Post-processing
  beam.calculateReaction();
  beam.compareAnalytic();
  beam.exportResult("fem1D_6nodes.txt");

  std::cout << "\n=== Time elapsed: " << timer.elapsed() << " seconds ===\n";
  return 0;
}
