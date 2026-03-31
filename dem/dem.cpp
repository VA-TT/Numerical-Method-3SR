#include "../library/dem.h"
#include "iostram"

int main() {
  constexpr Index ngh = 5;
  constexpr Index ngx = 5;
  double rmin = 0.5e-3;    // smallest diameter is 1 mm
  double rmax = 1e-3;      // largest diameter is 2 mm
  double density = 2700.0; // density of sand particles (2700 kg/m^3)
  double xmax = 0.0;       // right limit of container (width)
  double ymax = 0.0;       // top limit of container (height)
  double ymax0 = 0.0; // initial height of container (for strain estimation)
  double vRand;       // maximum random velocity (beginning of deposit)

  // kinematic of the walls
  double topf = 0.0;
  double topa = 0.0;
  double topv = 0.0;
  double rightf = 0.0;
  double righta = 0.0;
  double rightv = 0.0;
}