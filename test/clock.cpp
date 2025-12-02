#include "../library/clock.h"
#include <iostream>

int main() {
  Timer t;

  // Code to time goes here

  std::cout << "Time elapsed: " << t.elapsed() << " seconds\n";

  return 0;
}
