#include "Vector.h"
#include <iostream>
#include <vector>
int main() {
  int nx = 5;
  int ny = 3;
  double B = 5;
  double H = 3;
  double lx = B / nx;
  double ly = H / ny;
  Vector<Index> eleOrigin;
  Vector<Index> eleEnd;
  Vector<double> nodes{};
  Vector<Vector<double>> element(nElements);
  Vector<double> k(nElements);
  Vector<double> length(nElements);
  Index first_node{0};         // index of the last node
  Index last_node{nNodes - 1}; // index of the last node
}