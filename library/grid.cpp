#include "Matrix.h"
#include "Vector.h"
#include "squareRoot.h"
#include <iostream>
#include <vector>

struct Point {
  double x{};
  double y{};
};

struct Mesh {
  int dimension{};
  double L{}, H{}, B{};
  Index nx{}, ny{}, nz{};
  Index nNodes{}, nElements{};
  Vector<Index> eleOrigin;
  Vector<Index> eleEnd;
  double lx{}, ly{};
  Vector<Point> nodes(nNodes);
  Index first_node{0}; // index of the last node
  void generateMesh(int dimension) {
    switch (dimension) {
    case 1: {
      nNodes = nx;
      nElements = nNodes - 1;
      lx = L / (nx - 1);
      Index last_node{nNodes - 1}; // index of the last node
      for (Index i{0}; i < nNodes; ++i) {
        nodes[i].x = lx * i;
      }

    } break;
    case 2: {
      nNodes = nx * ny;
      nElements = (nx - 1) * (ny - 1);
      lx = L / (nx - 1);
      ly = H / (ny - 1);
      Index last_node{nNodes - 1}; // index of the last node
      for (Index i{0}; i < nNodes; ++i) {
        for (Index j{0}; j < ny; j++) {
        }
        for (Index k{0}; k < nx; k++) {
        }
      }
    } break;
    default:
      std::cout << "Not yet implemented this dimension";
      // Implement 3D
      break;
    }
  }
};

Vector<double> generateMesh(int dimension) {
  for (Index i{0}; i < nNodes; ++i) {
    nodes.push_back((b - a) / (nNodes - 1) * i + a);
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
  }
  return
}
Vector<double> generateMesh(int dimension) {
  using namespace modelParameters;
  for (Index i{0}; i < nNodes; ++i) {
    nodes.push_back((b - a) / (nNodes - 1) * i + a);
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
  }
}

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
  nodes = generateMesh(a, b, nNodes);
  std::cout << "Coordinates of nodes x_i: \n" << nodes << '\n';
  // Generate element connectivity automatically from number of nodes
  eleOrigin = Vector<Index>();
  eleEnd = Vector<Index>();
  for (Index e = 0; e < nElements; ++e) {
    eleOrigin.push_back(e);
    eleEnd.push_back(e + 1);
  }
}