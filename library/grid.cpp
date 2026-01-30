#include "Vector.h"
#include <iostream>
#include <vector>
namespace modelParameters {
////////////////////////////////////////////////////////////
///////////////////// Input parameters /////////////////////
////////////////////////////////////////////////////////////
// Problem's domain (a,b)
double a{0.0};
double b{1.0};
// Dirichlet and Neumann boundary values
constexpr double g{0.0};  // Dirichlet at x=0 (u(0)=g)
constexpr double h{10.0}; // Neumann at x=1 (N(1)=h)
// Uniform distributed load . If spatially varying load is wanted, replace
// rhsFunction accordingly (e.g. [](auto x){ return x; }).
constexpr double uniform_load{5.0};
constexpr Index nNodes{6};             // Numbers of nodes
constexpr Index nElements{nNodes - 1}; // Numbers of elements
// Start and end of elements, will be generated automatically in main
Vector<Index> eleOrigin;
Vector<Index> eleEnd;
constexpr double EA = 10.0;

// Right hand side function: default is uniform load `uniform_load`.
auto rhsFunction = [](auto x) { return uniform_load; };

////////////////////////////////////////////////////////////
//////////////////////////// END ///////////////////////////
////////////////////////////////////////////////////////////

// Initiate needed containers
Vector<Vector<std::function<double(double)>>> N(nNodes), N_x(nNodes);
Matrix<double, nNodes, nNodes> K{};
Matrix<double, nNodes, 1> U{};
Matrix<double, nNodes, 1> F{};
Matrix<double, nNodes, 1> R{};
Vector<double> nodes{};
Vector<Vector<double>> element(nElements);
Vector<double> k(nElements);
Vector<double> length(nElements);
Index first_node{0};         // index of the last node
Index last_node{nNodes - 1}; // index of the last node

} // namespace modelParameters

Vector<double> generateMesh(double a, double b, Index n) {
  Vector<double> nodes;
  for (Index i{0}; i < n; ++i) {
    nodes.push_back((b - a) / (n - 1) * i + a);
  }
  return nodes;
}

// calculating shape functions and its derivative
// (obviously redundance but objectif here is to show the algorithm
// in non-shorcut form and test out implemented libraries)
void shapeFunction() {
  // Equally divied
  using namespace modelParameters;
  int local_i = 0;
  int local_j = 1;
  for (Index e{0}; e < nElements; ++e) {
    int i = eleOrigin[e];
    int j = eleEnd[e];
    double x_i = nodes[i];
    double x_j = nodes[j];
    element[e] = {x_i, x_j};
    length[e] = constexpr_fabs(x_j - x_i);
    k[e] = EA / length[e];

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
        N_x[kk].push_back(dShape_i);
      } else if (kk == j) {
        N[kk].push_back(shapefunction_j);
        N_x[kk].push_back(dShape_j);
      } else {
        N[kk].push_back([](double x) { return 0.0; });
        N_x[kk].push_back([](double x) { return 0.0; });
      }
    }
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