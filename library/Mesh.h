#include "Vector.h"
#include "parentElement.h"
#include <cassert>
#include <iostream>

// 1D Mesh for line elements
template <typename T> class Mesh1D {
private:
  T m_length{};
  Index m_nNodes{}, m_nElements{};
  Vector<T> m_nodes{};
  Vector<Vector<Index>> m_connectivity{}; // [node_i, node_j]
  Vector<Vector<T>> m_elements{};         // [xi, xj]

public:
  // Constructor
  Mesh1D(T length, Index nNodes)
      : m_length{length}, m_nNodes{nNodes}, m_nElements{nNodes - 1} {
    assert(m_nElements > 0 && "Number of elements must be positive");
    assert(length > 0 && "Domain length must be positive");

    // Resize vectors before using them
    m_nodes.resize(m_nNodes);
    m_connectivity.resize(m_nElements);
    m_elements.resize(m_nElements);

    T lx{m_length / m_nElements};
    for (Index i{0}; i < m_nNodes; ++i) {
      m_nodes[i] = lx * i; // Coordinates of nodes
    }
    for (Index e{0}; e < m_nElements; ++e) {
      m_connectivity[e] = {e, e + 1};               // Node-index of elements
      m_elements[e] = {m_nodes[e], m_nodes[e + 1]}; // Coordinates of elements
    }
  };
  // Other defaults
  Mesh1D() = default;
  Mesh1D(const Mesh1D &) = default;
  Mesh1D(Mesh1D &&) = default;
  Mesh1D &operator=(const Mesh1D &) = default;
  Mesh1D &operator=(Mesh1D &&) = default;
  ~Mesh1D() = default;

  // Getter
  Index getNumNodes() const { return m_nNodes; }
  Index getNumElements() const { return m_nElements; }
  const Vector<T> &nodeCoords() const { return m_nodes; }
  T getLengthEle(Index e) const {
    assert(e >= 0 && e < m_nElements && "Invalid element ID");
    return std::abs(m_elements[e][0] - m_elements[e][1]);
  };
  std::pair<T, T> getElementNodes(Index e) const {
    assert(e >= 0 && e < m_nElements && "Invalid element ID");
    return {m_nodes[e], m_nodes[e + 1]};
  }
  // Print mesh info
  void print() const {
    std::cout << "=== 1D Mesh Info ===" << '\n';
    std::cout << "Number of nodes: " << m_nNodes << '\n';
    std::cout << "Number of elements: " << m_nElements << '\n';
    std::cout << "Node coordinates: " << m_nodes << '\n';
    std::cout << "Connectivity:\n";
    for (Index e = 0; e < m_nElements; ++e) {
      std::cout << "  Element " << e << ": [" << m_connectivity[e][0] << ", "
                << m_connectivity[e][1] << "]\n";
    }
  }
  // Setter
  void setNumElements(Index number) {
    m_nElements = number;
    m_nNodes = number + 1;
  }
  void setNumNodes(Index number) {
    m_nNodes = number;
    m_nElements = number - 1;
  }
  void setLength(T length) { m_length = length; }
  void regenerateMesh() { *this = Mesh1D(m_length, m_nNodes); }

  // MPM helper function
  Index findCageID(T x) const {
    for (Index e{0}; e < m_nElements; ++e) {
      if (m_elements[e][0] <= x && m_elements[e][1] >= x)
        return e;
    }
    return -1; // Not found (outside domain)
  }
};

template <typename T> class Mesh2D {
private:
  T m_length{}, m_height{};
  Index m_nx{}, m_ny{}, m_nNodes{}, m_nElements{};
  Vector<std::pair<T, T>> m_nodes{}; // All node coordinates as (x, y) pairs
  Vector<Vector<Index>> m_connectivity{};
  Vector<T> m_x_coords; // Node x-coordinates vector
  Vector<T> m_y_coords; // Node y-coordinates vector

public:
  // Constructor
  Mesh2D(T length, T height, Index nx, Index ny)
      : m_length{length}, m_height{height}, m_nx{nx}, m_ny{ny},
        m_nNodes{nx * ny}, m_nElements{(nx - 1) * (ny - 1)} {
    assert(nx > 1 && ny > 1 && "Need at least 2 nodes in each direction");
    assert(length > 0 && height > 0 && "Domain dimensions must be positive");

    // Resize vectors before using them
    m_nodes.resize(m_nNodes);
    m_connectivity.resize(m_nElements);
    m_x_coords.resize(m_nNodes);
    m_y_coords.resize(m_nNodes);

    T dx = length / (nx - 1);
    T dy = height / (ny - 1);
    // Node numbering: row-major (i + j*nx)
    // y ^
    //   |  6--7--8
    //   |  3--4--5
    //   |  0--1--2
    //   +---------> x
    Index nodeID{0};
    for (Index j{0}; j < ny; j++) {
      for (Index i{0}; i < nx; i++) {
        nodeID = j * nx + i;
        T x = dx * i;
        T y = dy * j;
        m_x_coords[nodeID] = x;
        m_y_coords[nodeID] = y;
        m_nodes[nodeID] = {x, y}; // Store as pair
      }
    }

    // Generate connectivity for Q4 elements
    // Element node ordering (counterclockwise from bottom-left):
    //   n4 ------ n3
    //   |          |
    //   |          |
    //   n1 ------ n2
    Index elemID{0};
    Vector<Index> eleConnectivity(4);
    for (Index j{0}; j < ny - 1; j++) {
      for (Index i{0}; i < nx - 1; i++) {
        eleConnectivity[0] = i + j * nx;
        eleConnectivity[1] = (i + 1) + j * nx;
        eleConnectivity[2] = (i + 1) + (j + 1) * nx;
        eleConnectivity[3] = i + (j + 1) * nx;
        m_connectivity[elemID] = eleConnectivity;
        ++elemID;
      }
    }
  };
  // Other defaults
  Mesh2D() = default;
  Mesh2D(const Mesh2D &) = default;
  Mesh2D(Mesh2D &&) = default;
  Mesh2D &operator=(const Mesh2D &) = default;
  Mesh2D &operator=(Mesh2D &&) = default;
  ~Mesh2D() = default;

  // Getter
  Index getNumNodes() const { return m_nNodes; }
  Index getNumElements() const { return m_nElements; }
  Index nx() const { return m_nx; }
  Index ny() const { return m_ny; }
  const Vector<std::pair<T, T>> &getNodes() const { return m_nodes; }
  const Vector<T> &getXCoords() const { return m_x_coords; }
  const Vector<T> &getYCoords() const { return m_y_coords; }
  // Get element connectivity
  const Vector<Index> &getConnectivity(Index elemID) const {
    assert(elemID >= 0 && elemID < m_nElements && "Invalid element ID");
    return m_connectivity[elemID];
  }

  // Get node coordinates by ID
  std::pair<T, T> getNode(Index nodeID) const {
    assert(nodeID >= 0 && nodeID < m_nNodes && "Invalid node ID");
    return {m_x_coords[nodeID], m_y_coords[nodeID]};
  }

  // Print mesh info
  void print() const {
    std::cout << "=== 2D Mesh Info ===" << '\n';
    std::cout << "Number of nodes: " << m_nNodes << " (" << m_nx << " x "
              << m_ny << ")\n";
    std::cout << "Number of elements: " << m_nElements << '\n';

    std::cout << "\nFirst 5 nodes:\n";
    for (Index i = 0; i < std::min(m_nNodes, Index(5)); ++i) {
      std::cout << "  Node " << i << ": (" << m_x_coords[i] << ", "
                << m_y_coords[i] << ")\n";
    }

    std::cout << "\nFirst 3 elements:\n";
    for (Index e = 0; e < std::min(m_nElements, Index(3)); ++e) {
      std::cout << "  Element " << e << ": [";
      for (Index i = 0; i < 4; ++i) {
        std::cout << m_connectivity[e][i];
        if (i < 3)
          std::cout << ", ";
      }
      std::cout << "]\n";
    }
  }
  // Setter
  void setXnodes(Index nx) { m_nx = nx; }
  void setYnodes(Index ny) { m_ny = ny; }
  void setLength(T length) { m_length = length; }
  void setHeight(T height) { m_height = height; }
  void regenerateMesh() { *this = Mesh2D(m_length, m_height, m_nx, m_ny); }

  // MPM helper function
  bool isPointInElement(T x, T y, Index e) const {
    assert(e >= 0 && e < m_nElements && "Invalid element ID");
    const auto &localNodes = m_connectivity[e];
    Vector<T> x_nodes(4), y_nodes(4);
    for (Index i{0}; i < 4; i++) {
      x_nodes[i] = m_x_coords[localNodes[i]];
      y_nodes[i] = m_y_coords[localNodes[i]];
    }
    auto [xi, eta] = parentCoor2D(x, y, x_nodes, y_nodes);
    return (xi >= -1.0 && xi <= 1.0 && eta >= -1.0 && eta <= 1.0);
  }
  Index findCageID(T x, T y, Index lastElement = -1) const {
    // Optimized search: start from last known element
    if (lastElement >= 0 && lastElement < m_nElements) {
      if (isPointInElement(x, y, lastElement)) {
        return lastElement;
      }
    }
    for (Index e = 0; e < m_nElements; ++e) {
      if (isPointInElement(x, y, e)) {
        return e;
      }
    }
    return -1; // Not found (outside domain)
  }
};
