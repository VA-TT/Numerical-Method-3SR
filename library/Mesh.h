#ifndef DISCRETIZING_RECTANGULAR_MESH_H
#define DISCRETIZING_RECTANGULAR_MESH_H

#include "Vector.h"
#include "parentElement.h"
#include <cassert>
#include <iostream>

// 1D Mesh for line elements
template <typename T> class Mesh1D {
private:
  T m_length{};
  Index m_nNodes{}, m_nElements{}, m_nMPperEle{}, m_nMPs{};
  Vector<T> m_nodes{};
  Vector<char> m_activeNodes{}; // Use char instead of bool to avoid
                                // std::vector<bool> issues
  Vector<T> m_MPs{};
  Vector<Index> m_mpElementId{}; // Cached element ID for each MP
  Vector<T> m_nodes_initial{};   // Store initial configuration for reset
  Vector<Vector<Index>> m_connectivity{}; // [node_i, node_j]

public:
  // Constructor
  Mesh1D(T length, Index nNodes, Index nMPperEle = 0)
      : m_length{length}, m_nNodes{nNodes}, m_nElements{nNodes - 1},
        m_nMPperEle{nMPperEle}, m_nMPs{nMPperEle * (nNodes - 1)} {
    assert(m_nElements > 0 && "Number of elements must be positive");
    assert(length > 0 && "Domain length must be positive");

    // Resize vectors before using them
    m_nodes.resize(m_nNodes);
    m_nodes_initial.resize(m_nNodes);
    m_activeNodes.resize(m_nNodes);
    m_connectivity.resize(m_nElements);

    T lx{m_length / m_nElements};
    for (Index i{0}; i < m_nNodes; ++i) {
      m_nodes[i] = lx * i; // Coordinates of nodes
    }
    // Copy to initial configuration
    m_nodes_initial = m_nodes;

    // Initialize active node flags to 0 (inactive)
    for (Index i{0}; i < m_nNodes; ++i) {
      m_activeNodes[i] = 0;
    }

    for (Index e{0}; e < m_nElements; ++e) {
      m_connectivity[e] = {e, e + 1}; // Node indices of elements
    }

    // Initialize Material Points (MPs) if needed
    if (m_nMPperEle > 0) {
      m_MPs.reserve(m_nMPs);
      for (Index e{0}; e < m_nElements; ++e) {
        T x_start = m_nodes[e];
        T le = getLengthEle(e);
        for (Index p{0}; p < m_nMPperEle; ++p) {
          m_MPs.push_back(x_start + (p + 1) * le / (m_nMPperEle + 1));
        }
      }
    }

    // Initialize cached element IDs for MPs
    m_mpElementId.resize(static_cast<std::size_t>(m_nMPs), Index{-1});
    for (Index p{0}; p < static_cast<Index>(m_MPs.size()); ++p) {
      m_mpElementId[p] = findCageID(m_MPs[p]);
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
  Index getNumMPs() const { return m_nMPs; }
  const Vector<char> &getActiveNodes() const { return m_activeNodes; }
  const Vector<T> &nodeCoords() const { return m_nodes; }
  const Vector<T> &getMPCoords() const { return m_MPs; }
  const Vector<Index> &getMPElementIds() const { return m_mpElementId; }
  T getMPCoord(Index p) const {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    return m_MPs[p];
  }
  Index getMPElementId(Index p) const {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    return m_mpElementId[p];
  }

  void setMPCoords(const Vector<T> &mp_positions) {
    assert(mp_positions.size() == static_cast<std::size_t>(m_nMPs) &&
           "Size mismatch: mp_positions must match mesh MP count");
    m_MPs = mp_positions;
  }

  void setMPCoord(Index p, T x) {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    m_MPs[p] = x;
    // Keep cached element id consistent
    if (m_mpElementId.size() != m_nMPs) {
      m_mpElementId.resize(static_cast<std::size_t>(m_nMPs), Index{-1});
    }
    m_mpElementId[p] = findCageID(x);
  }

  // Update cached element IDs using the mesh-internal MP coordinates.
  void updateMPElementIds() {
    if (m_mpElementId.size() != m_nMPs) {
      m_mpElementId.resize(static_cast<std::size_t>(m_nMPs), Index{-1});
    }
    for (Index p{0}; p < m_nMPs; ++p) {
      m_mpElementId[p] = findCageID(m_MPs[p]);
    }
  }
  const Vector<Vector<Index>> &getConnectivity() const {
    return m_connectivity;
  }

  // Get connectivity of a single element
  const Vector<Index> &getEleConnectivity(Index e) const {
    assert(e >= 0 && e < m_nElements && "Invalid element ID");
    return m_connectivity[e];
  }

  T getLengthEle(Index e) const {
    assert(e >= 0 && e < m_nElements && "Invalid element ID");
    Index n1 = m_connectivity[e][0];
    Index n2 = m_connectivity[e][1];
    return std::abs(m_nodes[n2] - m_nodes[n1]);
  }

  std::pair<T, T> getElementNodes(Index e) const {
    assert(e >= 0 && e < m_nElements && "Invalid element ID");
    Index n1 = m_connectivity[e][0];
    Index n2 = m_connectivity[e][1];
    return {m_nodes[n1], m_nodes[n2]};
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
  void regenerateMesh() { *this = Mesh1D(m_length, m_nNodes, m_nMPperEle); }
  void setNumElements(Index number) {
    m_nElements = number;
    m_nNodes = number + 1;
    regenerateMesh();
  }
  void setNumNodes(Index number) {
    m_nNodes = number;
    m_nElements = number - 1;
    regenerateMesh();
  }
  void setNumMPs(Index number) {
    m_nMPperEle = number;
    regenerateMesh();
  }
  void setLength(T length) { m_length = length; }

  // Update mesh for time-stepping (e.g., Updated Lagrangian FEM)
  void updateNodePosition(Index nodeID, T new_x) {
    assert(nodeID >= 0 && nodeID < m_nNodes && "Invalid node ID");
    m_nodes[nodeID] = new_x;
  }

  void updateAllNodes(const Vector<T> &new_positions) {
    assert(new_positions.size() == m_nNodes && "Size mismatch");
    m_nodes = new_positions;
  }

  // Reset to initial configuration
  void nodalReset() { m_nodes = m_nodes_initial; } // Excluding MPs

  const Vector<T> &getInitialNodes() const { return m_nodes_initial; }

  // MPM helper function: Find element containing position x
  Index findCageID(T x) const {
    for (Index e{0}; e < m_nElements; ++e) {
      Index n1 = m_connectivity[e][0];
      Index n2 = m_connectivity[e][1];
      T x1 = m_nodes[n1];
      T x2 = m_nodes[n2];
      if (x1 <= x && x2 >= x)
        return e;
    }
    return -1; // Not found (outside domain)
  }

  // Activate nodes that contain Material Points
  void activateNodes() {
    // Defensive: ensure storage is correctly sized
    if (m_activeNodes.size() != m_nNodes) {
      m_activeNodes.resize(static_cast<std::size_t>(m_nNodes));
    }
    // Reset all nodes to inactive
    for (Index i{0}; i < m_nNodes; ++i) {
      m_activeNodes[i] = 0;
    }

    // Activate nodes of elements containing MPs
    for (Index p{0}; p < m_nMPs; ++p) {
      Index elemID = m_mpElementId[p];
      if (elemID != -1) {
        // Get node indices (not coordinates!)
        Index n1 = m_connectivity[elemID][0];
        Index n2 = m_connectivity[elemID][1];
        m_activeNodes[n1] = 1;
        m_activeNodes[n2] = 1;
      }
    }
  }

  // Check if a node is active (contains MPs)
  bool isActiveNode(Index nodeID) const {
    assert(nodeID >= 0 && nodeID < m_nNodes && "Invalid node ID");
    return m_activeNodes[nodeID] != 0;
  }
};

template <typename T> class Mesh2D {
private:
  T m_length{}, m_height{};
  Index m_nx{}, m_ny{}, m_nNodes{}, m_nElements{}, m_nMPs{};
  Vector<std::pair<T, T>> m_nodes{}; // All node coordinates as (x, y) pairs
  Vector<std::pair<T, T>> m_nodes_initial{}; // Initial configuration for reset
  Vector<std::pair<T, T>> m_MPs{};           // Material Points coordinates
  Vector<char> m_activeNodes{}; // Track which nodes contain Material Points
  Vector<Vector<Index>> m_connectivity{};
  Vector<T> m_nodes_x;         // Node x-coordinates vector
  Vector<T> m_nodes_y;         // Node y-coordinates vector
  Vector<T> m_nodes_x_initial; // Initial nodal x-coordinates
  Vector<T> m_nodes_y_initial; // Initial nodal y-coordinates
  Vector<T> m_MP_x_initial;    // Initial MP x-coordinates
  Vector<T> m_MP_y_initial;    // Initial MP y-coordinates
  Vector<T> m_MP_x;            // MP x-coordinates
  Vector<T> m_MP_y;            // MP y-coordinates

public:
  // Constructor
  Mesh2D(T length, T height, Index nx, Index ny, T x0 = 0.0, T y0 = 0.0,
         T x1 = 0.0, T y1 = 0.0, T MP_size = 0.0)
      : m_length{length}, m_height{height}, m_nx{nx}, m_ny{ny},
        m_nNodes{nx * ny}, m_nElements{(nx - 1) * (ny - 1)} {
    assert(nx > 1 && ny > 1 && "Need at least 2 nodes in each direction");
    assert(length > 0 && height > 0 && "Domain dimensions must be positive");

    // Resize vectors before using them
    m_nodes.resize(m_nNodes);
    m_nodes_initial.resize(m_nNodes);
    m_activeNodes.resize(m_nNodes, 0);
    m_connectivity.resize(m_nElements);
    m_nodes_x.resize(m_nNodes);
    m_nodes_y.resize(m_nNodes);
    m_nodes_x_initial.resize(m_nNodes);
    m_nodes_y_initial.resize(m_nNodes);

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
        m_nodes_x[nodeID] = x;
        m_nodes_y[nodeID] = y;
        m_nodes[nodeID] = {x, y};
      }
    }
    // Copy to initial configuration (more efficient than per-element
    // assignment)
    m_nodes_x_initial = m_nodes_x;
    m_nodes_y_initial = m_nodes_y;
    m_nodes_initial = m_nodes;

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
    generateSquareGridMP(x0, y0, x1, y1, MP_size);
  } // End of Mesh2D constructor

  // MP grid generation methods
  void generateSquareGridMP(T x0, T y0, T x1, T y1, T MP_size) {
    Index npx = static_cast<Index>((x1 - x0) / MP_size);
    Index npy = static_cast<Index>((y1 - y0) / MP_size);
    m_nMPs = npx * npy;

    // Initialize Material Points (MPs) if needed
    if (m_nMPs > 0) {
      m_MPs.reserve(m_nMPs);
      m_MP_x.reserve(m_nMPs);
      m_MP_y.reserve(m_nMPs);
      for (Index j{0}; j < npy; ++j) {
        for (Index i{0}; i < npx; ++i) {
          T x = x0 + MP_size / 2 + i * MP_size;
          T y = y0 + MP_size / 2 + j * MP_size;
          m_MP_x.push_back(x);
          m_MP_y.push_back(y);
          m_MPs.push_back({x, y});
        }
      }
      m_MP_x_initial = m_MP_x;
      m_MP_y_initial = m_MP_y;
    }
  }
  void generateCircleMPgrid() {}
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
  Index getNumMPs() const { return m_nMPs; }
  const Vector<char> &getActiveNodes() const { return m_activeNodes; }
  Index nx() const { return m_nx; }
  Index ny() const { return m_ny; }
  const Vector<std::pair<T, T>> &getNodes() const { return m_nodes; }
  const Vector<std::pair<T, T>> &getMPCoords() const { return m_MPs; }
  std::pair<T, T> getMPCoord(Index p) const {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    return m_MPs[p];
  }
  const Vector<T> &getXCoords() const { return m_nodes_x; }
  const Vector<T> &getYCoords() const { return m_nodes_y; }

  // MP getters
  const Vector<T> &getMPXCoords() const { return m_MP_x; }
  const Vector<T> &getMPYCoords() const { return m_MP_y; }
  const Vector<T> &getInitialMPXCoords() const { return m_MP_x_initial; }
  const Vector<T> &getInitialMPYCoords() const { return m_MP_y_initial; }
  T getMPXCoord(Index p) const {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    return m_MP_x[p];
  }
  T getMPYCoord(Index p) const {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    return m_MP_y[p];
  }

  // MP setters
  void setMPCoord(Index p, T x, T y) {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    m_MP_x[p] = x;
    m_MP_y[p] = y;
    m_MPs[p] = {x, y};
  }

  void setMPCoords(const Vector<T> &mp_x, const Vector<T> &mp_y) {
    assert(mp_x.size() == static_cast<std::size_t>(m_nMPs) &&
           mp_y.size() == static_cast<std::size_t>(m_nMPs) &&
           "Size mismatch: mp coordinates must match mesh MP count");
    m_MP_x = mp_x;
    m_MP_y = mp_y;
    for (Index p{0}; p < m_nMPs; ++p) {
      m_MPs[p] = {mp_x[p], mp_y[p]};
    }
  }

  void setMPCoords(const Vector<std::pair<T, T>> &mp_positions) {
    assert(mp_positions.size() == static_cast<std::size_t>(m_nMPs) &&
           "Size mismatch: mp_positions must match mesh MP count");
    m_MPs = mp_positions;
    for (Index p{0}; p < m_nMPs; ++p) {
      m_MP_x[p] = mp_positions[p].first;
      m_MP_y[p] = mp_positions[p].second;
    }
  }

  void updateMPCoord(Index p, T x, T y) {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    m_MP_x[p] = x;
    m_MP_y[p] = y;
    m_MPs[p] = {x, y};
  }

  void updateAllMPs(const Vector<T> &new_mp_x, const Vector<T> &new_mp_y) {
    assert(new_mp_x.size() == m_nMPs && new_mp_y.size() == m_nMPs &&
           "Size mismatch");
    m_MP_x = new_mp_x;
    m_MP_y = new_mp_y;
    for (Index p{0}; p < m_nMPs; ++p) {
      m_MPs[p] = {new_mp_x[p], new_mp_y[p]};
    }
  }

  void resetMPs() {
    m_MP_x = m_MP_x_initial;
    m_MP_y = m_MP_y_initial;
    for (Index p{0}; p < m_nMPs; ++p) {
      m_MPs[p] = {m_MP_x_initial[p], m_MP_y_initial[p]};
    }
  }

  // Get element connectivity
  const Vector<Index> &getEleConnectivity(Index elemID) const {
    assert(elemID >= 0 && elemID < m_nElements && "Invalid element ID");
    return m_connectivity[elemID];
  }

  // Get node coordinates by ID
  std::pair<T, T> getNodeCoor(Index nodeID) const {
    assert(nodeID >= 0 && nodeID < m_nNodes && "Invalid node ID");
    return {m_nodes_x[nodeID], m_nodes_y[nodeID]};
  }

  // Print mesh info
  void print() const {
    std::cout << "=== 2D Mesh Info ===" << '\n';
    std::cout << "Number of nodes: " << m_nNodes << " (" << m_nx << " x "
              << m_ny << ")\n";
    std::cout << "Number of elements: " << m_nElements << '\n';

    std::cout << "\nFirst 5 nodes:\n";
    for (Index i = 0; i < std::min(m_nNodes, Index(5)); ++i) {
      std::cout << "  Node " << i << ": (" << m_nodes_x[i] << ", "
                << m_nodes_y[i] << ")\n";
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
  void setXnodes(Index nx) {
    m_nx = nx;
    regenerateMesh();
  }
  void setYnodes(Index ny) {
    m_ny = ny;
    regenerateMesh();
  }
  void setLength(T length) {
    m_length = length;
    regenerateMesh();
  }
  void setHeight(T height) {
    m_height = height;
    regenerateMesh();
  }
  void regenerateMesh() { *this = Mesh2D(m_length, m_height, m_nx, m_ny); }

  // Update mesh for time-stepping (e.g., Updated Lagrangian FEM)
  void updateNodePosition(Index nodeID, T new_x, T new_y) {
    assert(nodeID >= 0 && nodeID < m_nNodes && "Invalid node ID");
    m_nodes_x[nodeID] = new_x;
    m_nodes_y[nodeID] = new_y;
    m_nodes[nodeID] = {new_x, new_y};
  }

  void updateAllNodes(const Vector<T> &new_x, const Vector<T> &new_y) {
    assert(new_x.size() == m_nNodes && new_y.size() == m_nNodes &&
           "Size mismatch");
    m_nodes_x = new_x;
    m_nodes_y = new_y;
    for (Index i{0}; i < m_nNodes; ++i) {
      m_nodes[i] = {new_x[i], new_y[i]};
    }
  }

  void updateAllNodes(const Vector<std::pair<T, T>> &new_positions) {
    assert(new_positions.size() == m_nNodes && "Size mismatch");
    m_nodes = new_positions;
    for (Index i{0}; i < m_nNodes; ++i) {
      m_nodes_x[i] = new_positions[i].first;
      m_nodes_y[i] = new_positions[i].second;
    }
  }
  void createMeshMP() {}
  // Reset to initial configuration
  void resetMesh() {
    m_nodes = m_nodes_initial;
    m_nodes_x = m_nodes_x_initial;
    m_nodes_y = m_nodes_y_initial;
  }

  const Vector<std::pair<T, T>> &getInitialNodes() const {
    return m_nodes_initial;
  }
  const Vector<T> &getInitialXCoords() const { return m_nodes_x_initial; }
  const Vector<T> &getInitialYCoords() const { return m_nodes_y_initial; }

  // MPM helper function
  bool isPointInElement(T x, T y, Index e) const {
    assert(e >= 0 && e < m_nElements && "Invalid element ID");
    const auto &localNodes = m_connectivity[e];
    Vector<T> x_nodes(4), y_nodes(4);
    for (Index i{0}; i < 4; i++) {
      x_nodes[i] = m_nodes_x[localNodes[i]];
      y_nodes[i] = m_nodes_y[localNodes[i]];
    }
    auto [xi, eta] = parentCoor(x, y, x_nodes, y_nodes);
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

  //   Activate nodes that contain Material Points
  void activateNodes() {
    for (Index i{0}; i < m_nNodes; ++i) {
      m_activeNodes[i] = 0;
    }
    for (const auto &mp : m_MPs) {
      Index elemID = findCageID(mp.first, mp.second);
      if (elemID != -1) {
        const auto &conn = m_connectivity[elemID];
        for (Index i{0}; i < 4; ++i) {
          m_activeNodes[conn[i]] = 1;
        }
      }
    }
  }

  bool isActiveNode(Index nodeID) const {
    assert(nodeID >= 0 && nodeID < m_nNodes && "Invalid node ID");
    return m_activeNodes[nodeID] != 0;
  }
};
//////// WHAT HAPPENS IF THE MPs IS PERFECTLY LANDS ON THE NODE? ////////

#endif