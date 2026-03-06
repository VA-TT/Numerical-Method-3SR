#ifndef DISCRETIZING_RECTANGULAR_MESH_H
#define DISCRETIZING_RECTANGULAR_MESH_H

#include "Vector.h"
#include "comparison.h"
#include "parentElement.h"
#include <cassert>
#include <iostream>
#include <utility>

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
  Vector<Vector<Index>>
      m_connectivity{}; // [node_i, node_j]:representing element
  static constexpr Index idError = -1;

public:
  // Constructor
  Mesh1D(T length, Index nNodes, Index nMPperEle = 0)
      : m_length{length}, m_nNodes{nNodes}, m_nElements{nNodes - 1},
        m_nMPperEle{nMPperEle}, m_nMPs{nMPperEle * (nNodes - 1)} {
    assert(m_nElements > 0 && "Number of elements must be positive");
    assert(length > 0 && "Domain length must be positive");

    // Initiate vectors to 0
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

    for (Index e{0}; e < m_nElements; ++e) {
      m_connectivity[e] = {e, e + 1}; // Node index of each element
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
    m_mpElementId.resize(m_nMPs, idError);
    for (Index p{0}; p < m_MPs.size(); ++p) {
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
  const Vector<Index> &getMPelementIDs() const { return m_mpElementId; }
  const Vector<T> &getInitialNodes() const { return m_nodes_initial; }

  T getMPCoord(Index p) const {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    return m_MPs[p];
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
  const Vector<T> &getInitialNodes() const { return m_nodes_initial; }

  std::pair<T, T> getElementNodes(Index e) const {
    assert(e >= 0 && e < m_nElements && "Invalid element ID");
    Index n1 = m_connectivity[e][0];
    Index n2 = m_connectivity[e][1];
    return {m_nodes[n1], m_nodes[n2]};
  }

  Index getMPelementID(Index p) const {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    return m_mpElementId[p];
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
  void setLength(T length) {
    m_length = length;
    regenerateMesh();
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
    m_mpElementId[p] = findCageID(x);
  }

  // Update cached element IDs using the mesh-internal MP coordinates.
  void updateMPElementIds() {
    for (Index p{0}; p < m_nMPs; ++p) {
      m_mpElementId[p] = findCageID(m_MPs[p]);
    }
  }

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
    return idError; // Not found (outside domain)
  }

  // List of activated nodes that contain Material Points
  void activateNodes() {
    // Reset all nodes to inactive
    m_activeNodes.resetZero();
    // Activate nodes of elements containing MPs
    for (Index p{0}; p < m_nMPs; ++p) {
      Index elemID = m_mpElementId[p];
      if (elemID != idError) {
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
};

template <typename T> class Mesh2D {
private:
  T m_length{}, m_height{};
  T m_MP_size{}; // MP spacing / particle size used for square-grid MPs
  Index m_nx{}, m_ny{}, m_nNodes{}, m_nElements{}, m_nMPs{};
  Vector<std::pair<T, T>> m_nodes{}; // All node coordinates as (x, y) pairs
  Vector<std::pair<T, T>> m_nodes_initial{}; // Initial configuration for reset
  Vector<std::pair<T, T>> m_MPs{};           // Material Points coordinates
  Vector<char> m_activeNodes{};  // Track which nodes contain Material Points
  Vector<T> m_nodes_x;           // Node x-coordinates vector
  Vector<T> m_nodes_y;           // Node y-coordinates vector
  Vector<T> m_nodes_x_initial;   // Initial nodal x-coordinates
  Vector<T> m_nodes_y_initial;   // Initial nodal y-coordinates
  Vector<T> m_MP_x_initial;      // Initial MP x-coordinates
  Vector<T> m_MP_y_initial;      // Initial MP y-coordinates
  Vector<T> m_MP_x;              // MP x-coordinates
  Vector<T> m_MP_y;              // MP y-coordinates
  Vector<Index> m_mpElementId{}; // Cached element ID for each MP
  Vector<Vector<Index>> m_connectivity{}; // Elements representing
  static constexpr Index idError = -1;

public:
  // Constructor
  Mesh2D(T length, T height, Index nx, Index ny, T x0 = 0.0, T y0 = 0.0,
         T x1 = 0.0, T y1 = 0.0, T MP_size = 0.0)
      : m_length{length}, m_height{height}, m_nx{nx}, m_ny{ny},
        m_nNodes{nx * ny}, m_nElements{(nx - 1) * (ny - 1)},
        m_MP_size{MP_size} {
    assert(nx > 1 && ny > 1 && "Need at least 2 nodes in each direction");
    assert(length > 0 && height > 0 && "Domain dimensions must be positive");

    // Resize vectors before using them
    m_nodes.resize(m_nNodes);
    m_nodes_initial.resize(m_nNodes);
    m_activeNodes.resize(m_nNodes);
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

  // Convenience constructor: pass MP region as corners
  Mesh2D(const std::pair<T, T> &gridDimension,
         const std::pair<Index, Index> &gridnPoints,
         const std::pair<T, T> &minCorner, const std::pair<T, T> &maxCorner,
         T MP_size = 0.0)
      : Mesh2D(gridDimension.first, gridDimension.second, gridnPoints.first,
               gridnPoints.second, minCorner.first, minCorner.second,
               maxCorner.first, maxCorner.second, MP_size) {}

  // MP grid generation methods
  void generateSquareGridMP(T x0, T y0, T x1, T y1, T MP_size) {
    // Defensive: avoid division by 0 / invalid ranges
    if (!(MP_size > T{0}) || !(x1 > x0) || !(y1 > y0)) {
      m_nMPs = 0;
      m_MPs = Vector<std::pair<T, T>>{};
      m_MP_x = Vector<T>{};
      m_MP_y = Vector<T>{};
      m_MP_x_initial = Vector<T>{};
      m_MP_y_initial = Vector<T>{};
      m_mpElementId = Vector<Index>{};
      return;
    }

    // Regenerate MPs: clear previous state first
    m_MPs = Vector<std::pair<T, T>>{};
    m_MP_x = Vector<T>{};
    m_MP_y = Vector<T>{};

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

    // Initialize cached element IDs for MPs
    m_mpElementId.resize(m_nMPs, idError);
    for (Index p{0}; p < m_nMPs; ++p) {
      m_mpElementId[p] = findCageID(m_MP_x[p], m_MP_y[p], m_mpElementId[p]);
    }
  }

  // Convenience overload: pass MP region as corners
  void generateSquareGridMP(const std::pair<T, T> &minCorner,
                            const std::pair<T, T> &maxCorner, T MP_size) {
    generateSquareGridMP(minCorner.first, minCorner.second, maxCorner.first,
                         maxCorner.second, MP_size);
  }

  // To be implemented
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
  T getMPSize() const { return m_MP_size; }
  const Vector<char> &getActiveNodes() const { return m_activeNodes; }
  Index nx() const { return m_nx; }
  Index ny() const { return m_ny; }
  const Vector<std::pair<T, T>> &getNodes() const { return m_nodes; }
  const Vector<std::pair<T, T>> &getMPCoords() const { return m_MPs; }

  // Cached MP element IDs
  const Vector<Index> &getMPelementIDs() const { return m_mpElementId; }
  Index getMPelementID(Index p) const {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    return m_mpElementId[p];
  }
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

  // Boundary-contact helpers for MPs (axis-aligned rectangular domain)
  // Bitmask values: left=1, right=2, bottom=4, top=8
  static constexpr char MP_CONTACT_LEFT = 1 << 0;
  static constexpr char MP_CONTACT_RIGHT = 1 << 1;
  static constexpr char MP_CONTACT_BOTTOM = 1 << 2;
  static constexpr char MP_CONTACT_TOP = 1 << 3;

  // Returns a per-MP contact bitmask based on support radius R.
  // Contact condition:
  //   left:   x - R <= 0
  //   right:  x + R >= L
  //   bottom: y - R <= 0
  //   top:    y + R >= H
  Vector<char> getMPBoundaryContactMask(T R) const {
    assert(R >= T{0} && "Support radius R must be non-negative");
    Vector<char> mask(static_cast<std::size_t>(m_nMPs));
    for (Index p{0}; p < m_nMPs; ++p) {
      const T x = m_MP_x[p];
      const T y = m_MP_y[p];
      char m = 0;
      const T xMinus = x - R;
      const T xPlus = x + R;
      const T yMinus = y - R;
      const T yPlus = y + R;

      if (xMinus < T{0} ||
          approximatelyEqualAbsRel(static_cast<double>(xMinus), 0.0))
        m |= MP_CONTACT_LEFT;
      if (xPlus > m_length ||
          approximatelyEqualAbsRel(static_cast<double>(xPlus),
                                   static_cast<double>(m_length)))
        m |= MP_CONTACT_RIGHT;
      if (yMinus < T{0} ||
          approximatelyEqualAbsRel(static_cast<double>(yMinus), 0.0))
        m |= MP_CONTACT_BOTTOM;
      if (yPlus > m_height ||
          approximatelyEqualAbsRel(static_cast<double>(yPlus),
                                   static_cast<double>(m_height)))
        m |= MP_CONTACT_TOP;
      mask[p] = m;
    }
    return mask;
  }

  // Convenience overload: use particle radius R = MP_size/2
  Vector<char> getMPBoundaryContactMask() const {
    assert(m_MP_size > T{0} && "MP_size must be > 0 to infer R = MP_size/2");
    return getMPBoundaryContactMask(m_MP_size / static_cast<T>(2));
  }

  bool isMPContactBound(Index p, T R) const {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    return getMPBoundaryContactMask(R)[p] != 0;
  }

  // Convenience overload: use particle radius R = MP_size/2
  bool isMPContactBound(Index p) const {
    assert(m_MP_size > T{0} && "MP_size must be > 0 to infer R = MP_size/2");
    return isMPContactBound(p, m_MP_size / static_cast<T>(2));
  }

  // MP setters
  void setMPCoord(Index p, T x, T y) {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    m_MP_x[p] = x;
    m_MP_y[p] = y;
    m_MPs[p] = {x, y};
    if (m_mpElementId.size() != m_nMPs) {
      m_mpElementId.resize(m_nMPs, idError);
    }
    m_mpElementId[p] = findCageID(x, y, m_mpElementId[p]);
  }

  void setMPCoords(const Vector<T> &mp_x, const Vector<T> &mp_y) {
    assert(mp_x.size() == static_cast<std::size_t>(m_nMPs) &&
           mp_y.size() == static_cast<std::size_t>(m_nMPs) &&
           "Size mismatch: mp coordinates must match mesh MP count");
    m_MP_x = mp_x;
    m_MP_y = mp_y;
    if (m_mpElementId.size() != m_nMPs) {
      m_mpElementId.resize(m_nMPs, idError);
    }
    for (Index p{0}; p < m_nMPs; ++p) {
      m_MPs[p] = {mp_x[p], mp_y[p]};
      m_mpElementId[p] = findCageID(mp_x[p], mp_y[p], m_mpElementId[p]);
    }
  }

  void setMPCoords(const Vector<std::pair<T, T>> &mp_positions) {
    assert(mp_positions.size() == static_cast<std::size_t>(m_nMPs) &&
           "Size mismatch: mp_positions must match mesh MP count");
    m_MPs = mp_positions;
    if (m_mpElementId.size() != m_nMPs) {
      m_mpElementId.resize(m_nMPs, idError);
    }
    for (Index p{0}; p < m_nMPs; ++p) {
      m_MP_x[p] = mp_positions[p].first;
      m_MP_y[p] = mp_positions[p].second;
      m_mpElementId[p] = findCageID(m_MP_x[p], m_MP_y[p], m_mpElementId[p]);
    }
  }

  void updateMPCoord(Index p, T x, T y) {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    m_MP_x[p] = x;
    m_MP_y[p] = y;
    m_MPs[p] = {x, y};
    if (m_mpElementId.size() != m_nMPs) {
      m_mpElementId.resize(m_nMPs, idError);
    }
    m_mpElementId[p] = findCageID(x, y, m_mpElementId[p]);
  }

  void updateAllMPs(const Vector<T> &new_mp_x, const Vector<T> &new_mp_y) {
    assert(new_mp_x.size() == m_nMPs && new_mp_y.size() == m_nMPs &&
           "Size mismatch");
    m_MP_x = new_mp_x;
    m_MP_y = new_mp_y;
    if (m_mpElementId.size() != m_nMPs) {
      m_mpElementId.resize(m_nMPs, idError);
    }
    for (Index p{0}; p < m_nMPs; ++p) {
      m_MPs[p] = {new_mp_x[p], new_mp_y[p]};
      m_mpElementId[p] = findCageID(new_mp_x[p], new_mp_y[p], m_mpElementId[p]);
    }
  }

  void resetMPs() {
    m_MP_x = m_MP_x_initial;
    m_MP_y = m_MP_y_initial;
    if (m_mpElementId.size() != m_nMPs) {
      m_mpElementId.resize(m_nMPs, idError);
    }
    for (Index p{0}; p < m_nMPs; ++p) {
      m_MPs[p] = {m_MP_x_initial[p], m_MP_y_initial[p]};
      m_mpElementId[p] = findCageID(m_MP_x[p], m_MP_y[p], m_mpElementId[p]);
    }
  }

  // Update cached element IDs using the mesh-internal MP coordinates.
  void updateMPElementIds() {
    if (m_mpElementId.size() != m_nMPs) {
      m_mpElementId.resize(m_nMPs, idError);
    }
    for (Index p{0}; p < m_nMPs; ++p) {
      m_mpElementId[p] = findCageID(m_MP_x[p], m_MP_y[p], m_mpElementId[p]);
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

  // Boundary node sets (useful for applying BCs)
  // y < eps, y > H-eps, x < eps, x > L-eps
  Vector<Index> bottomNodes(double absEps = 1e-8, double relEps = 1e-8) const {
    Vector<Index> ids;
    for (Index i{0}; i < m_nNodes; ++i) {
      const double y = static_cast<double>(m_nodes_y[i]);
      if (y < absEps || approximatelyEqualAbsRel(y, 0.0, absEps, relEps))
        ids.push_back(i);
    }
    return ids;
  }

  Vector<Index> topNodes(double absEps = 1e-8, double relEps = 1e-8) const {
    Vector<Index> ids;
    const double H = static_cast<double>(m_height);
    for (Index i{0}; i < m_nNodes; ++i) {
      const double y = static_cast<double>(m_nodes_y[i]);
      if (y > H - absEps || approximatelyEqualAbsRel(y, H, absEps, relEps))
        ids.push_back(i);
    }
    return ids;
  }

  Vector<Index> leftNodes(double absEps = 1e-8, double relEps = 1e-8) const {
    Vector<Index> ids;
    for (Index i{0}; i < m_nNodes; ++i) {
      const double x = static_cast<double>(m_nodes_x[i]);
      if (x < absEps || approximatelyEqualAbsRel(x, 0.0, absEps, relEps))
        ids.push_back(i);
    }
    return ids;
  }

  Vector<Index> rightNodes(double absEps = 1e-8, double relEps = 1e-8) const {
    Vector<Index> ids;
    const double L = static_cast<double>(m_length);
    for (Index i{0}; i < m_nNodes; ++i) {
      const double x = static_cast<double>(m_nodes_x[i]);
      if (x > L - absEps || approximatelyEqualAbsRel(x, L, absEps, relEps))
        ids.push_back(i);
    }
    return ids;
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

  Index findCageID(T x, T y, Index lastElement = idError) const {
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
    return idError; // Not found (outside domain)
  }

  //   Activate nodes that contain Material Points
  void activateNodes() {
    if (m_activeNodes.size() != m_nNodes) {
      m_activeNodes.resize(m_nNodes);
    }
    m_activeNodes.resetZero();

    if (m_mpElementId.size() != m_nMPs) {
      m_mpElementId.resize(m_nMPs, idError);
    }
    for (Index p{0}; p < m_nMPs; ++p) {
      Index elemID = m_mpElementId[p];
      if (elemID == idError) {
        elemID = findCageID(m_MP_x[p], m_MP_y[p], elemID);
        m_mpElementId[p] = elemID;
      }
      if (elemID != idError) {
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