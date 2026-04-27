#ifndef DISCRETIZING_RECTANGULAR_MESH_H
#define DISCRETIZING_RECTANGULAR_MESH_H

#include "ParentElement.h"
#include "Particle-Node.h"
#include "Vector.h"
#include "comparison.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <span>
#include <utility>

// 1D Mesh for line elements
template <typename T> class Mesh1D {
private:
  // length of domain
  T m_length{};

  // Number of nodes, elements, and material points
  Index m_nNodes{}, m_nElements{}, m_nMPs{}, m_nMPperEle{};

  // Nodes
  DynamicVector<Node1D<T>> m_nodes{};

  // Elements via Connectivity: [node_i, node_j]
  DynamicVector<ElementL2<T>> m_elements{};

  // Material Points
  DynamicVector<Particle1D<T>> m_MPs{};

public:
  // Constructor
  Mesh1D(T length, Index nNodes, Index nMPperEle = 0)
      : m_length{length}, m_nNodes{nNodes}, m_nElements{nNodes - 1},
        m_nMPperEle{nMPperEle}, m_nMPs{nMPperEle * (nNodes - 1)} {
    // Checking input
    assert(m_nElements > 0 && "Number of elements must be positive");
    assert(length > 0 && "Domain length must be positive");

    // Initiate vectors to 0
    m_nodes.resize(m_nNodes);
    m_elements.resize(m_nElements);

    // Position of nodes
    T lx{m_length / m_nElements};
    for (Index i{0}; i < nNodes; ++i) {
      m_nodes[i].pos = lx * i;             // Coordinates of nodes
      m_nodes[i].posInit = m_nodes[i].pos; // Saving intial configuration
    }

    // Node index of each element
    for (Index e{0}; e < m_nElements; ++e) {
      m_elements[e] = ElementL2<T>(e, e + 1, m_nodes);
    }

    // Generate Material Points (MPs) if needed
    if (m_nMPperEle > 0) {
      m_MPs.resize(m_nMPs);
      Index mpId{0};
      for (Index e{0}; e < m_nElements; ++e) {
        T x_start = m_elements[e].x1();
        T le = getLengthEle(e);
        for (Index p{0}; p < nMPperEle; ++p) {
          m_MPs[mpId].pos = x_start + (p + 1) * le / (nMPperEle + 1);
          ++mpId;
        }
      }
    }
    // Initialize cached element IDs for MPs
    activateNodesAndElements();
  }

  // Other defaults
  Mesh1D() = default;
  Mesh1D(const Mesh1D &) = default;
  Mesh1D(Mesh1D &&) = default;
  Mesh1D &operator=(const Mesh1D &) = default;
  Mesh1D &operator=(Mesh1D &&) = default;
  ~Mesh1D() = default;

  // Getters
  T getGridLength() const { return m_length; }
  Index getNumNodes() const { return m_nNodes; }
  Index getNumElements() const { return m_nElements; }
  Index getNumMPs() const { return m_nMPs; }
  Index getNumMPperEle() const { return m_nMPperEle; }

  T getLengthEle(Index e) const {
    assert(e >= 0 && e < m_nElements && "Invalid element ID");
    return std::abs(m_elements[e].x2() - m_elements[e].x1());
  }

  const ElementL2<T> &getElement(Index e) const {
    assert(e >= 0 && e < m_nElements && "Invalid element ID");
    return m_elements[e];
  } // No need to open element non-const ref

  const DynamicVector<ElementL2<T>> &getAllElements() const {
    return m_elements;
  } // No need to open element non-const ref

  const DynamicVector<Node1D<T>> &getAllNodes() const { return m_nodes; }
  DynamicVector<Node1D<T>> &getAllNodes() { return m_nodes; }

  const Node1D<T> &getNode(Index i) const {
    assert(i >= 0 && i < m_nNodes && "Invalid nodal index");
    return m_nodes[i];
  }
  Node1D<T> &getNode(Index i) {
    assert(i >= 0 && i < m_nNodes && "Invalid nodal index");
    return m_nodes[i];
  }
  DynamicVector<Particle1D<T>> &getAllMPs() { return m_MPs; }
  const DynamicVector<Particle1D<T>> &getAllMPs() const { return m_MPs; }
  const Particle1D<T> &getMP(Index p) const {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    return m_MPs[p];
  }
  Particle1D<T> &getMP(Index p) {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    return m_MPs[p];
  }

  // Get connectivity of a single element
  StaticVector<Index, 2> getEleConnectivity(Index e) const {
    assert(e >= 0 && e < m_nElements && "Invalid element ID");
    return {m_elements[e].n1, m_elements[e].n2};
  }

  StaticVector<T, 2> getEleCoords(Index e) const {
    assert(e >= 0 && e < m_nElements && "Invalid element ID");
    return {m_elements[e].x1(), m_elements[e].x2()};
  }

  DynamicVector<T> getActiveNodes() const {
    DynamicVector<T> activeNodes{};
    for (Index i{0}; i < m_nNodes; ++i) {
      if (m_nodes[i].isActive)
        activeNodes.push_back(i);
    }
    return activeNodes;
  }

  // Setters
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
    m_nMPperEle = number / m_nElements;
    regenerateMesh();
  }
  void setLength(T length) {
    m_length = length;
    regenerateMesh();
  }

  // MPM helper function: Find element containing position x
  Index findCageID(T x) const {
    for (Index e{0}; e < m_nElements; ++e) {
      const T x1 = m_elements[e].x1();
      const T x2 = m_elements[e].x2();
      if (x1 <= x && x2 >= x)
        return e;
    }
    return idError; // Not found (outside domain)
  }

  // Update cached element IDs using the mesh-internal MP coordinates.
  void activateNodesAndElements() {
    // Reset all nodes and elements to inactive
    for (Index n{0}; n < m_nNodes; ++n) {
      m_nodes[n].isActive = false;
      m_nodes[n].eleID = idError;
    }
    for (Index e{0}; e < m_nElements; ++e) {
      m_elements[e].isActive = false;
    }
    for (Index p{0}; p < m_nMPs; ++p) {
      Index activeEleID = findCageID(m_MPs[p].pos);
      if (activeEleID != idError) {
        m_MPs[p].eleID = activeEleID;
        m_elements[activeEleID].isActive = true;
        Index n1 = m_elements[activeEleID].n1;
        Index n2 = m_elements[activeEleID].n2;
        m_nodes[n1].isActive = m_nodes[n2].isActive = true;
        m_nodes[n1].eleID = m_nodes[n2].eleID = activeEleID;
      }
    }
  }

  // Update mesh for time-stepping (e.g., Updated Lagrangian FEM)
  void setNodeCoord(Index nodeID, T new_x) {
    assert(nodeID >= 0 && nodeID < m_nNodes && "Invalid node ID");
    m_nodes[nodeID].pos = new_x;
    activateNodesAndElements();
  }

  void setAllNodesCoords(const DynamicVector<T> &new_positions) {
    assert(new_positions.size() == m_nNodes && "Size mismatch");
    for (Index i{0}; i < m_nNodes; ++i) {
      m_nodes[i].pos = new_positions[i];
    }
    activateNodesAndElements();
  }

  // Reset to initial configuration
  void nodeReset() {
    for (Index i{0}; i < m_nNodes; ++i) {
      m_nodes[i].pos = m_nodes[i].posInit;
    }
    activateNodesAndElements();
  } // Excluding MPs

  // MP setters
  void setMPCoord(Index p, T x) {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    m_MPs[p].pos = x;
    activateNodesAndElements();
  }

  void setAllMPs(const DynamicVector<Particle1D<T>> &particles) {
    assert(particles.size() == m_nMPs &&
           "Size mismatch: particles must match mesh MP count");
    m_MPs = particles;
    activateNodesAndElements();
  }

  void setAllMpCoords(const DynamicVector<T> &mp_pos) {
    assert(mp_pos.size() == m_nMPs &&
           "Size mismatch: mp coordinates must match mesh MP count");
    for (Index p{0}; p < m_nMPs; ++p) {
      m_MPs[p].pos = mp_pos[p];
    }
    activateNodesAndElements();
  }

  // Print mesh info
  void print() const {
    std::cout << "=== 1D Mesh Info ===" << '\n';
    std::cout << "Number of nodes: " << m_nNodes << '\n';
    std::cout << "Number of elements: " << m_nElements << '\n';
    std::cout << "Node coordinates: [";
    for (Index i = 0; i < m_nNodes; ++i) {
      std::cout << m_nodes[i].pos;
      if (i + 1 < m_nNodes)
        std::cout << ", ";
    }
    std::cout << "]\n";
    std::cout << "Connectivity:\n";
    for (Index e = 0; e < m_nElements; ++e) {
      std::cout << "  Element " << e << ": [" << m_elements[e].n1 << ", "
                << m_elements[e].n2 << "]\n";
    }
  }
};

template <typename T> class Mesh2D {
private:
  // domain
  T m_length{}, m_height{};

  // Number of nodes on x & y-axis
  Index m_nx{}, m_ny{}, m_nNodes{};

  // Numbers of elements
  Index m_nElements{};

  // Number of Material Points
  Index m_nMPs{}, m_nMPperEle{};

  // Nodal object
  DynamicVector<Node2D<T>> m_nodes{};

  // Elements connectivity
  DynamicVector<ElementQ4<T>> m_elements{};

  // Material Points / Particles DEM
  DynamicVector<Particle2D<T>> m_MPs{};

  // MP grid generation methods
  void generateSquareGridMP(T x0, T y0, T x1, T y1, T MP_size) {
    // Defensive: avoid division by 0 / invalid ranges
    if (!(MP_size > T{0}) || !(x1 > x0) || !(y1 > y0)) {
      return;
    }

    Index npx = static_cast<Index>((x1 - x0) / MP_size);
    Index npy = static_cast<Index>((y1 - y0) / MP_size);
    m_nMPs = npx * npy;

    // Initialize Material Points (MPs) if needed
    if (m_nMPs > 0) {
      m_MPs.resize(m_nMPs);
      Index mpId{0};
      for (Index j{0}; j < npy; ++j) {
        for (Index i{0}; i < npx; ++i) {
          T x = x0 + (static_cast<T>(i) + static_cast<T>(0.5)) * MP_size;
          T y = y0 + (static_cast<T>(j) + static_cast<T>(0.5)) * MP_size;
          m_MPs[mpId].pos.x() = x;
          m_MPs[mpId].pos.y() = y;
          m_MPs[mpId].radius = MP_size;
          ++mpId;
        }
      }
    }
  }

  void generateMPGridPerElement(Index nMPperEleSide) {
    if (nMPperEleSide <= 0) {
      return;
    }

    m_nMPperEle = nMPperEleSide * nMPperEleSide;
    m_nMPs = m_nElements * m_nMPperEle;

    // Force dx = dy = radius for all MPs
    const T radius = m_length / static_cast<T>((m_nx - 1) * nMPperEleSide);
    const T radius2 = m_height / static_cast<T>((m_ny - 1) * nMPperEleSide);
    assert(approximatelyEqualAbsRel(radius, radius2) &&
           "Only implemented for square/disk material point!");
    m_MPs.resize(m_nMPs);

    // Element ordering follows natural ordering: i + j*(nx-1)
    Index mpID{0};
    for (Index ey{0}; ey < m_ny - 1; ++ey) {
      for (Index ex{0}; ex < m_nx - 1; ++ex) {
        const Index elemID = ex + ey * (m_nx - 1);
        const T x_left = static_cast<T>(ex) * (radius * nMPperEleSide);
        const T y_bot = static_cast<T>(ey) * (radius * nMPperEleSide);
        for (Index j{0}; j < nMPperEleSide; ++j) {
          for (Index i{0}; i < nMPperEleSide; ++i) {
            const T x =
                x_left + (static_cast<T>(i) + static_cast<T>(0.5)) * radius;
            const T y =
                y_bot + (static_cast<T>(j) + static_cast<T>(0.5)) * radius;
            m_MPs[mpID].pos.x() = x;
            m_MPs[mpID].pos.y() = y;
            m_MPs[mpID].radius = radius;
            m_MPs[mpID].eleID = elemID;
            ++mpID;
          }
        }
      }
    }
  }

  // Convenience overload: pass MP region as corners
  void generateSquareGridMP(const StaticVector<T, 2> &minCorner,
                            const StaticVector<T, 2> &maxCorner, T MP_size) {
    generateSquareGridMP(minCorner.x(), minCorner.y(), maxCorner.x(),
                         maxCorner.y(), MP_size);
  }

  // Generate MPs on a grid and keep only points inside a circle.
  void generateCircleMPgrid(const StaticVector<T, 2> &center, T radius,
                            T MP_size) {
    if (!(MP_size > T{0}) || !(radius > T{0})) {
      return;
    }

    const auto floorToSpacing = [&](T v) -> T {
      const double dv = static_cast<double>(v);
      const double ds = static_cast<double>(MP_size);
      return static_cast<T>(std::floor(dv / ds) * ds);
    };

    const T gridOriginX = floorToSpacing(center.x()); // floor to grid node
    const T gridOriginY = floorToSpacing(center.y());

    const T xmin = gridOriginX - radius;
    const T ymin = gridOriginY - radius;
    const T xmax = gridOriginX + radius;
    const T ymax = gridOriginY + radius;

    // Regenerate MPs: clear previous state first
    m_MPs = DynamicVector<Particle2D<T>>{};

    const double dx = static_cast<double>(xmax - xmin);
    const double dy = static_cast<double>(ymax - ymin);
    const double ds = static_cast<double>(MP_size);
    if (!(dx >= 0.0) || !(dy >= 0.0)) {
      return;
    }

    const Index nxPts = static_cast<Index>(std::floor(dx / ds)) + 1;
    const Index nyPts = static_cast<Index>(std::floor(dy / ds)) + 1;

    const double approxCount = (constants::pi * static_cast<double>(radius) *
                                static_cast<double>(radius)) /
                               (ds * ds);
    const Index reserveCount =
        static_cast<Index>(std::max(0.0, approxCount * 1.1));
    m_MPs.reserve(reserveCount);

    const double cx = static_cast<double>(center.x());
    const double cy = static_cast<double>(center.y());
    const double r2 = static_cast<double>(radius) * static_cast<double>(radius);

    for (Index j{0}; j < nyPts; ++j) {
      const double y = static_cast<double>(ymin) + static_cast<double>(j) * ds;
      for (Index i{0}; i < nxPts; ++i) {
        const double x =
            static_cast<double>(xmin) + static_cast<double>(i) * ds;
        const double dxp = x - cx;
        const double dyp = y - cy;
        if (dxp * dxp + dyp * dyp <= r2) {
          Particle2D<T> mp{};
          mp.pos.x() = static_cast<T>(x);
          mp.pos.y() = static_cast<T>(y);
          mp.radius = MP_size;
          m_MPs.push_back(mp);
        }
      }
    }

    m_nMPs = static_cast<Index>(m_MPs.size());
  }

public:
  // ===== Constructors =====
  // 1) Explicit 4-parameter constructor (mesh only)
  Mesh2D(T length, T height, Index nx, Index ny)
      : m_length{length}, m_height{height}, m_nx{nx}, m_ny{ny},
        m_nNodes{nx * ny}, m_nElements{(nx - 1) * (ny - 1)} {
    assert(nx > 1 && ny > 1 && "Need at least 2 nodes in each direction");
    assert(length > 0 && height > 0 && "Domain dimensions must be positive");

    // Resize vectors before using them
    m_nodes.resize(m_nNodes);
    m_elements.resize(m_nElements);

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
        m_nodes[nodeID].pos.x() = dx * i;
        m_nodes[nodeID].pos.y() = dy * j;
        m_nodes[nodeID].posInit = m_nodes[nodeID].pos;
      }
    }

    // Generate connectivity for Q4 elements
    // Element node ordering (counterclockwise from bottom-left):
    //   n4 ------ n3
    //   |          |
    //   |          |
    //   n1 ------ n2

    Index elemID{0};
    for (Index j{0}; j < ny - 1; j++) {
      for (Index i{0}; i < nx - 1; i++) {
        Index n1 = i + j * nx;
        Index n2 = (i + 1) + j * nx;
        Index n3 = (i + 1) + (j + 1) * nx;
        Index n4 = i + (j + 1) * nx;
        m_elements[elemID] = ElementQ4<T>(n1, n2, n3, n4, m_nodes);
        ++elemID;
      }
    }
  }

  // 2) Explicit full-parameter constructor (mesh + MP region)
  Mesh2D(T length, T height, Index nx, Index ny, T x0, T y0, T x1, T y1,
         T MP_size)
      : Mesh2D(length, height, nx, ny) {
    if ((MP_size > T{0}) && (x1 > x0) && (y1 > y0)) {
      generateSquareGridMP(x0, y0, x1, y1, MP_size);
    }

    activateNodesAndElements();
    updateAllMasks();
  }

  // Pair-based convenience ctor: delegates to the full-parameter constructor.
  Mesh2D(const StaticVector<T, 2> &gridDimension,
         const StaticVector<Index, 2> &gridnPoints,
         const StaticVector<T, 2> &minCorner,
         const StaticVector<T, 2> &maxCorner, T MP_size)
      : Mesh2D(gridDimension.x(), gridDimension.y(), gridnPoints.x(),
               gridnPoints.y(), minCorner.x(), minCorner.y(), maxCorner.x(),
               maxCorner.y(), MP_size) {}

  // 3) MPs-per-element constructor
  Mesh2D(T length, T height, Index nx, Index ny, Index nMPperEle)
      : Mesh2D(length, height, nx, ny) {
    assert(nMPperEle > 0 && "nMPperEle must be > 0");
    generateMPGridPerElement(nMPperEle);

    activateNodesAndElements();
    updateAllMasks();
  }

  // 4) Explicit full-parameter circle constructor
  Mesh2D(T length, T height, Index nx, Index ny, const std::pair<T, T> &center,
         T radius, T MP_size)
      : Mesh2D(length, height, nx, ny) {
    if ((MP_size > T{0}) && (radius > T{0})) {
      generateCircleMPgrid(StaticVector<T, 2>{center.first, center.second},
                           radius, MP_size);

      activateNodesAndElements();
      updateAllMasks();
    }
  }

  // Rule-of-zero style defaults
  Mesh2D() = default;
  Mesh2D(const Mesh2D &) = default;
  Mesh2D(Mesh2D &&) = default;
  Mesh2D &operator=(const Mesh2D &) = default;
  Mesh2D &operator=(Mesh2D &&) = default;
  ~Mesh2D() = default;

  // ===== Getters =====
  // Geometry and state
  T getGridLength() const { return m_length; }
  T getGridHeight() const { return m_height; }
  Index getNumNodes() const { return m_nNodes; }
  Index getNumElements() const { return m_nElements; }
  Index getNumMPs() const { return m_nMPs; }

  Index nx() const { return m_nx; }
  Index ny() const { return m_ny; }

  const ElementQ4<T> &getElement(Index elemID) const {
    assert(elemID >= 0 && elemID < m_nElements && "Invalid element ID");
    return m_elements[elemID];
  }
  const DynamicVector<ElementQ4<T>> &getAllElements() const {
    return m_elements;
  }

  const DynamicVector<Node2D<T>> &getAllNodes() const { return m_nodes; }
  DynamicVector<Node2D<T>> &getAllNodes() { return m_nodes; }

  const Node2D<T> &getNode(Index nodeID) const {
    assert(nodeID >= 0 && nodeID < m_nNodes && "Invalid node ID");
    return m_nodes[nodeID];
  }
  Node2D<T> &getNode(Index nodeID) {
    assert(nodeID >= 0 && nodeID < m_nNodes && "Invalid node ID");
    return m_nodes[nodeID];
  }
  DynamicVector<Particle2D<T>> &getAllMPs() { return m_MPs; }
  const DynamicVector<Particle2D<T>> &getAllMPs() const { return m_MPs; }
  // Outer layers
  const Particle2D<T> &getMP(Index p) const {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    return m_MPs[p];
  }
  Particle2D<T> &getMP(Index p) {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    return m_MPs[p];
  }

  StaticVector<Index, 4> getEleConnectivity(Index e) const {
    assert(e >= 0 && e < m_nElements && "Invalid element ID");
    return {m_elements[e].n1, m_elements[e].n2, m_elements[e].n3,
            m_elements[e].n4};
  }

  // Backward-compatible element-node accessor
  StaticVector<StaticVector<T, 4>, 2> getEleCoords(Index elemID) const {
    return {m_elements[elemID].getX_nodes(), m_elements[elemID].getY_nodes()};
  }

  // Compatibility helpers
  DynamicVector<T> getActiveNodes() const {
    DynamicVector<T> activeNodes{};
    for (Index i{0}; i < m_nNodes; ++i) {
      if (m_nodes[i].isActive)
        activeNodes.push_back(i);
    }
    return activeNodes;
  }

  DynamicVector<Index> getActiveElementsID() const {
    DynamicVector<Index> activeIDs{};
    for (Index e{0}; e < m_nElements; ++e) {
      if (m_elements[e].isActive)
        activeIDs.push_back(e);
    }
    return activeIDs;
  }

  DynamicVector<Index> getAllMPsInElement(Index elemID) const {
    assert(elemID >= 0 && elemID < m_nElements && "Invalid element ID");
    DynamicVector<Index> ids;
    for (Index p{0}; p < m_nMPs; ++p) {
      if (m_MPs[p].eleID == elemID) {
        ids.push_back(p);
      }
    }
    return ids;
  }

  // Boundary node sets filtered by active nodes
  // Useful when applying BCs that must only target nodes participating in
  // the current time step (i.e., nodes connected to MPs).
  // Contact-mask helpers
  void updateMask(Index p) {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    // radius is interpreted as support half-size (half side length).
    const T radius = m_MPs[p].radius;
    if (!(radius > T{0})) {
      m_MPs[p].mask = MPContact::None; // reset mask
      return;
    }

    const T x = m_MPs[p].pos.x();
    const T y = m_MPs[p].pos.y();

    MPContact mask = MPContact::None;
    if (x - radius <= T{0} ||
        approximatelyEqualAbsRel(static_cast<double>(x - radius), 0.0))
      mask |= MPContact::Left;
    if (x + radius >= m_length ||
        approximatelyEqualAbsRel(static_cast<double>(x + radius),
                                 static_cast<double>(m_length)))
      mask |= MPContact::Right;
    if (y - radius <= T{0} ||
        approximatelyEqualAbsRel(static_cast<double>(y - radius), 0.0))
      mask |= MPContact::Bottom;
    if (y + radius >= m_height ||
        approximatelyEqualAbsRel(static_cast<double>(y + radius),
                                 static_cast<double>(m_height)))
      mask |= MPContact::Top;
    m_MPs[p].mask = mask;
  }

  void updateAllMasks() {
    for (Index p{0}; p < m_nMPs; ++p) {
      updateMask(p);
    }
  }

  // Active-boundary node predicates
  bool isBottomActiveNode(Index i) const {
    if (!m_nodes[i].isActive)
      return false;
    const auto y = m_nodes[i].pos.y();
    return approximatelyEqualAbsRel(y, T{0});
  }

  bool isTopActiveNode(Index i) const {
    if (!m_nodes[i].isActive)
      return false;
    const auto y = m_nodes[i].pos.y();
    return approximatelyEqualAbsRel(y, m_height);
  }

  bool isLeftActiveNode(Index i) const {
    if (!m_nodes[i].isActive)
      return false;
    const auto x = m_nodes[i].pos.x();
    return approximatelyEqualAbsRel(x, T{0});
  }

  bool isRightActiveNode(Index i) const {
    if (!m_nodes[i].isActive)
      return false;
    const auto x = m_nodes[i].pos.x();
    return approximatelyEqualAbsRel(x, m_length);
  }

  DynamicVector<Index> bottomActiveNodes() const {
    DynamicVector<Index> ids;
    for (Index i{0}; i < m_nNodes; ++i) {
      if (isBottomActiveNode(i))
        ids.push_back(i);
    }
    return ids;
  }

  DynamicVector<Index> topActiveNodes() const {
    DynamicVector<Index> ids;
    for (Index i{0}; i < m_nNodes; ++i) {
      if (isTopActiveNode(i))
        ids.push_back(i);
    }
    return ids;
  }

  DynamicVector<Index> leftActiveNodes() const {
    DynamicVector<Index> ids;
    for (Index i{0}; i < m_nNodes; ++i) {
      if (isLeftActiveNode(i))
        ids.push_back(i);
    }
    return ids;
  }

  DynamicVector<Index> rightActiveNodes() const {
    DynamicVector<Index> ids;
    for (Index i{0}; i < m_nNodes; ++i) {
      if (isRightActiveNode(i))
        ids.push_back(i);
    }
    return ids;
  }

  // ===== Setters =====
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

  // Nodal updates for time-stepping (e.g., Updated Lagrangian FEM)
  void setNodeCoord(Index nodeID, StaticVector<T, 2> new_pos) {
    assert(nodeID >= 0 && nodeID < m_nNodes && "Invalid node ID");
    m_nodes[nodeID].pos.x() = new_pos.x();
    m_nodes[nodeID].pos.y() = new_pos.y();
    // activateNodesAndElements();
    // too expensive whenever update only 1 node
    // so, remember to update in the main algorithm instead of here
  }

  void setAllNodesCoords(const DynamicVector<Node2D<T>> &new_nodes) {
    assert(new_nodes.size() == m_nNodes && "Size mismatch");
    m_nodes = new_nodes;
    activateNodesAndElements();
  }

  // Reset to initial nodal configuration
  void nodeReset() {
    for (Index i{0}; i < m_nNodes; ++i) {
      m_nodes[i].pos = m_nodes[i].posInit;
    }
    activateNodesAndElements();
  } // Excluding MPs

  // MP setters
  void setMPCoord(Index p, T x, T y) {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    m_MPs[p].pos.x() = x;
    m_MPs[p].pos.y() = y;
    updateMask(p);
    // activateNodesAndElements();
  }

  void setAllMPs(const DynamicVector<Particle2D<T>> &particles) {
    assert(particles.size() == m_nMPs &&
           "Size mismatch: particles must match mesh MP count");
    m_MPs = particles;
    updateAllMasks();
    activateNodesAndElements();
  }

  void setAllMpCoords(const DynamicVector<T> &mp_x,
                      const DynamicVector<T> &mp_y) {
    assert(mp_x.size() == m_nMPs && mp_y.size() == m_nMPs &&
           "Size mismatch: mp coordinates must match mesh MP count");
    for (Index p{0}; p < m_nMPs; ++p) {
      m_MPs[p].pos.x() = mp_x[p];
      m_MPs[p].pos.y() = mp_y[p];
    }
    updateAllMasks();
    activateNodesAndElements();
  }

  // ===== Queries and activation =====
  // MPM helper function
  bool isPointInElement(T x, T y, Index e) const {
    assert(e >= 0 && e < m_nElements && "Invalid element ID");
    auto [xi, eta] = m_elements[e].parentCoor(x, y);
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

  // Update active nodes/elements and MP->element mapping
  void activateNodesAndElements() {
    // Reset all nodes and elements to inactive
    for (Index n{0}; n < m_nNodes; ++n) {
      m_nodes[n].isActive = false;
      m_nodes[n].eleID = idError;
    }
    for (Index e{0}; e < m_nElements; ++e) {
      m_elements[e].isActive = false;
    }
    for (Index p{0}; p < m_nMPs; ++p) {
      Index activeEleID = findCageID(m_MPs[p].pos.x(), m_MPs[p].pos.y());
      if (activeEleID != idError) {
        m_MPs[p].eleID = activeEleID;
        m_elements[activeEleID].isActive = true;
        auto [n1, n2, n3, n4] = getConnectivity(activeEleID);
        m_nodes[n1].isActive = m_nodes[n2].isActive = m_nodes[n3].isActive =
            m_nodes[n4].isActive = true;
        m_nodes[n1].eleID = m_nodes[n2].eleID = m_nodes[n3].eleID =
            m_nodes[n4].eleID = activeEleID;
      } else {
        m_MPs[p].eleID = idError;
      }
    }
  }

  // ===== Debug =====
  // Print mesh info
  void print() const {
    std::cout << "=== 2D Mesh Info ===" << '\n';
    std::cout << "Number of nodes: " << m_nNodes << " (" << m_nx << " x "
              << m_ny << ")\n";
    std::cout << "Number of elements: " << m_nElements << '\n';
  }
};
#endif
