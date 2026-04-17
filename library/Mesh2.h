#ifndef DISCRETIZING_RECTANGULAR_MESH_H
#define DISCRETIZING_RECTANGULAR_MESH_H

#include "Particle.h"
#include "Vector.h"
#include "comparison.h"
#include "parentElement.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>

// 1D Mesh for line elements
template <typename T, Index nNodes, Index nMPperEle = 0> class Mesh1D {
private:
  T m_length{};

  constexpr Index m_nElements{nNodes - 1}, m_nMPs{nMPperEle * (nNodes - 1)};
  static_assert(m_nElements > 0 && "Number of elements must be positive");

  // Nodes
  StaticVector<Node<1D>, nNodes> m_nodes{};

  // Elements via Connectivity: [node_i, node_j]
  StaticVector<StaticVector<Index, 2>, m_nElements> m_connectivity{};

  // Material Points / Particles DEM
  // (Must use Dynamic Vector as number of MPs could be 0)
  DynamicVector<Particle1D<T>> m_MPs{};
  mutable DynamicVector<T> m_mpCoordsCache{}; // Backward-compat coordinate view
  static constexpr Index idError = -1;        // default id error

public:
  // Constructor
  Mesh1D(T length) : m_length{length} {
    assert(length > 0 && "Domain length must be positive");

    T lx{m_length / m_nElements};
    for (Index i{0}; i < nNodes; ++i) {
      m_nodes[i].pos = lx * i;             // Coordinates of nodes
      m_nodes[i].posInit = m_nodes[i].pos; // Saving intial configuration
    }

    for (Index e{0}; e < m_nElements; ++e) {
      m_connectivity[e] = {e, e + 1}; // Node index of each element
    }

    // Initialize Material Points (MPs) if needed
    if (m_nMPperEle > 0) {
      m_MPs.resize(m_nMPs);
      Index mpId{0};
      for (Index e{0}; e < m_nElements; ++e) {
        T x_start = m_nodes[e].pos;
        T le = getLengthEle(e);
        for (Index p{0}; p < m_nMPperEle; ++p) {
          m_MPs[mpId].pos = x_start + (p + 1) * le / (m_nMPperEle + 1);
          ++mpId;
        }
      }
    }
    // Initialize cached element IDs for MPs
    updateMPsEleID();
    activateNodes();
  };
  // Other defaults
  Mesh1D() = default;
  Mesh1D(const Mesh1D &) = default;
  Mesh1D(Mesh1D &&) = default;
  Mesh1D &operator=(const Mesh1D &) = default;
  Mesh1D &operator=(Mesh1D &&) = default;
  ~Mesh1D() = default;

  // Getter
  Index getNumNodes() const { return nNodes; }
  Index getNumElements() const { return m_nElements; }
  Index getNumMPs() const { return m_nMPs; }
  const DynamicVector<char> &getActiveNodes() const { return m_activeNodes; }
  const DynamicVector<T> &nodeCoords() const { return m_nodes; }
  const DynamicVector<Particle1D<T>> &getMPs() const { return m_MPs; }

  const DynamicVector<T> &getMPCoords() const {
    m_mpCoordsCache.resize(m_nMPs);
    for (Index p{0}; p < m_nMPs; ++p) {
      m_mpCoordsCache[p] = m_MPs[p].pos;
    }
    return m_mpCoordsCache;
  }

  const DynamicVector<T> &getInitialNodes() const { return m_nodes_initial; }

  Particle1D<T> getMP(Index p) const {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    return m_MPs[p];
  }

  T getMPpos(Index p) const {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    return m_MPs[p].pos;
  }
  const DynamicVector<DynamicVector<Index>> &getConnectivity() const {
    return m_connectivity;
  }

  // Get connectivity of a single element
  const DynamicVector<Index> &getEleConnectivity(Index e) const {
    assert(e >= 0 && e < m_nElements && "Invalid element ID");
    return m_connectivity[e];
  }

  T getLengthEle(Index e) const {
    assert(e >= 0 && e < m_nElements && "Invalid element ID");
    Index n1 = m_connectivity[e][0];
    Index n2 = m_connectivity[e][1];
    return std::abs(m_nodes[n2].pos - m_nodes[n1].pos);
  }

  std::pair<T, T> getElementNodes(Index e) const {
    assert(e >= 0 && e < m_nElements && "Invalid element ID");
    Index n1 = m_connectivity[e][0];
    Index n2 = m_connectivity[e][1];
    return {m_nodes[n1], m_nodes[n2]};
  }

  Index getMPelementID(Index p) const {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    return m_MPs[p].eleID;
  }

  // Setter
  void regenerateMesh() { *this = Mesh1D(m_length, nNodes, m_nMPperEle); }
  void setNumElements(Index number) {
    m_nElements = number;
    nNodes = number + 1;
    regenerateMesh();
  }
  void setNumNodes(Index number) {
    nNodes = number;
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

  void setMPCoords(const DynamicVector<T> &mp_positions) {
    assert(mp_positions.size() == m_nMPs &&
           "Size mismatch: mp_positions must match mesh MP count");
    for (Index p{0}; p < m_nMPs; ++p) {
      m_MPs[p].pos = mp_positions[p];
    }
    for (Index p{0}; p < m_nMPs; ++p) {
      m_MPs[p].eleID = findCageID(m_MPs[p].pos);
    }
  }

  void setMPs(const DynamicVector<Particle1D<T>> &particles) {
    assert(particles.size() == m_nMPs &&
           "Size mismatch: particles must match mesh MP count");
    m_MPs = particles;
    updateMPsEleID();
  }

  void setMPCoord(Index p, T x) {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    m_MPs[p].pos = x;
    // Keep cached element id consistent
    m_MPs[p].eleID = findCageID(x);
  }

  // Update cached element IDs using the mesh-internal MP coordinates.
  void updateMPsEleID() {
    for (Index p{0}; p < m_nMPs; ++p) {
      m_MPs[p].eleID = findCageID(m_MPs[p].pos);
    }
  }

  // Update mesh for time-stepping (e.g., Updated Lagrangian FEM)
  void updateNodePosition(Index nodeID, T new_x) {
    assert(nodeID >= 0 && nodeID < nNodes && "Invalid node ID");
    m_nodes[nodeID] = new_x;
  }

  void updateAllNodes(const DynamicVector<T> &new_positions) {
    assert(new_positions.size() == nNodes && "Size mismatch");
    m_nodes = new_positions;
  }

  // Reset to initial configuration
  void nodalReset() { m_nodes = m_nodes_initial; } // Excluding MPs

  // MPM helper function: Find element containing position x
  Index findCageID(T x) const {
    for (Index e{0}; e < m_nElements; ++e) {
      Index n1 = m_connectivity[e][0];
      Index n2 = m_connectivity[e][1];
      T x1 = m_nodes[n1].pos;
      T x2 = m_nodes[n2].pos;
      if (x1 <= x && x2 >= x)
        return e;
    }
    return idError; // Not found (outside domain)
  }

  // List of activated nodes that contain Material Points
  void activateNodes() {
    // Reset all nodes to inactive
    for (Index n{0}; n < nNodes; ++n) {
      m_nodes[n].isActive = false;
    }
    // Update elements containing MPs
    updateMPsEleID();

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
    assert(nodeID >= 0 && nodeID < nNodes && "Invalid node ID");
    return m_activeNodes[nodeID] != 0;
  }

  void activate

  // Print mesh info
  void print() const {
    std::cout << "=== 1D Mesh Info ===" << '\n';
    std::cout << "Number of nodes: " << nNodes << '\n';
    std::cout << "Number of elements: " << m_nElements << '\n';
    std::cout << "Node coordinates: " << m_nodes << '\n';
    std::cout << "Connectivity:\n";
    for (Index e = 0; e < m_nElements; ++e) {
      std::cout << "  Element " << e << ": [" << m_connectivity[e][0] << ", "
                << m_connectivity[e][1] << "]\n";
    }
  }
};

//////// WHAT HAPPENS IF THE MPs IS PERFECTLY LANDS ON THE NODE? ////////

#endif