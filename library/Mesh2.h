#ifndef DISCRETIZING_RECTANGULAR_MESH_H
#define DISCRETIZING_RECTANGULAR_MESH_H

#include "ParentElement.h"
#include "Particle.h"
#include "Vector.h"
#include "comparison.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
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

  // Material Points / Particles DEM
  // (Must use Dynamic Vector as number of MPs could be 0)
  DynamicVector<Particle1D<T>> m_MPs{};

public:
  // Constructor
  Mesh1D(T length, Index nNodes, Index nMPperEle = 0)
      : m_length{length}, m_nNodes{nNodes}, m_nElements{m_nNodes - 1},
        m_nMPperEle{nMPperEle}, m_nMPs{nMPperEle * (m_nNodes - 1)} {
    assert(m_nElements > 0 && "Number of elements must be positive");
    assert(length > 0 && "Domain length must be positive");

    // reset
    m_nodes.resize(m_nNodes);
    m_elements.resize(m_nElements);

    // Position of nodes
    T lx{m_length / m_nElements};
    for (Index i{0}; i < m_nNodes; ++i) {
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
    activateNodeAndElement();
  }

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
  Index getNumMPperEle() const { return m_nMPperEle; di lam nuoc ngoai co }

  T getNodePos(Index i) const {
    assert(i >= 0 && i < m_nMPs && "Invalid node's index");
    return m_MPs[p].pos;
  }
  T getMPpos(Index p) const {
    assert(p >= 0 && p < m_nMPs && "Invalid MP's index");
    return m_MPs[p].pos;
  }

  T getLengthEle(Index e) const {
    assert(e >= 0 && e < m_nElements && "Invalid element ID");
    return std::abs(m_elements[e].x2() - m_elements[e].x1());
  }

  std::pair<T, T> getElementNodes(Index e) const {
    assert(e >= 0 && e < m_nElements && "Invalid element ID");
    return {m_elements.n1, m_elements.n2};
  }

  Index getMPelementID(Index p) const {
    assert(p >= 0 && p < m_nMPs && "Invalid MP index");
    return m_MPs[p].eleID;
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
    m_nMPperEle = number / m_nElements;
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
  }

  // Update cached element IDs using the mesh-internal MP coordinates.
  void activateNodeAndElement() {
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
        m_nodes[n1].isActive = true;
        m_nodes[n2].isActive = true;
        m_nodes[n1].eleID = activeEleID;
        m_nodes[n2].eleID = activeEleID;
      }
    }
  }

  // Update mesh for time-stepping (e.g., Updated Lagrangian FEM)
  void updateNodePosition(Index nodeID, T new_x) {
    assert(nodeID >= 0 && nodeID < m_nNodes && "Invalid node ID");
    m_nodes[nodeID] = new_x;
  }

  void updateAllNodes(const DynamicVector<T> &new_positions) {
    assert(new_positions.size() == m_nNodes && "Size mismatch");
    for (Index i{0}; i < m_nNodes; ++i) {
      m_nodes[i].pos = new_positions;
    }
    activateNodeAndElement();
  }

  // Reset to initial configuration
  void nodalReset() {
    m_nodes = m_nodes_initial;
    activateNodeAndElement();
  } // Excluding MPs

  // MPM helper function: Find element containing position x
  Index findCageID(T x) const {
    for (Index e{0}; e < m_nElements; ++e) {
      T x1 = m_nodes[m_elements.n1].pos;
      T x2 = m_nodes[m_elements.n2].pos;
      if (x1 <= x && x2 >= x)
        return e;
    }
    return idError; // Not found (outside domain)
  }

  // Print mesh info
  void print() const {
    std::cout << "=== 1D Mesh Info ===" << '\n';
    std::cout << "Number of nodes: " << m_nNodes << '\n';
    std::cout << "Number of elements: " << m_nElements << '\n';
    std::cout << "Node coordinates: " << m_nodes << '\n';
    std::cout << "Connectivity:\n";
    for (Index e = 0; e < m_nElements; ++e) {
      std::cout << "  Element " << e << ": [" << m_elements[e].n1 << ", "
                << m_elements[e].n2 << "]\n";
    }
  }
};

//////// WHAT HAPPENS IF THE MPs IS PERFECTLY LANDS ON THE NODE? ////////

#endif