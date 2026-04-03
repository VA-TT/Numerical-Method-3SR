#ifndef DISCRETE_ELEMENT_METHOD_H
#ifndef DISCRETE_ELEMENT_METHOD_H

#include "Mesh.h"
#include "Particle.h"
#include "Vector.h"
#include "cassert"
#include "ioDirectory.h"
#include "physicConstants.h"
#include "random.h"
#include "signFunction.h"
#include <iostream>

template <typename T> struct Particle2D {
  // position
  StaticVector<T, 2> pos{};
  StaticVector<T, 2> vel{};
  StaticVector<T, 2> acc{};

  // Rotation
  StaticVector<T, 2> rot{};
  StaticVector<T, 2> vrot{};
  StaticVector<T, 2> arot{};

  // Properties
  T radius{};
  T mass{};
  T inertia{};

  // Raideur
  T kn{}, kt{};

  // force et moment de résultat (Cartesian)
  StaticVector<T, 2> force{};  // fx & fy
  StaticVector<T, 2> moment{}; // frot
}

struct Interaction {
  Index i{};
  Index j{};
  T dn{}; //  Overlaps
  bool touch { false; }

  T fn{};
  T ft{};
  Interaction(Index i, Index j) {};
}

periodicBC(){};

reducedCoordinate(){};

template <typename T, Index ngl, Index ngh> class DEM2D {
private:
  constexpr Index m_npc{ngh * ngh};
  T m_legnth{}, m_height {}
  StaticVector<Particle2D<T>, m_npc> m_particles {}
  DynamicVector<Particle2D<T>> m_interactions {}
  DynamicVector<Particle2D<T>> m_neighbors {}
  T m_z{};             // Nombre de coordination
  T m_e{};             // Indice de vide
  T m_I{};             // Nombre d'inertie
  T m_kappa{};         // Niveau de raideur
  T m_dt{};            // Pas de temps
  T m_eta{};           // Amortissement (viscocité)
  T m_alpha {}         // Amortissement artificiel
  T m_solidFraction{}; // Solid Fraction
  // T m_p{1000};          // Isotropic compression
  // T m_v{};              // Velocité de cissailement
  T m_h{};               // matrix de condition limite
  T m_muy{};             // Coefficient de frottement
  T m_c{};               // Cohesion
  T m_dmax{};            // distance maximale pour etre voisin
  T m_Vsolid{};          // Somme des veloumes des perticules
  T m_Msolid{};          // Somme des masses des perticules
  Index m_stepUpdate{50} // nombre de pas de MAJ la liste de voisinage
  T m_g{};               // aceeleration gravitaire
  T m_fnMax{};           // max de fn
  T m_dt{};
  Index m_iStep{}; // compte de pas
  Index m_nStep{}; // nombre de pas

  // L'énergie de système
  T m_kinEnergy{}, m_potentEnergy{}, m_totalEnergy0{}, m_dissiEnergy{};

public:
  DEM2D(T rho, T rmin, T rmax, T kn, T kt, T muy, T c, T eta, T dt,
        T duration = 10000.0)
      : m_muy{muy}, m_c{c}, m_dt{dt}, m_dmax{0.95 * rmin},
        m_nStep{static_cast<Index>{m_duration / m_nStep}} {

    T totalM{};
    T totalR{};
    // Initiate particles' position
    T legnth = rmax * T{2} * ngl;
    T height = rmax * T{2} * ngh;
    m_mesh = Mesh2D<T>(length, height, ngl + 1, ngw + 1);
    Index indexEle{0};
    T vRand{rmin / T{(50) * m_dt}}; // Emperical value
    for (auto &p : m_particles) {
      p.radius = Random::get(rmin, rmax);
      p.mass = rho * constants::pi * p.radius * p.radius;
      p.inertia = T{0.5} * p.mass * p.radis * p.radius;
      p.pos.x() = m_mesh.getEleCenter(Index index);
      p.pos.y() = m_mesh.getEleCenter(Index index);
      p.vel.x() = Random::get(-vRand, vRand);
      p.vel.y() = Random::get(-vRand, vRand);
      p.kn = kn;
      p.kt = kt;
      totalR += p.radius;
      m_Msolid += p.mass;
      m_Vsolid += constant::pi * p.radius * p.radius;
      ++indexEle;
    }
    T meanMass = totalM / static_cast<T>(npc);
    T meanRadius = totalR / static_cast<T>(npc);
    T meanDiameter = meanRadius * T{2};

    T dt_c = constants::pi * std::sqrt(meanMass / kn) / T{20};
    assert(m_dt <= dt_c && "Time steps must be smaller than the critical one.");
    // D'autre condition de côntrole (cisaillement)
    //  T kappa = kn / (m_p);
    //  assert(kappa >= 1000 && "Kappa must >= 1000.");
    //  double I = m_v / m_H0 * std::sqrt(meanMass / (m_p * meanRadius * T{2}));
    //  assert(I <= 1e-4 && "Inertia number should be <1e-4 to atteint
    //  quasi-static regime.");

    for (Index i{0}; i < m_npc; ++i) {
      m_totalEnergy0 += m_G * m_particles[i].mass * m_particles[i].pos.y() +
                        T{0.5} * m_particles[i].mass *
                            dotProduct(m_particles[i].vel, m_particles[i].vel);
    }
  }

  DEM2D() = default;
  DEM2D(const DEM2D &) = default;
  DEM2D(DEM2D &&) = default;
  DEM2D &operator=(const DEM2D &) = default;
  DEM2D &operator=(DEM2D &&) = default;
  ~DEM2D() = default;

  // Getters

  // Setters
  void setG() {
    for (auto &p : m_particles) {
      p.acc.y() += -constants::gravity; // Not require to compute the gravity in
                                        // the force balance equation
    }
  }
  void setAlpha(T alpha) {
    assert(T{0} < alpha && alpha < T{1} && "Alpha must be in range [0,1]");
    m_alpha = alpha;
  }
  // Main algorithms
  void computeEnergy() {
    m_potentEnergy = T{};
    m_kinEnergy = T{};
    for (Index i{0}; i < m_npc; ++i) {
      m_potentEnergy += m_G * m_particles[i].mass * m_particles[i].pos.y();
      m_kinEnergy += T{0.5} * m_particles[i].mass *
                     dotProduct(m_particles[i].vel, m_particles[i].vel);
    }
    m_dissiEnergy = m_totalEnergy0 - m_potentEnergy - m_kinEnergy;
  }

  void verifyContact() {
    for (Index i{0}; i < m_npc - 1; ++i) {
      for (Index j{i + 1}; j < m_npc; ++j) {
        StaticVector<T, 2> r_ij{m_particles[j].pos - m_particles[i].pos};
        if (magnitude(r_ij) - m_particles[i].radius - m_particles[j].pos <
            T{0}) {
          m_interactions.push_back(Interaction(i, j));
        }
      }
    }
  }
  void updateNeighborList() {}
  void timeIntegration() {}
  void fowardEuler() {
    for (auto &p : m_particles) {
      p.pos += m_dt * p.vel + T{0.5} * p.acc * m_dt;
      p.vel += m_dt * p.acc;
      p.vel *= m_alpha; // Amortissement artificiel

      p.rot += m_dt * p.vrot + T{0.5} * p.arot * m_dt;
      p.vrot += m_dt * p.arot;
      p.vrot *= m_alpha; // Amortissement artificiel
    }
  }
  void artificialDissipation() {
    for (Index i{0}; i < m_npc; ++i) {
      m_particles[i].vel *= (T{1} - alpha);
      m_particles[i].vrot *= (T{1} - alpha);
    }
  }
  void velocityVerlet() {
    computeHalfVel();
    computeHalfVel();
    computeAcc();
    computeFullVel();
  }
  void computeHalfVel() {}
  void updatePosition() {}
  void computeAcc() {}
  void computeFullVel() {}
  void coulombFriction() {}
  void updateNeighbotList() {}
}
#endif