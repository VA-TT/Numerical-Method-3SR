#ifndef DEM_LIBRARY_H
#define DEM_LIBRARY_H

#include "Particle.h"
#include "Vector.h"
#include "cassert"
#include "ioDirectory.h"
#include "physicConstants.h"
#include "random.h"
#include "signFunction.h"
#include <iostream>

template <typename T> struct Particle2D {
  StaticVector<T, 2> pos{};
  StaticVector<T, 2> vel{};
  StaticVector<T, 2> acc{};

  // quaternion Q; //only for 3D
  StaticVector<T, 2> vrot{};
  StaticVector<T, 2> arot{};

  T radius{};
  T mass{};
  T inertia{};

  T kn{}, kt{};

  StaticVector<T, 2> force{};  // fx & fy
  StaticVector<T, 2> moment{}; // frot
}

struct Interaction {
  Index i{};
  Index j{};
  T fn{};
  T ft{};
  bool touch { false; }
  Interaction(Index i, Index j) {};
}

periodicBC(){};

reducedCoordinate(){};

template <typename T, Index ngw, Index ngh> class DEM2D {
private:
  constexpr Index m_npc{ngh * ngh};
  StaticVector<Particle2D<T>, m_npc> m_particles {}
  DynamicVector<Particle2D<T>> m_interactions {}
  T m_z{};       // Nombre de coordination
  T m_e{};       // Indice de vide
  T m_I{};       // Nombre d'inertie
  T m_kappa{};   // Niveau de raideur
  T m_dt{};      // Pas de temps
  T m_eta{0.95}; // Amortissement
  T m_phi{0.64}; // Solid Fraction
  T m_p{1000};   // Isotropic compression
  T m_H0{};      // hauteur initiale
  T m_v{};       // Velocité de cissailement
  T m_muy{};     // Coefficient de frottement
  T m_c{};       // Cohesion
  T m_dmax{};    // distance maximale pour etre voisin

public:
  DEM2D(T rho, T rmin, T rmax) {
    T totalM{};
    T totalR{};
    T totalV{};
    for (auto &p : particles) {
      p.radius = Random::get(rmin, rmax);
      p.mass = rho * constants::pi * p.radius * p.radius;
      p.inertia = T{0.5} * p.mass * p.radis * p.radius;
      totalM += p.mass;
      totalR += p.radius;
      totalV += constant::pi * p.radius * p.radius;
      //   p.pos = 0;
    }
    T meanMass = totalM / static_cast<T>(npc);
    T meanRadius = totalR / static_cast<T>(npc);
    T kn = m_particles[0].kn;
    T kappa = kn / (m_p);
    T dt_c = constants::pi * std::sqrt(meanMass / kn) / T{20};
    assert(kappa >= 1000 && "Kappa must >= 1000.");
    assert(m_dt <= dt_c && "Time steps must be smaller than the critical one.");
    double I = m_v / m_H0 * std::sqrt(meanMass / (m_p * meanRadius * T{2}));
    assert(I <= 1e-4 &&
           "Inertia number should be <1e-4 to atteint quasi-static regime.");
    m_dmax = 0.95 * rmin;
  }
  DEM2D() = default;
  DEM2D(const DEM2D &) = default;
  DEM2D(DEM2D &&) = default;
  DEM2D &operator=(const DEM2D &) = default;
  DEM2D &operator=(DEM2D &&) = default;
  ~DEM2D() = default;

  verifyContact() {}
  updateNeighborList() {}
  timeIntegration() {}
  velocityVerlet() {
    computeHalfVel();
    computeHalfVel();
    computeAcc();
    computeFullVel();
  }
  computeHalfVel() {}
  updatePosition() {}
  computeAcc() {}
  computeFullVel() {}
  coulombFriction() {}
  updateNeighbotList() {}
}

#endif