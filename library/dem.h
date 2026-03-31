#ifndef DEM_LIBRARY_H
#define DEM_LIBRARY_H

#include "Particle.h"
#include "Vector.h"
#include "cassert"
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

  StaticVector<T, 2> force{};
  StaticVector<T, 2> moment{};
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

template <typename T, Index npc> class DEM2D {
private:
  StaticVector<Particle2D<T>, npc> m_particles {}
  Vector<Particle2D<T>, nct> m_particles {}

public:
  DEM2D(T rho, T rmin, T rmax) {
    for (auto &p : particles) {
      p.radius = Random::get(rmin, rmax);
      p.mass = rho * constants::pi * p.radius * p.radius;
      p.inertia = T{0.5} * p.mass * p.radis * p.radius;
      p.pos = 0;
    }
  }
  DEM2D() = default;
  DEM2D(const DEM2D &) = default;
  DEM2D(DEM2D &&) = default;
  DEM2D &operator=(const DEM2D &) = default;
  DEM2D &operator=(DEM2D &&) = default;
  ~DEM2D() = default;
}

#endif