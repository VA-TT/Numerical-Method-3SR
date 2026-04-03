#ifndef MY_PARTICLE_H
#define MY_PARTICLE_H

#include "Vector.h"

template <typename T> struct Particle2D {
  // position
  StaticVector<T, 2> pos{};
  StaticVector<T, 2> vel{};
  StaticVector<T, 2> acc{};

  // Rotation
  T rot{};
  T vrot{};
  T arot{};

  // Properties
  T radius{};
  T mass{};
  T inertia{};

  // Raideur
  T kn{}, kt{};

  // force et moment de résultat (Cartesian)
  StaticVector<T, 2> force{}; // fx & fy
  T moment{};                 // torque
};

template <typename T> struct Particle3D {
  // position
  StaticVector<T, 3> pos{};
  StaticVector<T, 3> vel{};
  StaticVector<T, 3> acc{};

  // Rotation
  StaticVector<T, 3> rot{};
  StaticVector<T, 3> vrot{};
  StaticVector<T, 3> arot{};

  // Properties
  T radius{};
  T mass{};
  StaticVector<T, 3> inertia{};

  // Raideur
  T kn{}, kt{};
  T visco{};

  // force et moment de résultat (Cartesian)
  StaticVector<T, 3> force{};  // fx & fy & fz
  StaticVector<T, 3> moment{}; // frot

  // Quaternion q{};
};

template <typename T> struct Interaction {
  Index i{};
  Index j{};
  T dn{}; //  Overlaps
  bool touch{false};

  T fn{};
  T ft{};
  Interaction(Index i_, Index j_) : i(i_), j(j_) {}
};

#endif // MY_PARTICLE_H