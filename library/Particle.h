#ifndef MY_PARTICLE_H
#define MY_PARTICLE_H

#include "Vector.h"
static constexpr Index idError = -1; // default id error

template <typename T> struct Node1D {
  T mass{};
  T pos{}, posInit{};
  T vel{};
  T acc{};
  T momen{};
  bool velCon{0};
  bool accCon{0};
  bool MomenCon{0};
  bool forceCon{0};

  // Constraint values (so constraints can be set anytime and enforced later)
  T velConVal{};
  T accConVal{};
  T momenConVal{};
  T forceConVal{};
  T bodyForce{}, tractionForce{};
  // Nodal external forces
  T forceExt{}, forceInt{}, forceTot{};

  bool isActive{false};
  Index eleID{idError};
};

template <typename T> struct Node2D {
  T mass{};
  StaticVector<T, 2> pos{}, posInit{};
  StaticVector<T, 2> velocity{};
  StaticVector<T, 2> acceleration{};
  StaticVector<T, 2> nomentum{};
  bool velocityConstrained{};
  bool accelerationConstrained{};
  bool momentumConstrained{};
  bool forceConstrained{};

  // Constraint values (so constraints can be set anytime and enforced later)
  StaticVector<T, 2> velocityConstraintValue{};
  StaticVector<T, 2> accelerationConstraintValue{};
  StaticVector<T, 2> momentumConstraintValue{};
  StaticVector<T, 2> forceConstraintValue{};
  StaticVector<T, 2> bodyForce{}, tractionForce{};
  // Nodal external forces
  StaticVector<T, 2> forceExternal{}, forceInternal{}, forceTotal{};

  bool isActive{false};
  Index eleID{idError};
};

template <typename T> struct Particle1D {
  // position
  T pos{};
  T vel{};
  T acc{};

  // Properties
  T mass{};
  T volume{};

  // force
  T force{};

  // MPM 1D state
  T momentum{};
  T stress{};
  T strain{};
  T strainRate{};
  T dStrain{};

  Index eleID{idError}; // default: error
};

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
  T volume{};

  T kn{}, kt{}; // Raideur
  T vis{};      // viscocité

  // force et moment de résultat (Cartesian)
  StaticVector<T, 2> force{}; // fx & fy
  T moment{};                 // torque

  Index eleID{idError}; // default: error
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
  T volume{};
  StaticVector<T, 3> inertia{};

  // Raideur
  T kn{}, kt{};
  T visco{};

  // force et moment de résultat (Cartesian)
  StaticVector<T, 3> force{};  // fx & fy & fz
  StaticVector<T, 3> moment{}; // frot

  Index eleID{idError}; // default: error

  // Quaternion q{};
};

template <typename T> struct Interaction2D {
  Index i{};
  Index j{};
  T dn{}; //  Overlaps
  bool touch{false};

  T fn{};
  T ft{};
  Interaction2D(Index i_, Index j_) : i(i_), j(j_) {}
};

template <typename T> struct Interaction3D {
  Index i{};
  Index j{};
  T dn{}; //  Overlaps
  bool touch{false};

  T fn{};
  T ft{};
  Interaction3D(Index i_, Index j_) : i(i_), j(j_) {}
};

#endif // MY_PARTICLE_H