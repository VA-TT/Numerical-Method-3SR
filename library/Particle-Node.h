#ifndef MY_PARTICLE_AND_NODE_H
#define MY_PARTICLE_AND_NODE_H

#include "Vector.h"

// Boundary-contact helpers for MPs (axis-aligned rectangular domain)
// Bitmask values: left=1, right=2, bottom=4, top=8
enum class MPContact : unsigned char {
  None = 0,
  Left = 1u << 0,   // 0001
  Right = 1u << 1,  // 0010
  Bottom = 1u << 2, // 0100
  Top = 1u << 3     // 1000
};

inline MPContact operator|(MPContact a, MPContact b) {
  return static_cast<MPContact>(static_cast<unsigned char>(a) |
                                static_cast<unsigned char>(b));
}

inline MPContact operator&(MPContact a, MPContact b) {
  return static_cast<MPContact>(static_cast<unsigned char>(a) &
                                static_cast<unsigned char>(b));
}

inline MPContact &operator|=(MPContact &a, MPContact b) {
  a = a | b;
  return a;
}

inline bool hasFlag(MPContact mask, MPContact flag) {
  return static_cast<unsigned char>(mask & flag) != 0;
}

static constexpr Index idError = -1; // default id error

template <typename T> struct Node1D {
  T m{};
  T pos{}, posInit{};
  T v{};
  T a{};
  T P{};        // Momentum
  char vCon{0}; // character 1 for yes, 0 for no
  char aCon{0};
  char pCon{0};
  char fCon{0};

  // Constraint values (so constraints can be set anytime and enforced later)
  T vConVal{};
  T aConVal{};
  T pConVal{};
  T fConVal{};
  T bodyF{}, tracF{}; // boday force / traction force
  // Nodal external/internal forces
  T extF{}, intF{}, totF{};

  bool isActive{false};
  Index eleID{idError};
};

template <typename T> struct Node2D {
  T mass{};
  StaticVector<T, 2> pos{}, posInit{};
  StaticVector<T, 2> v{};
  StaticVector<T, 2> a{};
  StaticVector<T, 2> P{};
  StaticVector<char, 2> vCon{};
  StaticVector<char, 2> aCon{};
  StaticVector<char, 2> pCon{};
  StaticVector<char, 2> FCon{};

  // Constraint values (so constraints can be set anytime and enforced later)
  StaticVector<T, 2> vConVal{};
  StaticVector<T, 2> aConVal{};
  StaticVector<T, 2> pConVal{};
  StaticVector<T, 2> FConVal{};
  StaticVector<T, 2> bodyF{}, tracF{};
  // Nodal forces
  StaticVector<T, 2> extF{}, intF{}, totF{};

  bool isActive{false};
  Index eleID{idError};
};

template <typename T> struct Particle1D {
  // position
  T pos{};
  T v{}; // vel
  T a{};

  // Properties
  T m{};
  T V{}, V0{}; // Volume

  // force
  T f{};

  // MPM 1D state
  T P{};
  T sig{};
  T eps{};
  T epsDot{}; // strain rate
  T dEps{};

  Index eleID{idError}; // default: error
};

template <typename T> struct Particle2D {
  // position
  StaticVector<T, 2> pos{};
  StaticVector<T, 2> v{};
  StaticVector<T, 2> a{};

  // Rotation
  T rot{};
  T vrot{};
  T arot{};

  // Properties
  T R{}; // radius
  T m{};
  T I{};       // moment of inertia
  T V{}, V0{}; // volume

  T kn{}, kt{}; // Raideur
  T vis{};      // viscocité

  // force et moment de résultat (Cartesian)
  StaticVector<T, 2> f{}; // fx & fy
  T M{};                  // torque/Moment

  // MPM 2D state
  StaticVector<T, 2> P{};
  Matrix<T, 2, 2> sig{};
  Matrix<T, 2, 2> eps{};
  Matrix<T, 2, 2> epsDot{}; // strain rate
  Matrix<T, 2, 2> dEps{};
  Matrix<T, 2, 2> F{}; // deformation gradient

  Index eleID{idError}; // default: error
  MPContact mask{MPContact::None};
};

// template <typename T> struct Particle3D {
//   // position
//   StaticVector<T, 3> pos{};
//   StaticVector<T, 3> v{};
//   StaticVector<T, 3> a{};

//   // Rotation
//   StaticVector<T, 3> rot{};
//   StaticVector<T, 3> vrot{};
//   StaticVector<T, 3> arot{};

//   // Properties
//   T R{};
//   T m{};
//   T V{};
//   StaticVector<T, 3> I{};

//   // Raideur
//   T kn{}, kt{};
//   T vis{};

//   // force et moment de résultat (Cartesian)
//   StaticVector<T, 3> F{}; // fx & fy & fz
//   StaticVector<T, 3> M{}; // frot

//   Index eleID{idError}; // default: error
//   MPContact mask{MPContact::None};

//   // Quaternion q{};
// };

// template <typename T> struct Interaction3D {
//   Index i{};
//   Index j{};
//   T dn{}; //  Overlaps
//   bool touch{false};

//   T fn{};
//   T ft{};
//   Interaction3D(Index i_, Index j_) : i(i_), j(j_) {}
// };

#endif // MY_PARTICLE_H