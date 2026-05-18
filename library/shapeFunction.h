#ifndef SHAPE_FUNCTION_MPM_H
#define SHAPE_FUNCTION_MPM_H

#include "Mesh.h"
#include "Vector.h"
#include "interpolate.h" // cubicBSpline1..4, cubicBSpline1D
#include "signFunction.h"
#include <cmath>
#include <iostream>

// max nodes is the number of biggest possibly number of affected nodes
// e.g : hat function 1D: MaxNodes = 2, cubicBSpline 2D: MaxNodes = 16;
template <typename T, int MaxNodes> struct shapeFunction {
  StaticVector<T, MaxNodes> N{};
  StaticVector<T, MaxNodes> dN_dx{};
  StaticVector<T, MaxNodes> dN_dy{}; // only use in 2D
  StaticVector<Index, MaxNodes> nodeIdx{};
  Index count{0}; // real number of influenced nodes
};

enum class shapePolicy : unsigned char {
  linear,
  quadraticBSpline,
  cubicBSpline,
  Bezier,
  policyCount,
};

namespace bspline_detail {

// Modified B-Splines 3rd order polynomial (k=3 and C^2 continuous)
//  => 4 types of shape functions
// Supports of node I spans 4 elements round it: -2 <= r <= 2 so that S_I != 0
// r = (x_p _ x_I) / h

// Type 1: xI = xB (node I located at the boundary)
template <typename T> T cubicBSpline1(T r) {
  if (T{-2} <= r && r <= T{-1})
    return T{1} / T{6} * r * r * r + r * r + T{2} * r + T{4} / T{3};
  else if (T{-1} <= r && r <= T{0})
    return T{-1} / T{6} * r * r * r + r + T{1};
  else if (T{0} <= r && r <= T{1})
    return T{1} / T{6} * r * r * r - r + T{1};
  else if (T{1} <= r && r <= T{2})
    return T{-1} / T{6} * r * r * r + r * r - T{2} * r + T{4} / T{3};
  else
    return T{0};
}
template <typename T> T dCubicBSpline1_dr(T r) {
  if (T{-2} <= r && r < T{-1})
    return T{0.5} * r * r + T{2} * r + T{2};
  else if (T{-1} <= r && r < T{0})
    return T{-0.5} * r * r + T{1};
  else if (T{0} <= r && r < T{1})
    return T{0.5} * r * r - T{1};
  else if (T{1} <= r && r <= T{2})
    return T{-0.5} * r * r + T{2} * r - T{2};
  else
    return T{0};
}

// Type 2: xI = xB + h (node I located on the right side of the closest boundary
// + 1 cell away)
template <typename T> T cubicBSpline2(T r) {
  if (T{-1} <= r && r <= T{0})
    return (T{-1} / T{3}) * r * r * r - r * r + (T{2} / T{3});
  else if (T{0} <= r && r <= T{1})
    return (T{1} / T{2}) * r * r * r - r * r + (T{2} / T{3});
  else if (T{1} <= r && r <= T{2})
    return (T{-1} / T{6}) * r * r * r + r * r - T{2} * r + (T{4} / T{3});
  else
    return T{0};
}

template <typename T> T dCubicBSpline2_dr(T r) {
  if (T{-1} <= r && r < T{0})
    return -r * r - T{2} * r;
  else if (T{0} <= r && r < T{1})
    return T{1.5} * r * r - T{2} * r;
  else if (T{1} <= r && r <= T{2})
    return T{-0.5} * r * r + T{2} * r - T{2};
  else
    return T{0};
}

// Type 3: xI >= xB + 2h (node I located at atleast 2 cells away from any
// boundary)
template <typename T> T cubicBSpline3(T r) {
  if (T{-2} <= r && r <= T{-1})
    return T{1} / T{6} * r * r * r + r * r + T{2} * r + T{4} / T{3};
  else if (T{-1} <= r && r <= T{0})
    return (T{-1} / T{2}) * r * r * r - r * r + (T{2} / T{3});
  else if (T{0} <= r && r <= T{1})
    return (T{1} / T{2}) * r * r * r - r * r + (T{2} / T{3});
  else if (T{1} <= r && r <= T{2})
    return T{-1} / T{6} * r * r * r + r * r - T{2} * r + T{4} / T{3};
  else
    return T{0};
}
template <typename T> T dCubicBSpline3_dr(T r) {
  if (T{-2} <= r && r < T{-1})
    return T{0.5} * r * r + T{2} * r + T{2};
  else if (T{-1} <= r && r < T{0})
    return T{-1.5} * r * r - T{2} * r;
  else if (T{0} <= r && r < T{1})
    return T{1.5} * r * r - T{2} * r;
  else if (T{1} <= r && r <= T{2})
    return T{-0.5} * r * r + T{2} * r - T{2};
  else
    return T{0};
}

// Type 4: xI = xB - h (node I located on the left side of the closest
// boundary + 1 cell away)
template <typename T> T cubicBSpline4(T r) {
  if (T{-2} <= r && r <= T{-1})
    return T{1} / T{6} * r * r * r + r * r + T{2} * r + T{4} / T{3};
  else if (T{-1} <= r && r <= T{0})
    return T{-1} / T{2} * r * r * r - r * r + T{2} / T{3};
  else if (T{0} <= r && r <= T{1})
    return T{-1} / T{3} * r * r * r - r * r + T{2} / T{3};
  else
    return T{0};
}

template <typename T> T dCubicBSpline4_dr(T r) {
  if (T{-2} <= r && r < T{-1})
    return T{0.5} * r * r + T{2} * r + T{2};
  else if (T{-1} <= r && r < T{0})
    return T{-1.5} * r * r - T{2} * r;
  else if (T{0} <= r && r < T{1})
    return -r * r - T{2} * r;
  else
    return T{0};
}

// Dispatch by ShapeType
template <typename T> T evalN(ShapeType type, T r) {
  switch (type) {
  case ShapeType::boundary:
    return cubicBSpline1(r);
  case ShapeType::nearLeft:
    return cubicBSpline2(r);
  case ShapeType::nearRight:
    return cubicBSpline4(r);
  case ShapeType::middle:
    return cubicBSpline3(r);
  default:
    return T{0};
  }
}

template <typename T> T evaldN(ShapeType type, T r) {
  switch (type) {
  case ShapeType::boundary:
    return dCubicBSpline1_dr(r);
  case ShapeType::nearLeft:
    return dCubicBSpline2_dr(r);
  case ShapeType::nearRight:
    return dCubicBSpline4_dr(r);
  case ShapeType::middle:
    return dCubicBSpline3_dr(r);
  default:
    return T{0};
  }
}
}; // namespace bspline_detail

// -------------------------------------------------------
// computeSF — 1D
// Returns shapeFunction<T, 4> (max: cubic spline = 4 nodes)
// Actual used nodes stored in result.count
// -------------------------------------------------------
template <typename T>
shapeFunction<T, 4> computeSF(T xp, const Mesh1D<T> &mesh,
                              shapePolicy sfPolicy) {
  const T h = mesh.getLengthEle(0); // uniform
  const Index nNodes = mesh.getNumNodes();
  shapeFunction<T, 4> result{}; // Initialize once, reuse for all paths
  result.count = 0;

  if (sfPolicy == shapePolicy::linear) {
    // ---- Linear hat: 2 nodes affected
    assert(nNodes >= 2 && "Mesh must be discretized with more nodes!");

    Index e = mesh.findCageID(xp);
    if (e == idError) {
      std::cout << "Particle at " << xp << " move out of domain!\n";
      return result;
    }

    for (Index k = 0; k < 2; ++k) {
      T xI = mesh.getNode(e + k).pos;
      result.nodeIdx[k] = e + k;
      T dist = std::abs(xp - xI);
      if (dist / h <= T{1}) {
        result.N[k] = T{1} - dist / h;
        result.dN_dx[k] = -sgn(xp - xI) / h;

      } else {
        result.N[k] = T{0};
        result.dN_dx[k] = T{0};
      }
    }
    result.count = 2; // always = 2
    return result;
  } else if (sfPolicy == shapePolicy::cubicBSpline) {
    // ---- Cubic Spline: 4 nodes affected
    assert(nNodes >= 4 && "Mesh must be discretized with more nodes!");

    Index e = mesh.findCageID(xp);
    if (e == idError) {
      std::cout << "Particle at " << xp << " move out of domain!\n";
      return result;
    }
    Index start = e - 1; // index start of supports
    Index count = 0;

    for (Index k = start; k <= start + 3; ++k) {
      // Skip the node that gone out of domain
      if (k < 0 || k >= mesh.getNumNodes())
        continue;

      T xI = mesh.getNode(k).pos;
      T r = (xp - xI) / h; // normalized distance

      ShapeType type = mesh.getNode(k).shapeType;

      result.N[count] = bspline_detail::evalN(type, r);
      // dN_dx = dN_dr * dr_dx = dN_dr * 1/h
      result.dN_dx[count] = bspline_detail::evaldN(type, r) / h;
      result.nodeIdx[count] = k;
      ++count;
    }
    result.count = count; //<=4
    return result;
  } else {
    return result; // fallback: empty
  }
}

// -------------------------------------------------------
// computeSF — 2D (tensor product)
// N_IJ(x,y)   = N_I(x) * N_J(y)
// dN_IJ/dx    = dN_I/dx * N_J(y) + 0
// dN_IJ/dy    = N_I(x)  * dN_J/dy + 0
// Returns shapeFunction<T, 16> (max possible nodes: 4x4 for cubic)
// Actual used nodes stored in result.count
// -------------------------------------------------------
template <typename T>
shapeFunction<T, 16> computeSF(T xp, T yp, const Mesh2D<T> &mesh,
                               shapePolicy sfPolicy) {
  const Index nx = mesh.nx();
  const Index ny = mesh.ny();
  const T hx = mesh.getGridLength() / static_cast<T>(nx - 1);
  const T hy = mesh.getGridHeight() / static_cast<T>(ny - 1);
  shapeFunction<T, 16> result{}; // Initialize once
  result.count = 0;

  if (sfPolicy == shapePolicy::linear) {
    // Linear: 2x2 = 4 nodes
    assert(nx >= 2 && ny >= 2 && "Mesh must be discretized with more nodes!");

    const Index elemID = mesh.findCageID(xp, yp);
    if (elemID == idError) {
      std::cout << "Particle at " << xp << " move out of domain!\n";
      return result;
    }

    const Index ix = mesh.getElement(elemID).idx.x();
    const Index iy = mesh.getElement(elemID).idx.y();
    int count = 0;

    for (int j = iy; j <= iy + 1; ++j) {
      T yJ = mesh.getNode(j * nx).pos.y();
      T ry = (yp - yJ) / hy;
      T Ny = T{1} - std::abs(ry);
      T dNy = (ry >= T{0}) ? -T{1} / hy : T{1} / hy;

      for (int i = ix; i <= ix + 1; ++i) {
        T xI = mesh.getNode(i).pos.x();
        T rx = (xp - xI) / hx;
        T Nx = T{1} - std::abs(rx);
        T dNx = (rx >= T{0}) ? -T{1} / hx : T{1} / hx;

        result.N[count] = Nx * Ny;
        result.dN_dx[count] = dNx * Ny;
        result.dN_dy[count] = Nx * dNy;
        result.nodeIdx[count] = j * nx + i;
        ++count;
      }
    }
    result.count = count; // <= 4
    return result;

  } else if (sfPolicy == shapePolicy::cubicBSpline) {
    // ---- Cubic B-Spline:  4x4 = 16 nodes ----
    assert(nx >= 4 && ny >= 4 && "Mesh must be discretized with more nodes!");

    const Index elemID = mesh.findCageID(xp, yp);
    if (elemID == idError)
      return result;

    const Index ix = mesh.getElement(elemID).idx.x();
    const Index iy = mesh.getElement(elemID).idx.y();
    Index count = 0;

    for (Index j = iy - 1; j <= iy + 2; ++j) {
      if (j < 0 || j >= ny)
        continue;

      T yJ = mesh.getNode(j * nx).pos.y();
      T ry = (yp - yJ) / hy;
      ShapeType typeY = mesh.getNode(j * nx).shapeType.y();
      T Ny = bspline_detail::evalN(typeY, ry);
      T dNy = bspline_detail::evaldN(typeY, ry) / hy;

      for (int i = ix - 1; i <= ix + 2; ++i) {
        if (i < 0 || i >= nx)
          continue;

        T xI = mesh.getNode(i).pos.x();
        T rx = (xp - xI) / hx;
        ShapeType typeX = mesh.getNode(i).shapeType.x();
        T Nx = bspline_detail::evalN(typeX, rx);
        T dNx = bspline_detail::evaldN(typeX, rx) / hx;

        result.N[count] = Nx * Ny;
        result.dN_dx[count] = dNx * Ny;
        result.dN_dy[count] = Nx * dNy;
        result.nodeIdx[count] = j * nx + i;
        ++count;
      }
    }
    result.count = count; // <= 16
    return result;

  } else {
    return result; // fallback: empty
  }
}

#endif
