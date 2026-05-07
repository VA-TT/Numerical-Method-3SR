#ifndef SHAPE_FUNCTION_MPM_H
#define SHAPE_FUNCTION_MPM_H

#include "Node-Particle.h"
#include "interpolate.h"
#include <array>
#include <cmath>

// -------------------------------------------------------
// Kết quả shape function tại 1 particle
// MaxNodes = 2 (linear 1D), 4 (cubic 1D / linear 2D), 16 (cubic 2D)
// -------------------------------------------------------
template <typename T, int MaxNodes> struct SFResult {
  std::array<T, MaxNodes> N{};
  std::array<T, MaxNodes> dN_dx{};
  std::array<T, MaxNodes> dN_dy{}; // chỉ dùng trong 2D
  std::array<int, MaxNodes> nodeIdx{};
  int count{0};
};

struct LinearSF {};
struct CubicBSplineSF {};
struct QuadraticSF {};
struct BezierSF {};

// -------------------------------------------------------
// Derivative of 4 cubic B-Spline types
// r = (xp - xI) / h
// -------------------------------------------------------
namespace bspline_detail {

template <typename T> T dCubicBSpline1_dr(T r) {
  if (T{-2} <= r && r < T{-1})
    return T{0.5} * r * r + T{2} * r + T{2};
  if (T{-1} <= r && r < T{0})
    return T{-0.5} * r * r + T{1};
  if (T{0} <= r && r < T{1})
    return T{0.5} * r * r - T{1};
  if (T{1} <= r && r <= T{2})
    return T{-0.5} * r * r + T{2} * r - T{2};
  return T{0};
}

template <typename T> T dCubicBSpline2_dr(T r) {
  if (T{-1} <= r && r < T{0})
    return -r * r - T{2} * r;
  if (T{0} <= r && r < T{1})
    return T{1.5} * r * r - T{2} * r;
  if (T{1} <= r && r <= T{2})
    return T{-0.5} * r * r + T{2} * r - T{2};
  return T{0};
}

template <typename T> T dCubicBSpline3_dr(T r) {
  if (T{-2} <= r && r < T{-1})
    return T{0.5} * r * r + T{2} * r + T{2};
  if (T{-1} <= r && r < T{0})
    return T{-1.5} * r * r - T{2} * r;
  if (T{0} <= r && r < T{1})
    return T{1.5} * r * r - T{2} * r;
  if (T{1} <= r && r <= T{2})
    return T{-0.5} * r * r + T{2} * r - T{2};
  return T{0};
}

template <typename T> T dCubicBSpline4_dr(T r) {
  if (T{-2} <= r && r < T{-1})
    return T{0.5} * r * r + T{2} * r + T{2};
  if (T{-1} <= r && r < T{0})
    return T{-1.5} * r * r - T{2} * r;
  if (T{0} <= r && r < T{1})
    return -r * r - T{2} * r;
  return T{0};
}

// Dispatch theo ShapeType đã có sẵn trong Node của bạn
template <typename T> T evalN(T r, ShapeType type) {
  switch (type) {
  case ShapeType::boundary:
    return cubicBSpline1(r); // node 0 hoặc n-1
  case ShapeType::nearLeft:
    return cubicBSpline2(r); // node 1
  case ShapeType::nearRight:
    return cubicBSpline4(r); // node n-2
  case ShapeType::middle:
    return cubicBSpline3(r); // interior
  default:
    return T{0};
  }
}

template <typename T> T evalDN(T r, ShapeType type) {
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

} // namespace bspline_detail

// -------------------------------------------------------
// Interface chính — 1D
// Nhận Mesh1D (hoặc bất kỳ grid nào có getNode/getNumNodes/h)
// -------------------------------------------------------
template <typename SFPolicy, typename T, typename MeshType>
auto getShapeFunction1D(T xp, const MeshType &mesh) {
  const int nNodes = static_cast<int>(mesh.getNumNodes());
  const T h = mesh.getLengthEle(0); // uniform → tất cả đều bằng nhau

  if constexpr (std::is_same_v<SFPolicy, LinearSF>) {
    // Support = 2 nodes
    SFResult<T, 2> result;
    // Tìm cell chứa particle
    int i0 = static_cast<int>((xp - mesh.getNode(0).pos) / h);
    i0 = std::clamp(i0, 0, nNodes - 2);

    const int nodes[2] = {i0, i0 + 1};
    for (int k = 0; k < 2; ++k) {
      T xI = mesh.getNode(nodes[k]).pos;
      T xi = (xp - xI) / h; // [-1, 1]
      // Hat function
      result.N[k] = T{1} - std::abs(xi);
      result.dN_dx[k] = (xi >= T{0}) ? -T{1} / h : T{1} / h;
      result.nodeIdx[k] = nodes[k];
    }
    result.count = 2;
    return result;

  } else if constexpr (std::is_same_v<SFPolicy, CubicBSplineSF>) {
    // Support = 4 nodes (i0-1 .. i0+2)
    SFResult<T, 4> result;
    int i0 = static_cast<int>((xp - mesh.getNode(0).pos) / h);
    i0 = std::clamp(i0, 0, nNodes - 2);
    int start = i0 - 1; // cubic: cần node trước cell

    int count = 0;
    for (int k = start; k <= start + 3; ++k) {
      if (k < 0 || k >= nNodes)
        continue;

      T xI = mesh.getNode(k).pos;
      T r = (xp - xI) / h; // normalized, không cần parent coord

      // Lấy ShapeType từ Node — đã được Mesh set sẵn trong constructor!
      ShapeType type = mesh.getNode(k).shapeType;

      result.N[count] = bspline_detail::evalN(r, type);
      result.dN_dx[count] = bspline_detail::evalDN(r, type) / h;
      result.nodeIdx[count] = k;
      ++count;
    }
    result.count = count;
    return result;
  }
}

// -------------------------------------------------------
// Interface chính — 2D (tensor product)
// -------------------------------------------------------
template <typename SFPolicy, typename T, typename MeshType>
auto getShapeFunction2D(T xp, T yp, const MeshType &mesh) {
  const int nx = static_cast<int>(mesh.nx());
  const int ny = static_cast<int>(mesh.ny());
  const T hx = mesh.getGridLength() / static_cast<T>(nx - 1);
  const T hy = mesh.getGridHeight() / static_cast<T>(ny - 1);

  if constexpr (std::is_same_v<SFPolicy, CubicBSplineSF>) {
    SFResult<T, 16> result; // 4x4 max
    int ix = std::clamp(static_cast<int>((xp) / hx), 0, nx - 2);
    int iy = std::clamp(static_cast<int>((yp) / hy), 0, ny - 2);

    int count = 0;
    for (int j = iy - 1; j <= iy + 2; ++j) {
      if (j < 0 || j >= ny)
        continue;

      // ShapeType theo y — lấy từ node (0, j)
      ShapeType typeY = mesh.getNode(j * nx).shapeType.y();
      T yJ = mesh.getNode(j * nx).pos.y();
      T ry = (yp - yJ) / hy;
      T Ny = bspline_detail::evalN(ry, typeY);
      T dNy = bspline_detail::evalDN(ry, typeY) / hy;

      for (int i = ix - 1; i <= ix + 2; ++i) {
        if (i < 0 || i >= nx)
          continue;

        ShapeType typeX = mesh.getNode(i).shapeType.x();
        T xI = mesh.getNode(i).pos.x();
        T rx = (xp - xI) / hx;
        T Nx = bspline_detail::evalN(rx, typeX);
        T dNx = bspline_detail::evalDN(rx, typeX) / hx;

        result.N[count] = Nx * Ny;
        result.dN_dx[count] = dNx * Ny;
        result.dN_dy[count] = Nx * dNy;
        result.nodeIdx[count] = j * nx + i;
        ++count;
      }
    }
    result.count = count;
    return result;
  }
}