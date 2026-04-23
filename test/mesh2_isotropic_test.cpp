// Test đơn giản kiểm tra isotropic support cho Mesh2.h
#include "../library/Mesh2.h"
#include <cassert>
#include <iostream>

int main() {
  // Thông số test
  double radius = 0.5;
  int nx = 4, ny = 3;
  double x0 = 0.0, y0 = 0.0;
  double lx = nx * radius, ly = ny * radius;

  // Tạo mesh với MPs vuông cạnh radius
  Mesh2 mesh(x0, y0, lx, ly, nx, ny, radius);

  // Kiểm tra số lượng MPs
  size_t nMP = mesh.m_mp.size();
  std::cout << "So luong MP: " << nMP << std::endl;
  assert(nMP == nx * ny);

  // Kiểm tra từng MP
  for (const auto &mp : mesh.m_mp) {
    // Kiểm tra radius
    assert(mp.radius == radius);
    // Kiểm tra vị trí nằm trong domain
    assert(mp.x >= x0 && mp.x <= x0 + lx);
    assert(mp.y >= y0 && mp.y <= y0 + ly);
  }

  // Kiểm tra mask (nếu có hàm updateMask)
  if constexpr (requires(Mesh2 m) { m.updateMask(); }) {
    mesh.updateMask();
    std::cout << "updateMask() goi thanh cong." << std::endl;
  }

  std::cout << "Test isotropic support PASSED!\n";
  return 0;
}
