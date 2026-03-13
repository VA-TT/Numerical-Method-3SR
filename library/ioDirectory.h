#pragma once

#include <filesystem>
#include <iostream>
#include <system_error>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace ioFile {

inline std::filesystem::path executableDir() {
#ifdef _WIN32
  {
    std::wstring buffer(MAX_PATH, L'\0');
    for (;;) {
      const DWORD len = GetModuleFileNameW(nullptr, buffer.data(),
                                           static_cast<DWORD>(buffer.size()));
      if (len == 0) {
        break;
      }
      if (len < buffer.size() - 1) {
        buffer.resize(len);
        return std::filesystem::path(buffer).parent_path();
      }
      buffer.resize(buffer.size() * 2);
    }
  }
#endif
#ifdef __linux__
  {
    std::error_code ec;
    const auto exe = std::filesystem::read_symlink("/proc/self/exe", ec);
    if (!ec) {
      return exe.parent_path();
    }
  }
#endif
  return std::filesystem::current_path();
}

// Returns <executableDir>/<subdir> and creates the directory.
inline std::filesystem::path vtkOutputDir(const std::string &subdir,
                                          bool log = true) {
  const auto dir = executableDir() / subdir;
  std::error_code ec;
  std::filesystem::create_directories(dir, ec);
  if (log) {
    if (ec) {
      std::cerr << "WARNING: cannot create VTK output dir: " << dir.string()
                << " (" << ec.message() << ")\n";
    } else {
      std::cout << "VTK output dir: " << dir.string() << "\n";
    }
  }
  return dir;
}

} // namespace ioFile
