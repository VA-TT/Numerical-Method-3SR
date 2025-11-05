#ifndef GLOBAL_CONSTANTS_HPP
#define GLOBAL_CONSTANTS_HPP

#include <limits>

namespace constants {
// Mathematical constant pi
[[maybe_unused]] inline constexpr double pi{3.14159265358979323846};
[[maybe_unused]] inline constexpr double avogadro{6.02214076e23};
[[maybe_unused]] inline constexpr double speed_of_light{299792458.0};
[[maybe_unused]] inline constexpr double planck{6.62607015e-34};
[[maybe_unused]] inline constexpr double gravity{9.80665};
// Numeric (machine) epsilon for double precision
[[maybe_unused]] inline constexpr double machine_epsilon{
    std::numeric_limits<double>::epsilon()};
[[maybe_unused]] inline constexpr double relEps{1e-8};
[[maybe_unused]] inline constexpr double absEps{1e-11};

} // namespace constants

#endif