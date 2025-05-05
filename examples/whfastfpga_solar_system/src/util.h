#ifndef UTIL_H
#define UTIL_H

#include <array>

struct Body {
    std::array<double, 3> pos;
    std::array<double, 3> vel;
    double mass;
};

constexpr std::size_t N_BODIES = 9;

void move_to_center_of_mass(std::array<Body, N_BODIES>& bodies);

// Transforms positions and velocities to democratic heliocentric coordinates
void inertial_to_democraticheliocentric_posvel(std::array<Body, N_BODIES>& bodies, Body *com);

#endif // UTIL_H