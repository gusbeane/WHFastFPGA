#ifndef UTIL_H
#define UTIL_H

#include <array>
#include <immintrin.h>
#include <iomanip>
#include <iostream>
#include <sstream>

std::string double_to_hex(double d);

// Add prototype for converting a long to hex string
std::string long_to_hex(long l);

struct Body {
    std::array<double, 3> pos;
    std::array<double, 3> vel;
    double mass;
};

constexpr std::size_t N_BODIES = 9;
constexpr std::size_t N_PLANETS = N_BODIES - 1;

struct Body compute_com(std::array<Body, N_BODIES>& bodies);

void move_to_center_of_mass(std::array<Body, N_BODIES>& bodies);

// Transforms positions and velocities to democratic heliocentric coordinates
void inertial_to_democraticheliocentric_posvel(std::array<Body, N_BODIES> &bodies, Body *com,
                                               double *x_vec, double *y_vec, double *z_vec,
                                               double *vx_vec, double *vy_vec, double *vz_vec,
                                               double *m_vec);

// Transforms positions and velocities from democratic heliocentric to inertial coordinates
void democraticheliocentric_to_inertial_posvel(std::array<Body, N_BODIES> &bodies, Body *com,
                                               double *x_vec, double *y_vec, double *z_vec,
                                               double *vx_vec, double *vy_vec, double *vz_vec,
                                               double *m_vec);

#endif // UTIL_H