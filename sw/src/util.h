#ifndef UTIL_H
#define UTIL_H

#include <array>
#include <immintrin.h>

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
                                               __m512d *x_vec, __m512d *y_vec, __m512d *z_vec,
                                               __m512d *vx_vec, __m512d *vy_vec, __m512d *vz_vec,
                                               __m512d *m_vec);

// Transforms positions and velocities from democratic heliocentric to inertial coordinates
void democraticheliocentric_to_inertial_posvel(std::array<Body, N_BODIES> &bodies, Body *com,
                                               __m512d x_vec, __m512d y_vec, __m512d z_vec,
                                               __m512d vx_vec, __m512d vy_vec, __m512d vz_vec,
                                               __m512d m_vec);

#endif // UTIL_H