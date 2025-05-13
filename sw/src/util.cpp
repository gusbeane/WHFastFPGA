#include <array>
#include <cmath>
#include <cstdio>
#include <immintrin.h>
#include "util.h"

struct Body compute_com(std::array<Body, N_BODIES> &bodies)
{
    // This is defined in a weird way just to replicate the
    // implementation in rebound
    struct Body com = {0};
    for (const auto &b : bodies)
    {
        for (int i = 0; i < 3; ++i)
        {
            com.pos[i] = com.pos[i] * com.mass + b.mass * b.pos[i];
            com.vel[i] = com.vel[i] * com.mass + b.mass * b.vel[i];
        }
        com.mass += b.mass;
        for (int i = 0; i < 3; ++i)
        {
            com.pos[i] /= com.mass;
            com.vel[i] /= com.mass;
        }
    }
    return com;
}

// Moves all bodies so that the system's center of mass is at the origin
void move_to_center_of_mass(std::array<Body, N_BODIES> &bodies)
{
    struct Body com = compute_com(bodies);
    for (auto &b : bodies)
    {
        for (int i = 0; i < 3; ++i)
        {
            b.pos[i] -= com.pos[i];
            b.vel[i] -= com.vel[i];
        }
    }
}

// Transforms positions and velocities from barycentric inertial to democratic heliocentric coordinates
void inertial_to_democraticheliocentric_posvel(std::array<Body, N_BODIES> &bodies, Body *com,
                                               __m512d *x_vec, __m512d *y_vec, __m512d *z_vec,
                                               __m512d *vx_vec, __m512d *vy_vec, __m512d *vz_vec,
                                               __m512d *m_vec)
{
    // Assume bodies[0] is the central star

    *com = {};

    // Compute center of mass velocity of solar system
    for (std::size_t i = 0; i < N_BODIES; ++i)
    {
        com->mass += bodies[i].mass;
        for (int j = 0; j < 3; ++j)
        {
            com->vel[j] += bodies[i].mass * bodies[i].vel[j];
            com->pos[j] += bodies[i].mass * bodies[i].pos[j];
        }
    }
    for (int j = 0; j < 3; ++j)
    {
        com->vel[j] /= com->mass;
        com->pos[j] /= com->mass;
    }

    // Shift planets to heliocentric coordinates (relative to star)
    for (std::size_t i = 1; i < N_BODIES; i++)
    {
        for (int j = 0; j < 3; ++j)
        {
            bodies[i].pos[j] -= bodies[0].pos[j];
        }
    }
    for (int j = 0; j < 3; j++)
    {
        bodies[0].pos[j] -= bodies[0].pos[j]; // star position is zero
    }

    // Make solar system at rest
    for (std::size_t i = 0; i < N_BODIES; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            bodies[i].vel[j] -= com->vel[j];
        }
    }

    *m_vec = _mm512_setr_pd(
        bodies[1].mass, bodies[2].mass, bodies[3].mass, bodies[4].mass,
        bodies[5].mass, bodies[6].mass, bodies[7].mass, bodies[8].mass);
    *x_vec = _mm512_setr_pd(
        bodies[1].pos[0], bodies[2].pos[0], bodies[3].pos[0], bodies[4].pos[0],
        bodies[5].pos[0], bodies[6].pos[0], bodies[7].pos[0], bodies[8].pos[0]);
    *y_vec = _mm512_setr_pd(
        bodies[1].pos[1], bodies[2].pos[1], bodies[3].pos[1], bodies[4].pos[1],
        bodies[5].pos[1], bodies[6].pos[1], bodies[7].pos[1], bodies[8].pos[1]);
    *z_vec = _mm512_setr_pd(
        bodies[1].pos[2], bodies[2].pos[2], bodies[3].pos[2], bodies[4].pos[2],
        bodies[5].pos[2], bodies[6].pos[2], bodies[7].pos[2], bodies[8].pos[2]);
    *vx_vec = _mm512_setr_pd(
        bodies[1].vel[0], bodies[2].vel[0], bodies[3].vel[0], bodies[4].vel[0],
        bodies[5].vel[0], bodies[6].vel[0], bodies[7].vel[0], bodies[8].vel[0]);
    *vy_vec = _mm512_setr_pd(
        bodies[1].vel[1], bodies[2].vel[1], bodies[3].vel[1], bodies[4].vel[1],
        bodies[5].vel[1], bodies[6].vel[1], bodies[7].vel[1], bodies[8].vel[1]);
    *vz_vec = _mm512_setr_pd(
        bodies[1].vel[2], bodies[2].vel[2], bodies[3].vel[2], bodies[4].vel[2],
        bodies[5].vel[2], bodies[6].vel[2], bodies[7].vel[2], bodies[8].vel[2]);
}

void democraticheliocentric_to_inertial_posvel(std::array<Body, N_BODIES> &bodies, Body *com,
                                               __m512d x_vec, __m512d y_vec, __m512d z_vec,
                                               __m512d vx_vec, __m512d vy_vec, __m512d vz_vec,
                                               __m512d m_vec)
{
    // First store into normal arrays
    double x[N_BODIES - 1], y[N_BODIES - 1], z[N_BODIES - 1];
    double vx[N_BODIES - 1], vy[N_BODIES - 1], vz[N_BODIES - 1];
    double m[N_BODIES - 1];
    _mm512_storeu_pd(x, x_vec);
    _mm512_storeu_pd(y, y_vec);
    _mm512_storeu_pd(z, z_vec);
    _mm512_storeu_pd(vx, vx_vec);
    _mm512_storeu_pd(vy, vy_vec);
    _mm512_storeu_pd(vz, vz_vec);
    _mm512_storeu_pd(m, m_vec);

    // Compute center of mass of planets
    struct Body com_planets = {0};
    for (int i = 0; i < N_BODIES - 1; i++)
    {
        com_planets.mass += m[i];
        com_planets.pos[0] += m[i] * x[i];
        com_planets.pos[1] += m[i] * y[i];
        com_planets.pos[2] += m[i] * z[i];
        com_planets.vel[0] += m[i] * vx[i];
        com_planets.vel[1] += m[i] * vy[i];
        com_planets.vel[2] += m[i] * vz[i];

        // Add com vel to set bodies velocities
        bodies[i + 1].vel[0] = vx[i] + com->vel[0];
        bodies[i + 1].vel[1] = vy[i] + com->vel[1];
        bodies[i + 1].vel[2] = vz[i] + com->vel[2];
    }

    for (int i = 0; i < 3; ++i)
    {
        // We divide by total mass of the system
        com_planets.pos[i] /= com->mass;
        // com_planets.vel[i] /= com->mass;
    }

    // Set solar position and velocity
    bodies[0].pos[0] = com->pos[0] - com_planets.pos[0];
    bodies[0].pos[1] = com->pos[1] - com_planets.pos[1];
    bodies[0].pos[2] = com->pos[2] - com_planets.pos[2];
    bodies[0].vel[0] = com->vel[0] - com_planets.vel[0];
    bodies[0].vel[1] = com->vel[1] - com_planets.vel[1];
    bodies[0].vel[2] = com->vel[2] - com_planets.vel[2];

    // Set planet positions
    for (int i = 0; i < N_BODIES - 1; i++)
    {
        bodies[i + 1].pos[0] = x[i] + bodies[0].pos[0];
        bodies[i + 1].pos[1] = y[i] + bodies[0].pos[1];
        bodies[i + 1].pos[2] = z[i] + bodies[0].pos[2];
    }

    return;
}