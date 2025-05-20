#include <array>
#include <cmath>
#include <cstdio>
#include <immintrin.h>
#include "util.h"

// Helper to get hex string of a double
std::string double_to_hex(double d) {
    union { double d; uint64_t u; } u;
    u.d = d;
    std::ostringstream oss;
    oss << std::hex << std::setw(16) << std::setfill('0') << u.u;
    return oss.str();
}

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
                                               double *x_vec, double *y_vec, double *z_vec,
                                               double *vx_vec, double *vy_vec, double *vz_vec,
                                               double *m_vec)
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

    for(std::size_t i=0; i<N_BODIES-1; i++)
    {
        m_vec[i] = bodies[i+1].mass;
        x_vec[i] = bodies[i+1].pos[0];
        y_vec[i] = bodies[i+1].pos[1];
        z_vec[i] = bodies[i+1].pos[2];
        vx_vec[i] = bodies[i+1].vel[0];
        vy_vec[i] = bodies[i+1].vel[1];
        vz_vec[i] = bodies[i+1].vel[2];
    }
}

void democraticheliocentric_to_inertial_posvel(std::array<Body, N_BODIES> &bodies, Body *com,
                                                double *x_vec,  double *y_vec,  double *z_vec,
                                                double *vx_vec,  double *vy_vec,  double *vz_vec,
                                                double *m_vec)
{
    // Compute center of mass of planets
    struct Body com_planets = {0};
    for (std::size_t i = 0; i < N_BODIES - 1; i++)
    {
        com_planets.mass += m_vec[i];
        com_planets.pos[0] += m_vec[i] * x_vec[i];
        com_planets.pos[1] += m_vec[i] * y_vec[i];
        com_planets.pos[2] += m_vec[i] * z_vec[i];
        com_planets.vel[0] += m_vec[i] * vx_vec[i];
        com_planets.vel[1] += m_vec[i] * vy_vec[i];
        com_planets.vel[2] += m_vec[i] * vz_vec[i];

        // Add com vel to set bodies velocities
        bodies[i + 1].vel[0] = vx_vec[i] + com->vel[0];
        bodies[i + 1].vel[1] = vy_vec[i] + com->vel[1];
        bodies[i + 1].vel[2] = vz_vec[i] + com->vel[2];
    }

    for (std::size_t i = 0; i < 3; ++i)
    {
        // We divide by total mass of the system
        com_planets.pos[i] /= com->mass;
    }

    // Set solar position and velocity
    bodies[0].pos[0] = com->pos[0] - com_planets.pos[0];
    bodies[0].pos[1] = com->pos[1] - com_planets.pos[1];
    bodies[0].pos[2] = com->pos[2] - com_planets.pos[2];
    bodies[0].vel[0] = com->vel[0] - com_planets.vel[0];
    bodies[0].vel[1] = com->vel[1] - com_planets.vel[1];
    bodies[0].vel[2] = com->vel[2] - com_planets.vel[2];

    // Set planet positions
    for (std::size_t i = 0; i < N_BODIES - 1; i++)
    {
        bodies[i + 1].pos[0] = x_vec[i] + bodies[0].pos[0];
        bodies[i + 1].pos[1] = y_vec[i] + bodies[0].pos[1];
        bodies[i + 1].pos[2] = z_vec[i] + bodies[0].pos[2];
    }

    return;
}