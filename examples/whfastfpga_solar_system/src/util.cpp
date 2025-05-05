#include <array>
#include <cmath>
#include <cstdio>
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
void inertial_to_democraticheliocentric_posvel(std::array<Body, N_BODIES> &bodies, Body *com)
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
}
