#include <array>
#include <cmath>
#include <cstdio>
#include "util.h"

// Moves all bodies so that the system's center of mass is at the origin
void move_to_center_of_mass(std::array<Body, N_BODIES>& bodies) {
    std::array<double, 3> com_pos = {0.0, 0.0, 0.0};
    std::array<double, 3> com_vel = {0.0, 0.0, 0.0};
    double total_mass = 0.0;

    for (const auto& b : bodies) {
        total_mass += b.mass;
        for (int i = 0; i < 3; ++i) {
            com_pos[i] += b.mass * b.pos[i];
            com_vel[i] += b.mass * b.vel[i];
        }
    }
    for (int i = 0; i < 3; ++i) {
        com_pos[i] /= total_mass;
        com_vel[i] /= total_mass;
    }
    for (auto& b : bodies) {
        for (int i = 0; i < 3; ++i) {
            b.pos[i] -= com_pos[i];
            b.vel[i] -= com_vel[i];
        }
    }
}

// Transforms positions and velocities from barycentric inertial to democratic heliocentric coordinates
void inertial_to_democraticheliocentric_posvel(std::array<Body, N_BODIES>& bodies, Body *com) {
    // Assume bodies[0] is the central star

    // Shift planets to heliocentric coordinates (relative to star)
    for (std::size_t i = 1; i < N_BODIES; i++) {
        for (int j = 0; j < 3; ++j) {
            bodies[i].pos[j] -= bodies[0].pos[j];
        }
    }
    for( int j = 0; j < 3; j++) {
        bodies[0].pos[j] -= bodies[0].pos[j]; // star position is unchanged
    }

    *com = {};
    
    // Compute center of mass velocity of solar system
    for (std::size_t i = 0; i < N_BODIES; ++i) {
        com->mass += bodies[i].mass;
        for (int j = 0; j < 3; ++j) {
            com->vel[j] += bodies[i].mass * bodies[i].vel[j];
            com->pos[j] += bodies[i].mass * bodies[i].pos[j];
        }
    }
    for (int j = 0; j < 3; ++j) {
        com->vel[j] /= com->mass;
        com->pos[j] /= com->mass;
    }

    // Make solar system at rest
    for (std::size_t i = 1; i < N_BODIES; ++i) {
        for (int j = 0; j < 3; ++j) {
            bodies[i].vel[j] -= com->vel[j];
        }
    }

//     // Debug: print body 0 and body 7 positions and velocities
//     printf("particle 1 pos=(%e, %e, %e) vel=(%e, %e, %e)\n",
//             bodies[1].pos[0], bodies[1].pos[1], bodies[1].pos[2],
//             bodies[1].vel[0], bodies[1].vel[1], bodies[1].vel[2]);
    
//     printf("particle 8 pos=(%e, %e, %e) vel=(%e, %e, %e)\n",
//             bodies[8].pos[0], bodies[8].pos[1], bodies[8].pos[2],
//             bodies[8].vel[0], bodies[8].vel[1], bodies[8].vel[2]);


// x[1] = -0.050606442387163809
// y[1] = -0.46199153044180374
// z[1] = -0.033112314887493777
// vx[1] = 1.2978664296086504
// vy[1] = -0.095245414385742383
// vz[1] = -0.1267757435000321
// x[8] = 29.796803635571621
// y[8] = -2.5135110786022508
// z[8] = -0.63492520913961792
// vx[8] = 0.014142948749760034
// vy[8] = 0.18292110904074921
// vz[8] = -0.0040928456197298641



}

