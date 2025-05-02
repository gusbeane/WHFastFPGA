#include <array>
#include <cmath>
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
void inertial_to_democraticheliocentric_posvel(std::array<Body, N_BODIES>& bodies) {
    // Assume bodies[0] is the central star
    double m0 = bodies[0].mass;

    // Shift planets to heliocentric coordinates (relative to star)
    for (std::size_t i = 1; i < N_BODIES; ++i) {
        for (int j = 0; j < 3; ++j) {
            bodies[i].pos[j] -= bodies[0].pos[j];
        }
    }
    
    // Compute center of mass velocity of solar system
    std::array<double, 3> com_vel = {0.0, 0.0, 0.0};
    for (std::size_t i = 1; i < N_BODIES; ++i) {
        for (int j = 0; j < 3; ++j) {
            com_vel[j] += bodies[i].mass * bodies[i].vel[j];
        }
    }
    for (int j = 0; j < 3; ++j) {
        com_vel[j] /= m0;
    }

    // Make solar system at rest
    for (std::size_t i = 1; i < N_BODIES; ++i) {
        for (int j = 0; j < 3; ++j) {
            bodies[i].vel[j] -= com_vel[j];
        }
    }
}

