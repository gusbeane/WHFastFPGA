#include <array>
#include <immintrin.h>
#include "util.h"
#include "whfast.h"
#include "whfast_kernel.h"
#include "whfast_constants.h"

// Prepare structure-of-arrays for vectorization; actual integration to be implemented
// Integrate up to time tmax using timestep dt; returns status
int whfast_integrate(std::array<Body, N_BODIES>& solarsystem, Body *com, double dt, long Nint)
{
    double x_vec[N_PLANETS], y_vec[N_PLANETS], z_vec[N_PLANETS];
    double vx_vec[N_PLANETS], vy_vec[N_PLANETS], vz_vec[N_PLANETS];
    double m_vec[N_PLANETS];

    inertial_to_democraticheliocentric_posvel(solarsystem, com, x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec);

    // Calculate necessary constants
    // Constructs a struct called kConsts with the constants
    initialize_constants(solarsystem[0].mass, m_vec);

    // Perform the initial half-drift step
    whfast_drift_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, com, dt / 2.0);

    // --------------------------------------
    // Now call the main kernel
    // --------------------------------------
    whfast_kernel(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, com, dt, Nint-1);

    // Final iteration happens outside loop to avoid branching
    // Perform the jump step (first half)
    whfast_jump_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt / 2.0);

    // Perform the interaction step
    whfast_interaction_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt);

    // Perform the jump step (second half)
    whfast_jump_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt / 2.0);

    // Perform the final half-drift step
    whfast_drift_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, com, dt / 2.0);

    // Convert back to inertial coordinates and store into solarsystem
    democraticheliocentric_to_inertial_posvel(solarsystem, com, x_vec, y_vec, z_vec,
        vx_vec, vy_vec, vz_vec, m_vec);

    return 0;
}
