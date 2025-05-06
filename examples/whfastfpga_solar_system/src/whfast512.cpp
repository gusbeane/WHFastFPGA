#include <array>
#include <immintrin.h>
#include "util.h"
#include "whfast512.h"
#include "whfast512_kernel.h"
#include "whfast512_constants.h"

// Prepare structure-of-arrays for vectorization; actual integration to be implemented
// Integrate up to time tmax using timestep dt; returns status
int whfast512_integrate(std::array<Body, N_BODIES>& solarsystem, Body *com, double dt, long Nint)
{
    // Prepare AVX-512 vectors for bodies 1-8 (ignore solarsystem[0])
    __m512d x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec;

    inertial_to_democraticheliocentric_posvel(solarsystem, com, &x_vec, &y_vec, &z_vec,
        &vx_vec, &vy_vec, &vz_vec, &m_vec);

    // Calculate necessary constants
    // Constructs a struct called kConsts with the constants
    initialize_constants(solarsystem[0].mass, m_vec);

    // Now call the main kernel
    whfast512_kernel(
        &x_vec, &y_vec, &z_vec,
        &vx_vec, &vy_vec, &vz_vec,
        m_vec, com, dt, Nint
    );

    // // Update positions and velocities of the bodies using a loop with AVX-512
    // for (int i = 1; i < N_BODIES; i++)
    // {
    //     solarsystem[i].pos[0] = _mm512_extractf64x2_pd(x_vec, (i - 1) / 4)[(i - 1) % 4];
    //     solarsystem[i].pos[1] = _mm512_extractf64x2_pd(y_vec, (i - 1) / 4)[(i - 1) % 4];
    //     solarsystem[i].pos[2] = _mm512_extractf64x2_pd(z_vec, (i - 1) / 4)[(i - 1) % 4];
    //     solarsystem[i].vel[0] = _mm512_extractf64x2_pd(vx_vec, (i - 1) / 4)[(i - 1) % 4];
    //     solarsystem[i].vel[1] = _mm512_extractf64x2_pd(vy_vec, (i - 1) / 4)[(i - 1) % 4];
    //     solarsystem[i].vel[2] = _mm512_extractf64x2_pd(vz_vec, (i - 1) / 4)[(i - 1) % 4];
    // }
    return 0;
}
