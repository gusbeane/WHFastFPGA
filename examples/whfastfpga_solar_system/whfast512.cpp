#include "whfast512_kernel.h"
#include <array>
#include "whfast512.h"
#include "util.h"

// Prepare structure-of-arrays for vectorization; actual integration to be implemented
// Integrate up to time tmax using timestep dt; returns status
int whfast512_integrate(std::array<Body, N_BODIES>& solarsystem, double dt, long Nint)
{
    // Prepare AVX-512 vectors for bodies 1-8 (ignore solarsystem[0])
    __m512d m_vec  = _mm512_setr_pd(
        solarsystem[1].mass, solarsystem[2].mass, solarsystem[3].mass, solarsystem[4].mass,
        solarsystem[5].mass, solarsystem[6].mass, solarsystem[7].mass, solarsystem[8].mass
    );
    __m512d x_vec  = _mm512_setr_pd(
        solarsystem[1].pos[0], solarsystem[2].pos[0], solarsystem[3].pos[0], solarsystem[4].pos[0],
        solarsystem[5].pos[0], solarsystem[6].pos[0], solarsystem[7].pos[0], solarsystem[8].pos[0]
    );
    __m512d y_vec  = _mm512_setr_pd(
        solarsystem[1].pos[1], solarsystem[2].pos[1], solarsystem[3].pos[1], solarsystem[4].pos[1],
        solarsystem[5].pos[1], solarsystem[6].pos[1], solarsystem[7].pos[1], solarsystem[8].pos[1]
    );
    __m512d z_vec  = _mm512_setr_pd(
        solarsystem[1].pos[2], solarsystem[2].pos[2], solarsystem[3].pos[2], solarsystem[4].pos[2],
        solarsystem[5].pos[2], solarsystem[6].pos[2], solarsystem[7].pos[2], solarsystem[8].pos[2]
    );
    __m512d vx_vec = _mm512_setr_pd(
        solarsystem[1].vel[0], solarsystem[2].vel[0], solarsystem[3].vel[0], solarsystem[4].vel[0],
        solarsystem[5].vel[0], solarsystem[6].vel[0], solarsystem[7].vel[0], solarsystem[8].vel[0]
    );
    __m512d vy_vec = _mm512_setr_pd(
        solarsystem[1].vel[1], solarsystem[2].vel[1], solarsystem[3].vel[1], solarsystem[4].vel[1],
        solarsystem[5].vel[1], solarsystem[6].vel[1], solarsystem[7].vel[1], solarsystem[8].vel[1]
    );
    __m512d vz_vec = _mm512_setr_pd(
        solarsystem[1].vel[2], solarsystem[2].vel[2], solarsystem[3].vel[2], solarsystem[4].vel[2],
        solarsystem[5].vel[2], solarsystem[6].vel[2], solarsystem[7].vel[2], solarsystem[8].vel[2]
    );
}
