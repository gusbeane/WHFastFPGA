#include "whfast512_kernel.h"
#include <array>
#include "whfast512.h"
#include "util.h"

// Prepare structure-of-arrays for vectorization; actual integration to be implemented
// Integrate up to time tmax using timestep dt; returns status
int whfast512_integrate(std::array<Body, N_BODIES>& solarsystem, double dt, long Nint)
{
    alignas(64) double m[N_BODIES], x[N_BODIES], y[N_BODIES], z[N_BODIES];
    alignas(64) double vx[N_BODIES], vy[N_BODIES], vz[N_BODIES];
    for (std::size_t i = 0; i < N_BODIES; ++i) {
        m[i]  = solarsystem[i].mass;
        x[i]  = solarsystem[i].pos[0];
        y[i]  = solarsystem[i].pos[1];
        z[i]  = solarsystem[i].pos[2];
        vx[i] = solarsystem[i].vel[0];
        vy[i] = solarsystem[i].vel[1];
        vz[i] = solarsystem[i].vel[2];
    }
    // TODO: implement WHFast512 integration here using m, x, y, z, vx, vy, vz arrays
    return 0;
}
