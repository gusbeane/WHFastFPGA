#include "interaction.h"

real_t get_invr3(real_t dx, real_t dy, real_t dz)
{
#pragma HLS inline
    real_t r2 = dx * dx + dy * dy + dz * dz;
    real_t r = sqrt(r2);
    return R(1.) / (r2 * r);
}

// Performs one full interaction step (scalar version)
struct bodies_t interaction_step(struct bodies_t ss, real_t M0, real_t dt)
{
#pragma HLS inline off
#pragma HLS pipeline II=10

    // claude suggested the following to reduce computation in the critical loop
       // Pre-compute position differences to reduce critical path
//     real_t dx_cache[N_PLANETS][N_PLANETS];
//     real_t dy_cache[N_PLANETS][N_PLANETS];
//     real_t dz_cache[N_PLANETS][N_PLANETS];
// #pragma HLS array_partition variable=dx_cache complete dim=0
// #pragma HLS array_partition variable=dy_cache complete dim=0
// #pragma HLS array_partition variable=dz_cache complete dim=0

//     // Pre-compute all position differences
// precompute_diffs:
//     for (int i = 0; i < N_PLANETS; ++i) {
// #pragma HLS UNROLL
//         for (int j = 0; j < i; ++j) {
// #pragma HLS UNROLL
//             dx_cache[i][j] = ss.x_vec[i] - ss.x_vec[j];
//             dy_cache[i][j] = ss.y_vec[i] - ss.y_vec[j];
//             dz_cache[i][j] = ss.z_vec[i] - ss.z_vec[j];
// #pragma HLS BIND_OP variable=dx_cache[i][j] op=dsub latency=6
// #pragma HLS BIND_OP variable=dy_cache[i][j] op=dsub latency=6
// #pragma HLS BIND_OP variable=dz_cache[i][j] op=dsub latency=6
//         }
//     }

    real_t star_dvx = R(0.0), star_dvy = R(0.0), star_dvz = R(0.0);

    for (int i = 0; i < N_PLANETS; ++i)
    {
#pragma HLS UNROLL
        real_t xi = ss.x_vec[i], yi = ss.y_vec[i], zi = ss.z_vec[i];
        real_t vxi = ss.vx_vec[i], vyi = ss.vy_vec[i], vzi = ss.vz_vec[i];
        real_t r2 = xi * xi + yi * yi + zi * zi;
        real_t r4 = r2 * r2;
        real_t invr4 = ss.gr_prefac * dt / r4;
        real_t dvx_i = invr4 * xi;
        real_t dvy_i = invr4 * yi;
        real_t dvz_i = invr4 * zi;
        ss.vx_vec[i] -= dvx_i;
        ss.vy_vec[i] -= dvy_i;
        ss.vz_vec[i] -= dvz_i;
        star_dvx += ss.gr_prefac2[i] * dvx_i;
        star_dvy += ss.gr_prefac2[i] * dvy_i;
        star_dvz += ss.gr_prefac2[i] * dvz_i;
    }

    for (int i = 0; i < N_PLANETS; ++i)
    {
#pragma HLS UNROLL
        for (int j = 0; j < i; ++j)
        {
#pragma HLS UNROLL
            real_t dx = ss.x_vec[i] - ss.x_vec[j];
            real_t dy = ss.y_vec[i] - ss.y_vec[j];
            real_t dz = ss.z_vec[i] - ss.z_vec[j];
            real_t invr3 = get_invr3(dx, dy, dz) * dt;
            real_t pref_i = ss.m_vec[j] * invr3;
            ss.vx_vec[i] -= pref_i * dx;
            ss.vy_vec[i] -= pref_i * dy;
            ss.vz_vec[i] -= pref_i * dz;

            real_t pref_j = ss.m_vec[i] * invr3;
            ss.vx_vec[j] += pref_j * dx;
            ss.vy_vec[j] += pref_j * dy;
            ss.vz_vec[j] += pref_j * dz;
        }
    }

    // Apply back-reaction onto star to all planets
    for (int i = 0; i < N_PLANETS; ++i)
    {
#pragma HLS UNROLL
        ss.vx_vec[i] -= star_dvx;
        ss.vy_vec[i] -= star_dvy;
        ss.vz_vec[i] -= star_dvz;
    }

    return ss;
}