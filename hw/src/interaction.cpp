#include "interaction.h"

real_t gravity_prefactor(real_t m, real_t dx, real_t dy, real_t dz)
{
#pragma HLS inline
    real_t r2 = dx * dx + dy * dy + dz * dz;
    real_t r = sqrt(r2);
    return m / (r2 * r);
}

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
#pragma HLS pipeline II=1
#pragma HLS allocation function instances=sqrt limit=4

    real_t c_val = R(10065.32);
    real_t gr_prefac  = R(6.)*M0*M0/(c_val*c_val);

    real_t star_dvx = R(0.0), star_dvy = R(0.0), star_dvz = R(0.0);

    // Note that we will apply, e.g., vx -= dvdx + star_dvx
    real_t dvx[N_PLANETS], dvy[N_PLANETS], dvz[N_PLANETS];
#pragma HLS array_partition variable=dvx complete
#pragma HLS array_partition variable=dvy complete
#pragma HLS array_partition variable=dvz complete

    for (int i = 0; i < N_PLANETS; ++i)
    {
        real_t gr_prefac2_i = ss.m_vec[i] / M0;
        real_t xi = ss.x_vec[i], yi = ss.y_vec[i], zi = ss.z_vec[i];
        real_t vxi = ss.vx_vec[i], vyi = ss.vy_vec[i], vzi = ss.vz_vec[i];
        real_t r2 = xi * xi + yi * yi + zi * zi;
        real_t r4 = r2 * r2;
        real_t invr4 = gr_prefac * dt / r4;
        real_t dvx_i = invr4 * xi;
        real_t dvy_i = invr4 * yi;
        real_t dvz_i = invr4 * zi;
        ss.vx_vec[i] -= dvx_i;
        ss.vy_vec[i] -= dvy_i;
        ss.vz_vec[i] -= dvz_i;
        star_dvx += gr_prefac2_i * dvx_i;
        star_dvy += gr_prefac2_i * dvy_i;
        star_dvz += gr_prefac2_i * dvz_i;
    }

    for (int i = 0; i < N_PLANETS; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
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
        ss.vx_vec[i] -= star_dvx;
        ss.vy_vec[i] -= star_dvy;
        ss.vz_vec[i] -= star_dvz;
    }

    return ss;
}