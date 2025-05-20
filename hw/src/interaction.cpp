#include "interaction.h"

real_t gravity_prefactor(real_t m, real_t dx, real_t dy, real_t dz)
{
    real_t r2 = dx * dx + dy * dy + dz * dz;
    real_t r = sqrt(r2);
    return m / (r2 * r);
}

// Performs one full interaction step (scalar version)
struct bodies_t interaction_step(struct bodies_t ss, real_t M0, real_t dt)
{
    real_t c_val = 10065.32;

    real_t star_dvx = 0.0, star_dvy = 0.0, star_dvz = 0.0;
    for (int i = 0; i < N_PLANETS; ++i)
    {
        real_t gr_prefac_i  = 6.*M0*M0/(c_val*c_val);
        real_t gr_prefac2_i = ss.m_vec[i]/M0;

        real_t xi = ss.x_vec[i], yi = ss.y_vec[i], zi = ss.z_vec[i];
        real_t vxi = ss.vx_vec[i], vyi = ss.vy_vec[i], vzi = ss.vz_vec[i];
        real_t r2 = xi * xi + yi * yi + zi * zi;
        real_t r4 = r2 * r2;
        real_t invr4 = gr_prefac_i * dt / r4;
        real_t dvx_i = invr4 * xi;
        real_t dvy_i = invr4 * yi;
        real_t dvz_i = invr4 * zi;
        vxi -= dvx_i;
        vyi -= dvy_i;
        vzi -= dvz_i;
        star_dvx += gr_prefac2_i * dvx_i;
        star_dvy += gr_prefac2_i * dvy_i;
        star_dvz += gr_prefac2_i * dvz_i;
        for (int j = 0; j < N_PLANETS; ++j)
        {
            if (i == j)
                continue;
            real_t dx = xi - ss.x_vec[j];
            real_t dy = yi - ss.y_vec[j];
            real_t dz = zi - ss.z_vec[j];
            real_t pref = gravity_prefactor(ss.m_vec[j], dx, dy, dz) * dt;
            vxi -= pref * dx;
            vyi -= pref * dy;
            vzi -= pref * dz;
        }
        ss.vx_vec[i] = vxi;
        ss.vy_vec[i] = vyi;
        ss.vz_vec[i] = vzi;
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