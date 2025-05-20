#include <stdio.h>
#include "whfast_kernel.h"
#include "whfast_constants.h"
#include <cmath>
#include <limits>

// Scalar Stiefel function for Halley's method, returning Gs0, Gs1, Gs2, and Gs3
inline void stiefel_Gs03(double *Gs0, double *Gs1, double *Gs2, double *Gs3, double beta, double X)
{
#ifdef PRINT_UTILITY
    if (beta < beta_min)
        beta_min = beta;
    if (beta > beta_max)
        beta_max = beta;
    if (X < X_min)
        X_min = X;
    if (X > X_max)
        X_max = X;
#endif // PRINT_UTILITY

    double X2 = X * X;
    double z = X2 * beta;
    // Use a truncated Taylor series for the Stumpff functions
    // nmax = 11 as in the vectorized version
    const int nmax = 11;
    double invfact[12] = {
        1.0,                    // 0!
        1.0,                    // 1!
        0.5,                    // 2!
        0.16666666666666666,    // 3!
        0.041666666666666664,   // 4!
        0.008333333333333333,   // 5!
        0.001388888888888889,   // 6!
        0.0001984126984126984,  // 7!
        2.48015873015873e-05,   // 8!
        2.7557319223985893e-06, // 9!
        2.755731922398589e-07,  // 10!
        2.505210838544172e-08   // 11!
    };
    *Gs3 = invfact[nmax];
    *Gs2 = invfact[nmax - 1];
    for (int np = nmax - 2; np >= 3; np -= 2)
    {
        *Gs3 = -z * (*Gs3) + invfact[np];
        *Gs2 = -z * (*Gs2) + invfact[np - 1];
    }
    *Gs0 = -z * (*Gs2) + 1.0;
    *Gs3 = (*Gs3) * X;
    *Gs1 = -z * (*Gs3) + X;
    *Gs3 = (*Gs3) * X2;
    *Gs2 = (*Gs2) * X2;
}

// Scalar Halley step function for Kepler iteration
inline void halley_step(double *X, double beta, double r0, double eta0, double zeta0, double dt)
{
#ifdef PRINT_UTILITY
    if (r0 < halley_r0_min)
        halley_r0_min = r0;
    if (r0 > halley_r0_max)
        halley_r0_max = r0;
    if (eta0 < halley_eta0_min)
        halley_eta0_min = eta0;
    if (eta0 > halley_eta0_max)
        halley_eta0_max = eta0;
    if (zeta0 < halley_zeta0_min)
        halley_zeta0_min = zeta0;
    if (zeta0 > halley_zeta0_max)
        halley_zeta0_max = zeta0;
#endif // PRINT_UTILITY

    double Gs0, Gs1, Gs2, Gs3;
    stiefel_Gs03(&Gs0, &Gs1, &Gs2, &Gs3, beta, *X);
    double f = r0 * (*X) - dt;
    f += eta0 * Gs2;
    f += zeta0 * Gs3;
    double fp = r0 + eta0 * Gs1 + zeta0 * Gs2;
    double fpp = eta0 * Gs0 + zeta0 * Gs1;
    double denom = fp * fp * 16.0 - 20.0 * f * fpp;
    denom = sqrt(denom);
    denom = fp + denom;
    *X = ((*X) * denom - 5.0 * f) / denom;
}

// Scalar Stiefel function for Newton's method, returning Gs1, Gs2, and Gs3
inline void stiefel_Gs13(double *Gs1, double *Gs2, double *Gs3, double beta, double X)
{
#ifdef PRINT_UTILITY
    if (beta < beta_min)
        beta_min = beta;
    if (beta > beta_max)
        beta_max = beta;
    if (X < X_min)
        X_min = X;
    if (X > X_max)
        X_max = X;
#endif // PRINT_UTILITY

    double X2 = X * X;
    double z = X2 * beta;
    // Use a truncated Taylor series for the Stumpff functions
    // nmax = 19 as in the vectorized version
    const int nmax = 19;
    double invfact[20] = {
        1.0,                    // 0!
        1.0,                    // 1!
        0.5,                    // 2!
        0.16666666666666666,    // 3!
        0.041666666666666664,   // 4!
        0.008333333333333333,   // 5!
        0.001388888888888889,   // 6!
        0.0001984126984126984,  // 7!
        2.48015873015873e-05,   // 8!
        2.7557319223985893e-06, // 9!
        2.755731922398589e-07,  // 10!
        2.505210838544172e-08,  // 11!
        2.08767569878681e-09,   // 12!
        1.6059043836821613e-10, // 13!
        1.1470745597729725e-11, // 14!
        7.647163731819816e-13,  // 15!
        4.779477332387385e-14,  // 16!
        2.8114572543455206e-15, // 17!
        1.5619206968586225e-16, // 18!
        8.22063524662433e-18    // 19!
    };
    *Gs3 = invfact[nmax];
    *Gs2 = invfact[nmax - 1];
    for (int np = nmax - 2; np >= 3; np -= 2)
    {
        *Gs3 = -z * (*Gs3) + invfact[np];
        *Gs2 = -z * (*Gs2) + invfact[np - 1];
    }
    *Gs3 = (*Gs3) * X;
    *Gs1 = -z * (*Gs3) + X;
    *Gs3 = (*Gs3) * X2;
    *Gs2 = (*Gs2) * X2;
}

// Scalar Newton step function for Kepler iteration
inline void newton_step(double *X, double beta, double r0, double eta0, double zeta0, double dt)
{
#ifdef PRINT_UTILITY
    if (r0 < newton_r0_min)
        newton_r0_min = r0;
    if (r0 > newton_r0_max)
        newton_r0_max = r0;
    if (eta0 < newton_eta0_min)
        newton_eta0_min = eta0;
    if (eta0 > newton_eta0_max)
        newton_eta0_max = eta0;
    if (zeta0 < newton_zeta0_min)
        newton_zeta0_min = zeta0;
    if (zeta0 > newton_zeta0_max)
        newton_zeta0_max = zeta0;
#endif // PRINT_UTILITY

    double Gs1, Gs2, Gs3;
    stiefel_Gs13(&Gs1, &Gs2, &Gs3, beta, *X);
    double eta0Gs1zeta0Gs2 = eta0 * Gs1 + zeta0 * Gs2;
    double ri = 1.0 / (r0 + eta0Gs1zeta0Gs2);
    *X = (*X) * eta0Gs1zeta0Gs2 - eta0 * Gs2 - zeta0 * Gs3 + dt;
    *X = ri * (*X);
}

void whfast_kepler_step(double *x_vec, double *y_vec, double *z_vec,
                        double *vx_vec, double *vy_vec, double *vz_vec, double dt)
{
    for (int i = 0; i < N_PLANETS; i++)
    {
        double r2, r0, r01, v2;
        double beta, eta0, zeta0;
        double Gs1, Gs2, Gs3, eta0Gs1zeta0Gs2, ri;
        r2 = x_vec[i] * x_vec[i] + y_vec[i] * y_vec[i] + z_vec[i] * z_vec[i];
        r0 = sqrt(r2);
        r01 = 1.0 / r0;

        v2 = vx_vec[i] * vx_vec[i] + vy_vec[i] * vy_vec[i] + vz_vec[i] * vz_vec[i];
        beta = 2.0 * kConsts->M0 * r01 - v2;

        eta0 = x_vec[i] * vx_vec[i] + y_vec[i] * vy_vec[i] + z_vec[i] * vz_vec[i];
        zeta0 = -beta * r0 + kConsts->M0;

        // Initial guess
        double dtr0i = dt * r01;
        double X = dtr0i * eta0;
        X = -0.5 * X * r0 + 1.0;
        X = dtr0i * X;

        // Iterations
        halley_step(&X, beta, r0, eta0, zeta0, dt);
        halley_step(&X, beta, r0, eta0, zeta0, dt);
        newton_step(&X, beta, r0, eta0, zeta0, dt);

        stiefel_Gs13(&Gs1, &Gs2, &Gs3, beta, X);
        eta0Gs1zeta0Gs2 = eta0 * Gs1 + zeta0 * Gs2;
        ri = 1.0 / (r0 + eta0Gs1zeta0Gs2);

        // f and g function
        double nf = kConsts->M0 * Gs2; // negative f
        nf = nf * r01;

        double g = -kConsts->M0 * Gs3 + dt;
        double nfd = kConsts->M0 * Gs1; // negative fd
        nfd = nfd * r01;
        nfd = nfd * ri;

        double ngd = kConsts->M0 * Gs2; // negative gd
        ngd = ngd * ri;

        double nx = -nf * x_vec[i] + x_vec[i];
        nx = nx + g * vx_vec[i];

        double ny = -nf * y_vec[i] + y_vec[i];
        ny = ny + g * vy_vec[i];

        double nz = -nf * z_vec[i] + z_vec[i];
        nz = nz + g * vz_vec[i];

        vx_vec[i] = -ngd * vx_vec[i] + vx_vec[i];
        vx_vec[i] = -nfd * x_vec[i] + vx_vec[i];
        vy_vec[i] = -ngd * vy_vec[i] + vy_vec[i];
        vy_vec[i] = -nfd * y_vec[i] + vy_vec[i];
        vz_vec[i] = -ngd * vz_vec[i] + vz_vec[i];
        vz_vec[i] = -nfd * z_vec[i] + vz_vec[i];

        x_vec[i] = nx;
        y_vec[i] = ny;
        z_vec[i] = nz;
    }
}

void whfast_com_step(Body *com, double dt)
{
    for (int i = 0; i < 3; i++)
        com->pos[i] += dt * com->vel[i];
}

void whfast_drift_step(double *x_vec, double *y_vec, double *z_vec,
                       double *vx_vec, double *vy_vec, double *vz_vec, Body *com, double dt)
{
    // We first do the kepler step, then the com step.

    whfast_kepler_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, dt);

    whfast_com_step(com, dt);
}

// Performs one complete jump step (scalar version)
void whfast_jump_step(double *x_vec, double *y_vec, double *z_vec,
                      double *vx_vec, double *vy_vec, double *vz_vec,
                      double *m_vec, double dt)
{
    // pf = dt / M0
    double pf = dt / kConsts->M0;
    double sumx = 0.0, sumy = 0.0, sumz = 0.0;
    for (int i = 0; i < N_PLANETS; i++)
    {
        sumx += m_vec[i] * vx_vec[i];
        sumy += m_vec[i] * vy_vec[i];
        sumz += m_vec[i] * vz_vec[i];
    }
    for (int i = 0; i < N_PLANETS; i++)
    {
        x_vec[i] += sumx * pf;
        y_vec[i] += sumy * pf;
        z_vec[i] += sumz * pf;
    }
}

// Scalar helper functions for the interaction step
static inline double gravity_prefactor_one(double dx, double dy, double dz)
{
    double r2 = dx * dx + dy * dy + dz * dz;
    double r = sqrt(r2);
    return 1.0 / (r2 * r);
}
static inline double gravity_prefactor(double m, double dx, double dy, double dz)
{
    double r2 = dx * dx + dy * dy + dz * dz;
    double r = sqrt(r2);
    return m / (r2 * r);
}

// Performs one full interaction step (scalar version)
void whfast_interaction_step(double *x_vec, double *y_vec, double *z_vec,
                             double *vx_vec, double *vy_vec, double *vz_vec,
                             double *m_vec, double dt)
{

    double star_dvx = 0.0, star_dvy = 0.0, star_dvz = 0.0;
    for (int i = 0; i < N_PLANETS; ++i)
    {
        double xi = x_vec[i], yi = y_vec[i], zi = z_vec[i];
        double vxi = vx_vec[i], vyi = vy_vec[i], vzi = vz_vec[i];
        double r2 = xi * xi + yi * yi + zi * zi;
        double r4 = r2 * r2;
        double invr4 = kConsts->gr_prefac[i] * dt / r4;
        double dvx_i = invr4 * xi;
        double dvy_i = invr4 * yi;
        double dvz_i = invr4 * zi;
        vxi -= dvx_i;
        vyi -= dvy_i;
        vzi -= dvz_i;
        star_dvx += kConsts->gr_prefac2[i] * dvx_i;
        star_dvy += kConsts->gr_prefac2[i] * dvy_i;
        star_dvz += kConsts->gr_prefac2[i] * dvz_i;
        for (int j = 0; j < N_PLANETS; ++j)
        {
            if (i == j)
                continue;
            double dx = xi - x_vec[j];
            double dy = yi - y_vec[j];
            double dz = zi - z_vec[j];
            double pref = gravity_prefactor(m_vec[j], dx, dy, dz) * dt;
            vxi -= pref * dx;
            vyi -= pref * dy;
            vzi -= pref * dz;
        }
        vx_vec[i] = vxi;
        vy_vec[i] = vyi;
        vz_vec[i] = vzi;
    }
    // Apply back-reaction onto star to all planets
    for (int i = 0; i < N_PLANETS; ++i)
    {
        vx_vec[i] -= star_dvx;
        vy_vec[i] -= star_dvy;
        vz_vec[i] -= star_dvz;
    }
}

void whfast_kernel(double *x_vec, double *y_vec, double *z_vec,
                   double *vx_vec, double *vy_vec, double *vz_vec,
                   double *m_vec, Body *com, double dt, long Nint)
{
    // Call the main integration routine
    // This is where the actual integration happens

    // Perform the initial half-drift step
    whfast_drift_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, com, dt / 2.0);

    for (int i = 0; i < Nint - 1; i++)
    {
        // Perform the jump step (first half)
        whfast_jump_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt / 2.0);

        // Perform the interaction step
        whfast_interaction_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt);

        // Perform the jump step (second half)
        whfast_jump_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt / 2.0);

        // Perform the drift step
        whfast_drift_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, com, dt);
    }

    // Final iteration happens outside loop to avoid branching

    // Perform the jump step (first half)
    whfast_jump_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt / 2.0);

    // Perform the interaction step
    whfast_interaction_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt);

    // Perform the jump step (second half)
    whfast_jump_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt / 2.0);

    // Perform the final half-drift step
    whfast_drift_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, com, dt / 2.0);
}