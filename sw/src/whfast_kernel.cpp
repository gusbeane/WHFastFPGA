#include <stdio.h>
#include "whfast_kernel.h"
#include "whfast_constants.h"
#include <cmath>
#include <limits>

// Convert AVX512 vectors to double arrays
#define CONVERT_PLANETS_AVX512_TO_DOUBLE() \
do { \
    _mm512_storeu_pd(x_vec_, *x_vec); \
    _mm512_storeu_pd(y_vec_, *y_vec); \
    _mm512_storeu_pd(z_vec_, *z_vec); \
    _mm512_storeu_pd(vx_vec_, *vx_vec); \
    _mm512_storeu_pd(vy_vec_, *vy_vec); \
    _mm512_storeu_pd(vz_vec_, *vz_vec); \
    _mm512_storeu_pd(m_vec_, m_vec); \
} while(0)

// Convert double arrays back to AVX512 vectors
#define CONVERT_PLANETS_DOUBLE_TO_AVX512() \
do { \
    *x_vec = _mm512_loadu_pd(x_vec_); \
    *y_vec = _mm512_loadu_pd(y_vec_); \
    *z_vec = _mm512_loadu_pd(z_vec_); \
    *vx_vec = _mm512_loadu_pd(vx_vec_); \
    *vy_vec = _mm512_loadu_pd(vy_vec_); \
    *vz_vec = _mm512_loadu_pd(vz_vec_); \
    m_vec = _mm512_loadu_pd(m_vec_); \
} while(0)

// Scalar Stiefel function for Halley's method, returning Gs0, Gs1, Gs2, and Gs3
inline void stiefel_Gs03(double* Gs0, double* Gs1, double* Gs2, double* Gs3, double beta, double X) {
#ifdef PRINT_UTILITY
    if (beta < beta_min) beta_min = beta;
    if (beta > beta_max) beta_max = beta;
    if (X < X_min) X_min = X;
    if (X > X_max) X_max = X;
#endif // PRINT_UTILITY

    double X2 = X * X;
    double z = X2 * beta;
    // Use a truncated Taylor series for the Stumpff functions
    // nmax = 11 as in the vectorized version
    const int nmax = 11;
    double invfact[12] = {
        1.0, // 0!
        1.0, // 1!
        0.5, // 2!
        0.16666666666666666, // 3!
        0.041666666666666664, // 4!
        0.008333333333333333, // 5!
        0.001388888888888889, // 6!
        0.0001984126984126984, // 7!
        2.48015873015873e-05, // 8!
        2.7557319223985893e-06, // 9!
        2.755731922398589e-07, // 10!
        2.505210838544172e-08 // 11!
    };
    *Gs3 = invfact[nmax];
    *Gs2 = invfact[nmax-1];
    for (int np = nmax-2; np >= 3; np -= 2) {
        *Gs3 = -z * (*Gs3) + invfact[np];
        *Gs2 = -z * (*Gs2) + invfact[np-1];
    }
    *Gs0 = -z * (*Gs2) + 1.0;
    *Gs3 = (*Gs3) * X;
    *Gs1 = -z * (*Gs3) + X;
    *Gs3 = (*Gs3) * X2;
    *Gs2 = (*Gs2) * X2;
}

// Scalar Halley step function for Kepler iteration
static inline void halley_step(double* X, double beta, double r0, double eta0, double zeta0, double dt) {
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
inline void stiefel_Gs13(double* Gs1, double* Gs2, double* Gs3, double beta, double X) {
    #ifdef PRINT_UTILITY
    if (beta < beta_min) beta_min = beta;
    if (beta > beta_max) beta_max = beta;
    if (X < X_min) X_min = X;
    if (X > X_max) X_max = X;
    #endif // PRINT_UTILITY

    double X2 = X * X;
    double z = X2 * beta;
    // Use a truncated Taylor series for the Stumpff functions
    // nmax = 19 as in the vectorized version
    const int nmax = 19;
    double invfact[20] = {
        1.0, // 0!
        1.0, // 1!
        0.5, // 2!
        0.16666666666666666, // 3!
        0.041666666666666664, // 4!
        0.008333333333333333, // 5!
        0.001388888888888889, // 6!
        0.0001984126984126984, // 7!
        2.48015873015873e-05, // 8!
        2.7557319223985893e-06, // 9!
        2.755731922398589e-07, // 10!
        2.505210838544172e-08, // 11!
        2.08767569878681e-09, // 12!
        1.6059043836821613e-10, // 13!
        1.1470745597729725e-11, // 14!
        7.647163731819816e-13, // 15!
        4.779477332387385e-14, // 16!
        2.8114572543455206e-15, // 17!
        1.5619206968586225e-16, // 18!
        8.22063524662433e-18 // 19!
    };
    *Gs3 = invfact[nmax];
    *Gs2 = invfact[nmax-1];
    for (int np = nmax-2; np >= 3; np -= 2) {
        *Gs3 = -z * (*Gs3) + invfact[np];
        *Gs2 = -z * (*Gs2) + invfact[np-1];
    }
    *Gs3 = (*Gs3) * X;
    *Gs1 = -z * (*Gs3) + X;
    *Gs3 = (*Gs3) * X2;
    *Gs2 = (*Gs2) * X2;
}

// Scalar Newton step function for Kepler iteration
static inline void newton_step(double* X, double beta, double r0, double eta0, double zeta0, double dt) {
    double Gs1, Gs2, Gs3;
    stiefel_Gs13(&Gs1, &Gs2, &Gs3, beta, *X);
    double eta0Gs1zeta0Gs2 = eta0 * Gs1 + zeta0 * Gs2;
    double ri = 1.0 / (r0 + eta0Gs1zeta0Gs2);
    *X = (*X) * eta0Gs1zeta0Gs2 - eta0 * Gs2 - zeta0 * Gs3 + dt;
    *X = ri * (*X);
}

void whfast_kepler_step(double *x_vec,  double *y_vec,  double *z_vec,
                          double *vx_vec, double *vy_vec, double *vz_vec,
                          double *m_vec, double dt)
{
    for(int i=0; i<N_PLANETS; i++)
    {   
        double r2, r0, r01, v2;
        double beta, eta0, zeta0;
        double Gs1, Gs2, Gs3, eta0Gs1zeta0Gs2, ri;
        r2 = x_vec[i]*x_vec[i] + y_vec[i]*y_vec[i] + z_vec[i]*z_vec[i];
        r0 = sqrt(r2);
        r01 = 1.0/r0;

        v2 = vx_vec[i]*vx_vec[i] + vy_vec[i]*vy_vec[i] + vz_vec[i]*vz_vec[i];
        beta = 2.0*kConsts->M0*r01 - v2;

        eta0 = x_vec[i]*vx_vec[i] + y_vec[i]*vy_vec[i] + z_vec[i]*vz_vec[i];
        zeta0 = -beta*r0 + kConsts->M0;
        
        // Initial guess
        double dtr0i = dt*r01;
        double X = dtr0i*eta0;
        X = -0.5*X*r0 + 1.0;
        X = dtr0i*X;

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

        double nx = -nf*x_vec[i] + x_vec[i];
        nx = nx + g*vx_vec[i];

        double ny = -nf*y_vec[i] + y_vec[i];
        ny = ny + g*vy_vec[i];

        double nz = -nf*z_vec[i] + z_vec[i];
        nz = nz + g*vz_vec[i];

        vx_vec[i] = -ngd*vx_vec[i] + vx_vec[i];
        vx_vec[i] = -nfd*x_vec[i] + vx_vec[i];
        vy_vec[i] = -ngd*vy_vec[i] + vy_vec[i];
        vy_vec[i] = -nfd*y_vec[i] + vy_vec[i];
        vz_vec[i] = -ngd*vz_vec[i] + vz_vec[i];
        vz_vec[i] = -nfd*z_vec[i] + vz_vec[i];
        
        x_vec[i] = nx;
        y_vec[i] = ny;
        z_vec[i] = nz;
    }
}

void whfast_com_step(Body *com, double dt)
{
    for(int i=0; i<3; i++)
        com->pos[i] += dt * com->vel[i];
}

void whfast_drift_step(__m512d *x_vec,  __m512d *y_vec,  __m512d *z_vec,
    __m512d *vx_vec, __m512d *vy_vec, __m512d *vz_vec,
    __m512d m_vec, Body *com, double dt)
{
// We first do the kepler step, then the com step.
__m512d xload = _mm512_loadu_pd(x_vec);

double x_vec_[N_PLANETS], y_vec_[N_PLANETS], z_vec_[N_PLANETS];
double vx_vec_[N_PLANETS], vy_vec_[N_PLANETS], vz_vec_[N_PLANETS];
double m_vec_[N_PLANETS];

    CONVERT_PLANETS_AVX512_TO_DOUBLE();
    whfast_kepler_step(x_vec_, y_vec_, z_vec_, vx_vec_, vy_vec_, vz_vec_, m_vec_, dt);
    CONVERT_PLANETS_DOUBLE_TO_AVX512();

whfast_com_step(com, dt);
}

// Performs one complete jump step
void whfast_jump_step(__m512d *x_vec, __m512d *y_vec, __m512d *z_vec,
                         __m512d *vx_vec, __m512d *vy_vec, __m512d *vz_vec,
                         __m512d m_vec, double dt)
{
    __m512d pf512 = _mm512_set1_pd(dt / kConsts->M0);

    __m512d sumx = _mm512_mul_pd(m_vec, *vx_vec);
    __m512d sumy = _mm512_mul_pd(m_vec, *vy_vec);
    __m512d sumz = _mm512_mul_pd(m_vec, *vz_vec);

    sumx = _mm512_add_pd(_mm512_shuffle_pd(sumx, sumx, 0x55), sumx); // Swapping neighbouring elements
    sumx = _mm512_add_pd(_mm512_permutex_pd(sumx, _MM_PERM_ABCD), sumx);
    sumx = _mm512_add_pd(_mm512_shuffle_f64x2(sumx, sumx, 78), sumx); // 78 is _MM_SHUFFLE(1,0,3,2), changed for icx

    sumy = _mm512_add_pd(_mm512_shuffle_pd(sumy, sumy, 0x55), sumy);
    sumy = _mm512_add_pd(_mm512_permutex_pd(sumy, _MM_PERM_ABCD), sumy);
    sumy = _mm512_add_pd(_mm512_shuffle_f64x2(sumy, sumy, 78), sumy);

    sumz = _mm512_add_pd(_mm512_shuffle_pd(sumz, sumz, 0x55), sumz);
    sumz = _mm512_add_pd(_mm512_permutex_pd(sumz, _MM_PERM_ABCD), sumz);
    sumz = _mm512_add_pd(_mm512_shuffle_f64x2(sumz, sumz, 78), sumz);

    *x_vec = _mm512_fmadd_pd(sumx, pf512, *x_vec);
    *y_vec = _mm512_fmadd_pd(sumy, pf512, *y_vec);
    *z_vec = _mm512_fmadd_pd(sumz, pf512, *z_vec);
}

// Helper functions for the interaction step
static __m512d inline gravity_prefactor_avx512_one(__m512d dx, __m512d dy, __m512d dz)
{
    __m512d r2 = _mm512_mul_pd(dx, dx);
    r2 = _mm512_fmadd_pd(dy, dy, r2);
    r2 = _mm512_fmadd_pd(dz, dz, r2);
    const __m512d r = _mm512_sqrt_pd(r2);
    const __m512d r3 = _mm512_mul_pd(r, r2);
    return _mm512_div_pd(kConsts->one, r3);
}

static __m512d inline gravity_prefactor_avx512(__m512d m, __m512d dx, __m512d dy, __m512d dz)
{
    __m512d r2 = _mm512_mul_pd(dx, dx);
    r2 = _mm512_fmadd_pd(dy, dy, r2);
    r2 = _mm512_fmadd_pd(dz, dz, r2);
    const __m512d r = _mm512_sqrt_pd(r2);
    const __m512d r3 = _mm512_mul_pd(r, r2);
    return _mm512_div_pd(m, r3);
}

// Performs one full interaction step
void whfast_interaction_step(__m512d *x_vec, __m512d *y_vec, __m512d *z_vec,
                                __m512d *vx_vec, __m512d *vy_vec, __m512d *vz_vec,
                                __m512d m_vec, double dt)
{
    __m512d x_j = *x_vec;
    __m512d y_j = *y_vec;
    __m512d z_j = *z_vec;
    __m512d dt512 = _mm512_set1_pd(dt);

    // General relativistic corrections
    if (true)
    {
        __m512d r2 = _mm512_mul_pd(x_j, x_j);
        r2 = _mm512_fmadd_pd(y_j, y_j, r2);
        r2 = _mm512_fmadd_pd(z_j, z_j, r2);
        const __m512d r4 = _mm512_mul_pd(r2, r2);
        __m512d prefac = _mm512_div_pd(kConsts->gr_prefac, r4);
        prefac = _mm512_mul_pd(prefac, dt512);
        __m512d dvx = _mm512_mul_pd(prefac, x_j);
        __m512d dvy = _mm512_mul_pd(prefac, y_j);
        __m512d dvz = _mm512_mul_pd(prefac, z_j);
        *vx_vec = _mm512_sub_pd(*vx_vec, dvx);
        *vy_vec = _mm512_sub_pd(*vy_vec, dvy);
        *vz_vec = _mm512_sub_pd(*vz_vec, dvz);

        // Calculate back reaction onto star and apply them to planets (heliocentric)
        dvx = _mm512_mul_pd(kConsts->gr_prefac2, dvx);
        dvy = _mm512_mul_pd(kConsts->gr_prefac2, dvy);
        dvz = _mm512_mul_pd(kConsts->gr_prefac2, dvz);

        dvx = _mm512_add_pd(_mm512_shuffle_pd(dvx, dvx, 0x55), dvx); // Swapping neighbouring elements
        dvx = _mm512_add_pd(_mm512_permutex_pd(dvx, _MM_PERM_ABCD), dvx);
        dvx = _mm512_add_pd(_mm512_shuffle_f64x2(dvx, dvx, 78), dvx);
        dvy = _mm512_add_pd(_mm512_shuffle_pd(dvy, dvy, 0x55), dvy);
        dvy = _mm512_add_pd(_mm512_permutex_pd(dvy, _MM_PERM_ABCD), dvy);
        dvy = _mm512_add_pd(_mm512_shuffle_f64x2(dvy, dvy, 78), dvy);
        dvz = _mm512_add_pd(_mm512_shuffle_pd(dvz, dvz, 0x55), dvz);
        dvz = _mm512_add_pd(_mm512_permutex_pd(dvz, _MM_PERM_ABCD), dvz);
        dvz = _mm512_add_pd(_mm512_shuffle_f64x2(dvz, dvz, 78), dvz);

        *vx_vec = _mm512_sub_pd(*vx_vec, dvx);
        *vy_vec = _mm512_sub_pd(*vy_vec, dvy);
        *vz_vec = _mm512_sub_pd(*vz_vec, dvz);
    }

    __m512d m_j = _mm512_mul_pd(m_vec, dt512);
    __m512d m_j_01234567 = m_j;

    {
        x_j = _mm512_permutex_pd(x_j, _MM_PERM_BACD); // within 256
        y_j = _mm512_permutex_pd(y_j, _MM_PERM_BACD);
        z_j = _mm512_permutex_pd(z_j, _MM_PERM_BACD);
        m_j = _mm512_permutex_pd(m_j, _MM_PERM_BACD);
        __m512d dx_j = _mm512_sub_pd(*x_vec, x_j);
        __m512d dy_j = _mm512_sub_pd(*y_vec, y_j);
        __m512d dz_j = _mm512_sub_pd(*z_vec, z_j);
        __m512d prefact = gravity_prefactor_avx512_one(dx_j, dy_j, dz_j);

        // 0123 4567
        // 3201 7645
        __m512d prefact1 = _mm512_mul_pd(prefact, m_j);
        *vx_vec = _mm512_fnmadd_pd(prefact1, dx_j, *vx_vec);
        *vy_vec = _mm512_fnmadd_pd(prefact1, dy_j, *vy_vec);
        *vz_vec = _mm512_fnmadd_pd(prefact1, dz_j, *vz_vec);

        dx_j = _mm512_permutex_pd(dx_j, _MM_PERM_ABDC); // within 256
        dy_j = _mm512_permutex_pd(dy_j, _MM_PERM_ABDC);
        dz_j = _mm512_permutex_pd(dz_j, _MM_PERM_ABDC);
        prefact = _mm512_permutex_pd(prefact, _MM_PERM_ABDC);
        m_j = _mm512_permute_pd(m_j, 0x55); // within 128

        // 0123 4567
        // 2310 6754
        __m512d prefact2 = _mm512_mul_pd(prefact, m_j);
        *vx_vec = _mm512_fmadd_pd(prefact2, dx_j, *vx_vec);
        *vy_vec = _mm512_fmadd_pd(prefact2, dy_j, *vy_vec);
        *vz_vec = _mm512_fmadd_pd(prefact2, dz_j, *vz_vec);
    }
    {
        x_j = _mm512_permutex_pd(x_j, _MM_PERM_BACD); // within 256
        y_j = _mm512_permutex_pd(y_j, _MM_PERM_BACD);
        z_j = _mm512_permutex_pd(z_j, _MM_PERM_BACD);
        m_j = _mm512_permutex_pd(m_j, _MM_PERM_ABDC);

        const __m512d dx_j = _mm512_sub_pd(*x_vec, x_j);
        const __m512d dy_j = _mm512_sub_pd(*y_vec, y_j);
        const __m512d dz_j = _mm512_sub_pd(*z_vec, z_j);

        // 0123 4567
        // 1032 5476
        const __m512d prefact = gravity_prefactor_avx512(m_j, dx_j, dy_j, dz_j);
        *vx_vec = _mm512_fnmadd_pd(prefact, dx_j, *vx_vec);
        *vy_vec = _mm512_fnmadd_pd(prefact, dy_j, *vy_vec);
        *vz_vec = _mm512_fnmadd_pd(prefact, dz_j, *vz_vec);
    }

    // //////////////////////////////////////
    // 256 bit lane crossing
    // //////////////////////////////////////

    __m512d dvx; // delta vx for 4567 1230
    __m512d dvy;
    __m512d dvz;

    {
        x_j = _mm512_permutexvar_pd(kConsts->so1, x_j); // accros 512
        y_j = _mm512_permutexvar_pd(kConsts->so1, y_j);
        z_j = _mm512_permutexvar_pd(kConsts->so1, z_j);
        m_j = _mm512_permutexvar_pd(kConsts->so1, m_j);

        __m512d dx_j = _mm512_sub_pd(*x_vec, x_j);
        __m512d dy_j = _mm512_sub_pd(*y_vec, y_j);
        __m512d dz_j = _mm512_sub_pd(*z_vec, z_j);
        __m512d prefact = gravity_prefactor_avx512_one(dx_j, dy_j, dz_j);

        // 0123 4567
        // 4567 1230
        __m512d prefact1 = _mm512_mul_pd(prefact, m_j);
        *vx_vec = _mm512_fnmadd_pd(prefact1, dx_j, *vx_vec);
        *vy_vec = _mm512_fnmadd_pd(prefact1, dy_j, *vy_vec);
        *vz_vec = _mm512_fnmadd_pd(prefact1, dz_j, *vz_vec);

        // 4567 1230
        // 0123 4567
        prefact = _mm512_mul_pd(prefact, m_j_01234567);
        dvx = _mm512_mul_pd(prefact, dx_j);
        dvy = _mm512_mul_pd(prefact, dy_j);
        dvz = _mm512_mul_pd(prefact, dz_j);
    }

    {
        x_j = _mm512_permutex_pd(x_j, _MM_PERM_ADCB); // within 256
        y_j = _mm512_permutex_pd(y_j, _MM_PERM_ADCB);
        z_j = _mm512_permutex_pd(z_j, _MM_PERM_ADCB);
        m_j = _mm512_permutex_pd(m_j, _MM_PERM_ADCB);

        __m512d dx_j = _mm512_sub_pd(*x_vec, x_j);
        __m512d dy_j = _mm512_sub_pd(*y_vec, y_j);
        __m512d dz_j = _mm512_sub_pd(*z_vec, z_j);
        __m512d prefact = gravity_prefactor_avx512_one(dx_j, dy_j, dz_j);

        // 0123 4567
        // 5674 2301
        __m512d prefact1 = _mm512_mul_pd(prefact, m_j);
        *vx_vec = _mm512_fnmadd_pd(prefact1, dx_j, *vx_vec);
        *vy_vec = _mm512_fnmadd_pd(prefact1, dy_j, *vy_vec);
        *vz_vec = _mm512_fnmadd_pd(prefact1, dz_j, *vz_vec);

        dx_j = _mm512_permutex_pd(dx_j, _MM_PERM_CBAD); // within 256
        dy_j = _mm512_permutex_pd(dy_j, _MM_PERM_CBAD);
        dz_j = _mm512_permutex_pd(dz_j, _MM_PERM_CBAD);
        prefact = _mm512_permutex_pd(prefact, _MM_PERM_CBAD);
        m_j_01234567 = _mm512_permutex_pd(m_j_01234567, _MM_PERM_CBAD);

        // 4567 1230
        // 3012 7456
        prefact = _mm512_mul_pd(prefact, m_j_01234567);
        dvx = _mm512_fmadd_pd(prefact, dx_j, dvx);
        dvy = _mm512_fmadd_pd(prefact, dy_j, dvy);
        dvz = _mm512_fmadd_pd(prefact, dz_j, dvz);
    }

    // //////////////////////////////////////
    // 256 bit lane crossing for final add
    // //////////////////////////////////////

    {
        dvx = _mm512_permutexvar_pd(kConsts->so2, dvx); // across 512
        dvy = _mm512_permutexvar_pd(kConsts->so2, dvy);
        dvz = _mm512_permutexvar_pd(kConsts->so2, dvz);

        *vx_vec = _mm512_add_pd(dvx, *vx_vec);
        *vy_vec = _mm512_add_pd(dvy, *vy_vec);
        *vz_vec = _mm512_add_pd(dvz, *vz_vec);
    }
}

void whfast_kernel(__m512d *x_vec, __m512d *y_vec, __m512d *z_vec,
                      __m512d *vx_vec, __m512d *vy_vec, __m512d *vz_vec,
                      __m512d m_vec, Body *com, double dt, long Nint)
{
    // Call the main integration routine
    // This is where the actual integration happens

    // Perform the initial half-drift step
    whfast_drift_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, com, dt / 2.0);

    for (int i = 0; i < Nint - 1; i++)
    {
        // Perform the jump step (first half)
        whfast_jump_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt / 2.0);

        // Perform the interaction step
        whfast_interaction_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt);

        // Perform the jump step (second half)
        whfast_jump_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt / 2.0);

        // Perform the drift step
        whfast_drift_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, com, dt);
    }

    // Final iteration happens outside loop to avoid branching

    // Perform the jump step (first half)
    whfast_jump_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt / 2.0);

    // Perform the interaction step
    whfast_interaction_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt);

    // Perform the jump step (second half)
    whfast_jump_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt / 2.0);

    // Perform the final half-drift step
    whfast_drift_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, com, dt / 2.0);
}