#include <stdio.h>
#include "whfast512_kernel.h"
#include "whfast512_constants.h"

// Stiefel function for Newton's method, returning Gs1, Gs2, and Gs3
static void inline mm_stiefel_Gs13_avx512(__m512d * Gs1, __m512d * Gs2, __m512d * Gs3, __m512d beta, __m512d X){
    __m512d X2 = _mm512_mul_pd(X,X); 
    __m512d z = _mm512_mul_pd(X2,beta); 

    // stumpff_cs. Note: assuming n = 0
    const int nmax = 19;
   *Gs3 = kConsts->invfactorial512[nmax]; 
   *Gs2 = kConsts->invfactorial512[nmax-1]; 

    for(int np=nmax-2;np>=3;np-=2){
        *Gs3 = _mm512_fnmadd_pd(z, *Gs3, kConsts->invfactorial512[np]);
        *Gs2 = _mm512_fnmadd_pd(z, *Gs2, kConsts->invfactorial512[np-1]);
    }
    *Gs3 = _mm512_mul_pd(*Gs3,X); 
    *Gs1 = _mm512_fnmadd_pd(z, *Gs3, X);
    *Gs3 = _mm512_mul_pd(*Gs3,X2); 
    *Gs2 = _mm512_mul_pd(*Gs2,X2); 
};

// Stiefel function for Halley's method, returning Gs0, Gs1, Gs2, and Gs3
static void inline mm_stiefel_Gs03_avx512(__m512d * Gs0, __m512d * Gs1, __m512d * Gs2, __m512d * Gs3, __m512d beta, __m512d X){
    __m512d X2 = _mm512_mul_pd(X,X); 
    __m512d z = _mm512_mul_pd(X2,beta); 

    // stumpff_cs. Note: assuming n = 0
    const int nmax = 11; // Note: reduced! needs to be improved with mm_stiefel_Gs13_avx512 on last step(s)
   *Gs3 = kConsts->invfactorial512[nmax]; 
   *Gs2 = kConsts->invfactorial512[nmax-1]; 

    for(int np=nmax-2;np>=3;np-=2){
        *Gs3 = _mm512_fnmadd_pd(z, *Gs3, kConsts->invfactorial512[np]);
        *Gs2 = _mm512_fnmadd_pd(z, *Gs2, kConsts->invfactorial512[np-1]);
    }
    *Gs0 = _mm512_fnmadd_pd(z, *Gs2, kConsts->one);
    *Gs3 = _mm512_mul_pd(*Gs3,X); 
    *Gs1 = _mm512_fnmadd_pd(z, *Gs3, X);
    *Gs3 = _mm512_mul_pd(*Gs3,X2); 
    *Gs2 = _mm512_mul_pd(*Gs2,X2); 
};

// Newton step function for Kepler iteration
static void inline newton_step(__m512d* X, __m512d beta, __m512d r0, __m512d eta0, __m512d zeta0, __m512d _dt) {
    __m512d Gs1, Gs2, Gs3;
    mm_stiefel_Gs13_avx512(&Gs1, &Gs2, &Gs3, beta, *X);
    
    __m512d eta0Gs1zeta0Gs2 = _mm512_mul_pd(eta0, Gs1);
    eta0Gs1zeta0Gs2 = _mm512_fmadd_pd(zeta0, Gs2, eta0Gs1zeta0Gs2);
    
    __m512d ri = _mm512_add_pd(r0, eta0Gs1zeta0Gs2);
    ri = _mm512_div_pd(kConsts->one, ri);
    
    *X = _mm512_mul_pd(*X, eta0Gs1zeta0Gs2);
    *X = _mm512_fnmadd_pd(eta0, Gs2, *X);
    *X = _mm512_fnmadd_pd(zeta0, Gs3, *X);
    *X = _mm512_add_pd(_dt, *X); // Using dt parameter
    *X = _mm512_mul_pd(ri, *X);
}

// Halley step function for Kepler iteration
static void inline halley_step(__m512d* X, __m512d beta, __m512d r0, __m512d eta0, __m512d zeta0, __m512d _dt) {
    __m512d Gs0, Gs1, Gs2, Gs3;
    mm_stiefel_Gs03_avx512(&Gs0, &Gs1, &Gs2, &Gs3, beta, *X);
    
    __m512d f = _mm512_fmsub_pd(r0, *X, _dt);
    f = _mm512_fmadd_pd(eta0, Gs2, f);
    f = _mm512_fmadd_pd(zeta0, Gs3, f);
    
    __m512d fp = _mm512_fmadd_pd(eta0, Gs1, r0);
    fp = _mm512_fmadd_pd(zeta0, Gs2, fp);
    
    __m512d fpp = _mm512_mul_pd(eta0, Gs0);
    fpp = _mm512_fmadd_pd(zeta0, Gs1, fpp);
    
    __m512d denom = _mm512_mul_pd(fp, fp);
    denom = _mm512_mul_pd(denom, kConsts->sixteen);
    
    denom = _mm512_fnmadd_pd(_mm512_mul_pd(f, fpp), kConsts->twenty, denom);
    /* not included: _mm512_abs_pd(denom) */
    denom = _mm512_sqrt_pd(denom);
    denom = _mm512_add_pd(fp, denom);
    
    *X = _mm512_fmsub_pd(*X, denom, _mm512_mul_pd(f, kConsts->five));
    *X = _mm512_div_pd(*X, denom);
}

void whfast512_kepler_step(__m512d *x_vec,  __m512d *y_vec,  __m512d *z_vec,
                          __m512d *vx_vec, __m512d *vy_vec, __m512d *vz_vec,
                          __m512d m_vec, double dt)
{
// Performs one full Kepler step
        __m512d _dt = _mm512_set1_pd(dt); 
            
        __m512d r2 = _mm512_mul_pd(*x_vec, *x_vec);
        r2 = _mm512_fmadd_pd(*y_vec, *y_vec, r2);
        r2 = _mm512_fmadd_pd(*z_vec, *z_vec, r2);
        __m512d r0 = _mm512_sqrt_pd(r2);
        __m512d r0i = _mm512_div_pd(kConsts->one,r0);
    
        __m512d v2 = _mm512_mul_pd(*vx_vec, *vx_vec);
        v2 = _mm512_fmadd_pd(*vy_vec, *vy_vec, v2);
        v2 = _mm512_fmadd_pd(*vz_vec, *vz_vec, v2);
        
        __m512d beta = _mm512_mul_pd(kConsts->two, kConsts->_M);
        beta = _mm512_fmsub_pd(beta, r0i, v2);
    
        __m512d eta0 = _mm512_mul_pd(*x_vec, *vx_vec);
        eta0 = _mm512_fmadd_pd(*y_vec, *vy_vec, eta0);
        eta0 = _mm512_fmadd_pd(*z_vec, *vz_vec, eta0);
    
        __m512d zeta0 = _mm512_fnmadd_pd(beta, r0, kConsts->_M);
    
        __m512d Gs1;
        __m512d Gs2;
        __m512d Gs3;
        __m512d eta0Gs1zeta0Gs2; 
        __m512d ri; 
    
        // Initial guess
        __m512d dtr0i = _mm512_mul_pd(_dt,r0i);
        __m512d X = _mm512_mul_pd(dtr0i,eta0);
        X = _mm512_mul_pd(X, kConsts->half);
        X = _mm512_fnmadd_pd(X,r0i,kConsts->one);
        X = _mm512_mul_pd(dtr0i,X);
    
        // Iterations
        halley_step(&X, beta, r0, eta0, zeta0, _dt);
        halley_step(&X, beta, r0, eta0, zeta0, _dt);
        newton_step(&X, beta, r0, eta0, zeta0, _dt);
        
        // Final Newton step (note: X not needed after this) 
        mm_stiefel_Gs13_avx512(&Gs1, &Gs2, &Gs3, beta, X);
        eta0Gs1zeta0Gs2 = _mm512_mul_pd(eta0, Gs1); 
        eta0Gs1zeta0Gs2 = _mm512_fmadd_pd(zeta0,Gs2, eta0Gs1zeta0Gs2); 
        ri = _mm512_add_pd(r0, eta0Gs1zeta0Gs2); 
        ri = _mm512_div_pd(kConsts->one, ri);
    
        // f and g function
    
        __m512d nf = _mm512_mul_pd(kConsts->_M,Gs2); //negative f
        nf = _mm512_mul_pd(nf,r0i);

        printf("nf[0] = %f\n", nf[0]);
    
        __m512d g = _mm512_fnmadd_pd(kConsts->_M, Gs3, _dt);

        printf("g[0] = %f\n", g[0]);
    
        __m512d nfd = _mm512_mul_pd(kConsts->_M, Gs1); // negative fd
        nfd = _mm512_mul_pd(nfd, r0i);
        nfd = _mm512_mul_pd(nfd, ri);

        printf("nfd[0] = %f\n", nfd[0]);
    
        __m512d ngd = _mm512_mul_pd(kConsts->_M, Gs2); // negative gd
        ngd = _mm512_mul_pd(ngd, ri);

        printf("ngd[0] = %f\n", ngd[0]);
    
        __m512d nx = _mm512_fnmadd_pd(nf, *x_vec, *x_vec);
        nx = _mm512_fmadd_pd(g, *vx_vec, nx);
        __m512d ny = _mm512_fnmadd_pd(nf, *y_vec, *y_vec);
        ny = _mm512_fmadd_pd(g, *vy_vec, ny);
        __m512d nz = _mm512_fnmadd_pd(nf, *z_vec, *z_vec);
        nz = _mm512_fmadd_pd(g, *vz_vec, nz);
    
        *vx_vec = _mm512_fnmadd_pd(ngd, *vx_vec, *vx_vec);
        *vx_vec = _mm512_fnmadd_pd(nfd, *x_vec, *vx_vec);
        *vy_vec = _mm512_fnmadd_pd(ngd, *vy_vec, *vy_vec);
        *vy_vec = _mm512_fnmadd_pd(nfd, *y_vec, *vy_vec);
        *vz_vec = _mm512_fnmadd_pd(ngd, *vz_vec, *vz_vec);
        *vz_vec = _mm512_fnmadd_pd(nfd, *z_vec, *vz_vec);
    
        *x_vec = nx;
        *y_vec = ny;
        *z_vec = nz;
}

void whfast512_com_step(Body *com, double dt)
{
    for(int i=0; i<3; i++)
        com->pos[i] += dt * com->vel[i];
}

void whfast512_drift_step(__m512d *x_vec,  __m512d *y_vec,  __m512d *z_vec,
    __m512d *vx_vec, __m512d *vy_vec, __m512d *vz_vec,
    __m512d m_vec, Body *com, double dt)
{
// We first do the kepler step, then the com step.
__m512d xload = _mm512_loadu_pd(x_vec);

whfast512_kepler_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt);

whfast512_com_step(com, dt);
}

// void whfast512_jump_step(__m512d *x_vec,  __m512d *y_vec,  __m512d *z_vec,
//     __m512d *vx_vec, __m512d *vy_vec, __m512d *vz_vec,
//     __m512d m_vec, double dt)
// {

// }

// Performs one complete jump step
void whfast512_jump_step(__m512d *x_vec, __m512d *y_vec, __m512d *z_vec,
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
void whfast512_interaction_step(__m512d *x_vec, __m512d *y_vec, __m512d *z_vec,
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

void whfast512_kernel(__m512d *x_vec, __m512d *y_vec, __m512d *z_vec,
                      __m512d *vx_vec, __m512d *vy_vec, __m512d *vz_vec,
                      __m512d m_vec, Body *com, double dt, long Nint)
{
    // Call the main integration routine
    // This is where the actual integration happens

    // We first do a half drift step
    whfast512_drift_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, com, dt/2.);

        // Now we enter main loop
    for (int i = 0; i < Nint - 1; i++)
    {
        // do first half-jump step
        whfast512_jump_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt / 2.);

        // do interaction step
        whfast512_interaction_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt);

        // do second half-jump step
        whfast512_jump_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt / 2.);

        // do drift step
        whfast512_drift_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, com, dt);
    }

    // // Do final step and synchronize
    // whfast512_jump_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt);
    // whfast512_interaction_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt);
    whfast512_drift_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, com, dt/2.);
}