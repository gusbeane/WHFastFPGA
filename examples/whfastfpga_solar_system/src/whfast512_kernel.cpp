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
static void inline newton_step(__m512d* X, __m512d beta, __m512d r0, __m512d eta0, __m512d zeta0) {
    __m512d Gs1, Gs2, Gs3;
    mm_stiefel_Gs13_avx512(&Gs1, &Gs2, &Gs3, beta, *X);
    
    __m512d eta0Gs1zeta0Gs2 = _mm512_mul_pd(eta0, Gs1);
    eta0Gs1zeta0Gs2 = _mm512_fmadd_pd(zeta0, Gs2, eta0Gs1zeta0Gs2);
    
    __m512d ri = _mm512_add_pd(r0, eta0Gs1zeta0Gs2);
    ri = _mm512_div_pd(kConsts->one, ri);
    
    *X = _mm512_mul_pd(*X, eta0Gs1zeta0Gs2);
    *X = _mm512_fnmadd_pd(eta0, Gs2, *X);
    *X = _mm512_fnmadd_pd(zeta0, Gs3, *X);
    *X = _mm512_add_pd(_mm512_set1_pd(1.0), *X); // Using dt parameter
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
        newton_step(&X, beta, r0, eta0, zeta0);
        
        // Final Newton step (note: X not needed after this) 
        mm_stiefel_Gs13_avx512(&Gs1, &Gs2, &Gs3, beta, X);
        eta0Gs1zeta0Gs2 = _mm512_mul_pd(eta0, Gs1); 
        eta0Gs1zeta0Gs2 = _mm512_fmadd_pd(zeta0,Gs2, eta0Gs1zeta0Gs2); 
        ri = _mm512_add_pd(r0, eta0Gs1zeta0Gs2); 
        ri = _mm512_div_pd(kConsts->one, ri);
    
        // f and g function
    
        __m512d nf = _mm512_mul_pd(kConsts->_M,Gs2); //negative f
        nf = _mm512_mul_pd(nf,r0i); 
    
        __m512d g = _mm512_fnmadd_pd(kConsts->_M, Gs3, _dt);
    
        __m512d nfd = _mm512_mul_pd(kConsts->_M, Gs1); // negative fd
        nfd = _mm512_mul_pd(nfd, r0i);
        nfd = _mm512_mul_pd(nfd, ri);
    
        __m512d ngd = _mm512_mul_pd(kConsts->_M, Gs2); // negative gd
        ngd = _mm512_mul_pd(ngd, ri);
    
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
    whfast512_kepler_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt);
    whfast512_com_step(com, dt);
}

// void whfast512_jump_step(__m512d *x_vec,  __m512d *y_vec,  __m512d *z_vec,
//     __m512d *vx_vec, __m512d *vy_vec, __m512d *vz_vec,
//     __m512d m_vec, double dt)
// {

// }

// void whfast512_interaction_step(__m512d *x_vec,  __m512d *y_vec,  __m512d *z_vec,
//     __m512d *vx_vec, __m512d *vy_vec, __m512d *vz_vec,
//     __m512d m_vec, double dt)
// {

// }

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
        // do jump step. this needs to be split when we add gr
        // whfast512_jump_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt);

        // do interaction step
        // whfast512_interaction_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt);

        // do drift step
        whfast512_drift_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, com, dt);
    }

    // // Do final step and synchronize
    // whfast512_jump_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt);
    // whfast512_interaction_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, dt);
    whfast512_drift_step(x_vec, y_vec, z_vec, vx_vec, vy_vec, vz_vec, m_vec, com, dt/2.);
}