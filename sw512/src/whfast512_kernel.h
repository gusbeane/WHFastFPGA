#ifndef WHFAST512_KERNEL_H
#define WHFAST512_KERNEL_H

#include <immintrin.h>
#include "util.h"
#include "whfast512_constants.h"

void whfast512_kernel(__m512d *x_vec, __m512d *y_vec, __m512d *z_vec,
                      __m512d *vx_vec, __m512d *vy_vec, __m512d *vz_vec,
                      __m512d m_vec, Body *com, double dt, long Nint);

void whfast512_kepler_step(__m512d *x_vec, __m512d *y_vec, __m512d *z_vec,
                           __m512d *vx_vec, __m512d *vy_vec, __m512d *vz_vec,
                           __m512d m_vec, double dt);

void whfast512_com_step(Body *com, double dt);

void whfast512_drift_step(__m512d *x_vec, __m512d *y_vec, __m512d *z_vec,
                          __m512d *vx_vec, __m512d *vy_vec, __m512d *vz_vec,
                          __m512d m_vec, Body *com, double dt);

void whfast512_jump_step(__m512d *x_vec, __m512d *y_vec, __m512d *z_vec,
                         __m512d *vx_vec, __m512d *vy_vec, __m512d *vz_vec,
                         __m512d m_vec, double dt);

void whfast512_interaction_step(__m512d *x_vec,  __m512d *y_vec,  __m512d *z_vec,
    __m512d *vx_vec, __m512d *vy_vec, __m512d *vz_vec,
    __m512d m_vec, double dt);

// Stiefel function for Newton's method, returning Gs1, Gs2, and Gs3
// Automatically inlined since it is in the header
static inline void mm_stiefel_Gs13_avx512(__m512d * Gs1, __m512d * Gs2, __m512d * Gs3, __m512d beta, __m512d X){
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
}

// Stiefel function for Halley's method, returning Gs0, Gs1, Gs2, and Gs3
// Automatically inlined since it is in the header
static inline void mm_stiefel_Gs03_avx512(__m512d * Gs0, __m512d * Gs1, __m512d * Gs2, __m512d * Gs3, __m512d beta, __m512d X){
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
}

#endif // WHFAST512_KERNEL_H