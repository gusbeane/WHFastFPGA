#include <immintrin.h>
#include "util.h"

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

// void whfast512_interaction_step(__m512d *x_vec,  __m512d *y_vec,  __m512d *z_vec,
//     __m512d *vx_vec, __m512d *vy_vec, __m512d *vz_vec,
//     __m512d m_vec, double dt);
