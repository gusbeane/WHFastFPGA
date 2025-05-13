#include <immintrin.h>
#include "util.h"

void whfast_kernel(__m512d *x_vec, __m512d *y_vec, __m512d *z_vec,
                      __m512d *vx_vec, __m512d *vy_vec, __m512d *vz_vec,
                      __m512d m_vec, Body *com, double dt, long Nint);

void whfast_kepler_step(__m512d *x_vec, __m512d *y_vec, __m512d *z_vec,
                           __m512d *vx_vec, __m512d *vy_vec, __m512d *vz_vec,
                           __m512d m_vec, double dt);

void whfast_com_step(Body *com, double dt);

void whfast_drift_step(__m512d *x_vec, __m512d *y_vec, __m512d *z_vec,
                          __m512d *vx_vec, __m512d *vy_vec, __m512d *vz_vec,
                          __m512d m_vec, Body *com, double dt);

void whfast_jump_step(__m512d *x_vec, __m512d *y_vec, __m512d *z_vec,
                         __m512d *vx_vec, __m512d *vy_vec, __m512d *vz_vec,
                         __m512d m_vec, double dt);

void whfast_interaction_step(__m512d *x_vec,  __m512d *y_vec,  __m512d *z_vec,
    __m512d *vx_vec, __m512d *vy_vec, __m512d *vz_vec,
    __m512d m_vec, double dt);
