#ifndef WHFAST_H
#define WHFAST_H

#include <ap_fixed.h>
#include <cmath>

#define N_PLANETS 8

// typedef ap_fixed<48, 10> real_t;
typedef double real_t;
#define R(x) real_t(x)

struct bodies_t
{
    real_t x_vec[N_PLANETS];
    real_t y_vec[N_PLANETS];
    real_t z_vec[N_PLANETS];
    real_t vx_vec[N_PLANETS];
    real_t vy_vec[N_PLANETS];
    real_t vz_vec[N_PLANETS];
    real_t m_vec[N_PLANETS];
};

struct bodies_t whfast_kernel(struct bodies_t ss, real_t M0, real_t dt, long Nint   );

#endif // WHFAST_H