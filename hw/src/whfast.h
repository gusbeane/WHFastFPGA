#ifndef WHFAST_H
#define WHFAST_H

#include <ap_fixed.h>
#include <cmath>

#define N_PLANETS 8

struct bodies_t
{
    double x_vec[N_PLANETS];
    double y_vec[N_PLANETS];
    double z_vec[N_PLANETS];
    double vx_vec[N_PLANETS];
    double vy_vec[N_PLANETS];
    double vz_vec[N_PLANETS];
    double m_vec[N_PLANETS];
};

// typedef ap_fixed<48, 10> real_t;
typedef double real_t;
#define R(x) real_t(x)

#endif // WHFAST_H