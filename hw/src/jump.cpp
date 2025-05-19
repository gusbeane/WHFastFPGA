#include "jump.h"

struct bodies_t jump_step(struct bodies_t ss, real_t M0, real_t dt)
{
    #pragma HLS pipeline II=1
    
    // pf = dt / M0
    double pf = dt / M0;
    double sumx = 0.0, sumy = 0.0, sumz = 0.0;
    for (int i = 0; i < N_PLANETS; i++)
    {
        #pragma HLS UNROLL factor=N_PLANETS
        sumx += ss.m_vec[i] * ss.vx_vec[i];
        sumy += ss.m_vec[i] * ss.vy_vec[i];
        sumz += ss.m_vec[i] * ss.vz_vec[i];
    }
    for (int i = 0; i < N_PLANETS; i++)
    {
        #pragma HLS UNROLL factor=N_PLANETS
        ss.x_vec[i] += sumx * pf;
        ss.y_vec[i] += sumy * pf;
        ss.z_vec[i] += sumz * pf;
    }
    
    return ss;
}