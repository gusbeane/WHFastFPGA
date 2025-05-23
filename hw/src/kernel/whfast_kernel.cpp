#include "whfast.h"
#include "jump.h"
#include "interaction.h"
#include "kepler.h"

struct bodies_t set_gr_prefac(struct bodies_t ss, real_t M0)
{
#pragma HLS inline
    real_t c_val = R(10065.32);
    ss.gr_prefac = R(6.) * M0 * M0 / (c_val * c_val);
    for (int i = 0; i < N_PLANETS; i++)
    {
        ss.gr_prefac2[i] = ss.m_vec[i] / M0;
    }

    return ss;
}

struct bodies_t whfast_kernel(struct bodies_t ss, real_t M0, real_t dt, long Nint)
{
    // Call the main integration routine
    // This is where the actual integration happens
    // This function assumes the input is _desynchronized_

    real_t dt_half1 = dt / 2.0;
    real_t dt_half2 = dt - dt_half1;

    ss = set_gr_prefac(ss, M0);

    #pragma HLS array_partition variable = ss.x_vec complete
    #pragma HLS array_partition variable = ss.y_vec complete
    #pragma HLS array_partition variable = ss.z_vec complete
    #pragma HLS array_partition variable = ss.vx_vec complete
    #pragma HLS array_partition variable = ss.vy_vec complete
    #pragma HLS array_partition variable = ss.vz_vec complete
    #pragma HLS array_partition variable = ss.m_vec complete
#pragma HLS array_partition variable = ss.gr_prefac2 complete

// #pragma HLS ALLOCATION function instances = jump_step limit = 2
// #pragma HLS ALLOCATION function instances = interaction_step limit = 1
// #pragma HLS ALLOCATION function instances = kepler_step limit = 1

whfast_kernel_loop:
    for (int i = 0; i < Nint; i++)
    {
#pragma HLS unroll off
#pragma HLS loop_tripcount min = 1 max = 3
        // Perform the jump step (first half)
        ss = jump_step(ss, M0, dt_half1);

        // Perform the interaction step
        ss = interaction_step(ss, M0, dt);

        // Perform the jump step (second half)
        ss = jump_step(ss, M0, dt_half2);

        // Perform the kepler step
        ss = kepler_step(ss, M0, dt);
    }

    return ss;
}