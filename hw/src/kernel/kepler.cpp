#include "kepler.h"

// Scalar Stiefel function for Halley's method, returning Gs0, Gs1, Gs2, and Gs3
void stiefel_Gs03(real_t* Gs0, real_t* Gs1, real_t* Gs2, real_t* Gs3, real_t beta, real_t X) {
#pragma HLS interface mode = ap_none register port = Gs0
#pragma HLS interface mode = ap_none register port = Gs1
#pragma HLS interface mode = ap_none register port = Gs2
#pragma HLS interface mode = ap_none register port = Gs3

// #pragma HLS pipeline II=1
#pragma HLS inline

    real_t X2 = X * X;
    real_t z = X2 * beta;
    // Use a truncated Taylor series for the Stumpff functions
    // nmax = 11 as in the vectorized version
    const int nmax = 11;
    real_t invfact[12] = {
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
    *Gs0 = -z * (*Gs2) + R(1.0);
    *Gs3 = (*Gs3) * X;
    *Gs1 = -z * (*Gs3) + X;
    *Gs3 = (*Gs3) * X2;
    *Gs2 = (*Gs2) * X2;
}

// Scalar Stiefel function for Newton's method, returning Gs1, Gs2, and Gs3
void stiefel_Gs13(real_t* Gs1, real_t* Gs2, real_t* Gs3, real_t beta, real_t X) {
#pragma HLS inline
// Change #1
// II 193 BRAM 2 DSP 19 FF 3642 LUT 3471
// II 124 BRAM 0 DSP 19 FF 2870 LUT 2217
#pragma HLS interface mode = ap_none register port = Gs1
#pragma HLS interface mode = ap_none register port = Gs2
#pragma HLS interface mode = ap_none register port = Gs3

// Change #2 (disabled due to inlining)
// II 124 BRAM 0 DSP 19  FF 2870  LUT 2217
// II 1   BRAM 0 DSP 227 FF 19815 LUT 15407
// #pragma HLS pipeline II=1

    real_t X2 = X * X;
    real_t z = X2 * beta;
    // Use a truncated Taylor series for the Stumpff functions
    // nmax = 19 as in the vectorized version
    const int nmax = 19;
    real_t invfact[20] = {
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

// Top‐level inlining + single‐cycle II
// static inline void stiefel_Gs13(
// void stiefel_Gs13(
//     real_t *Gs1,
//     real_t *Gs2,
//     real_t *Gs3,
//     real_t beta,
//     real_t X)
// {
// // #pragma HLS inline            // inline into caller
// // #pragma HLS pipeline II = 2 // Initiate every cycle
// #pragma HLS interface mode = ap_none register port = Gs1
// #pragma HLS interface mode = ap_none register port = Gs2
// #pragma HLS interface mode = ap_none register port = Gs3

//     // Step 1: precompute constants
//     real_t X2, z;
// // #pragma HLS bind_op variable=X2 op=mul impl=dsp latency=2
//     X2 = X * X;
//     z = X2 * beta;


//     static const real_t G3coeff[9] = {
//         0.16666666666666666,    //  1/3!
//         0.008333333333333333,   //  1/5!
//         0.0001984126984126984,  //  1/7!
//         2.7557319223985893e-06, //  1/9!
//         2.505210838544172e-08,  //  1/11!
//         1.6059043836821613e-10, //  1/13!
//         7.647163731819816e-13,  //  1/15!
//         2.8114572543455206e-15, //  1/17!
//         8.22063524662433e-18    //  1/19!
//     };

//     static const real_t G2coeff[9] = {
//         0.5,                    // 1/2!
//         0.041666666666666664,   // 1/4!
//         0.001388888888888889,   // 1/6!
//         2.48015873015873e-05,   // 1/8!
//         2.755731922398589e-07,  // 1/10!
//         2.08767569878681e-09,   // 1/12!
//         1.1470745597729725e-11, // 1/14!
//         4.779477332387385e-14,  // 1/16!
//         1.5619206968586225e-16  // 1/18!
//     };

// // Fully partition the array so each element is in its own register bank.
// // That lets you read any two coeffs simultaneously.
// #pragma HLS ARRAY_PARTITION variable = G3coeff complete dim = 1
// #pragma HLS ARRAY_PARTITION variable = G2coeff complete dim = 1

//     //  Reimplement using Estrin's method
//     real_t y = -z;
//     real_t y2 = y * y;
//     real_t y4 = y2 * y2;
//     real_t y8 = y4 * y4;

//     // Computing G3
//     real_t G3_p0 = G3coeff[0] + G3coeff[1] * y;
//     real_t G3_p1 = G3coeff[2] + G3coeff[3] * y;
//     real_t G3_p2 = G3coeff[4] + G3coeff[5] * y;
//     real_t G3_p3 = G3coeff[6] + G3coeff[7] * y;
//     real_t G3_p4 = G3coeff[8]; // the lonely top term

//     real_t G3_q0 = G3_p0 + G3_p1 * y2;
//     real_t G3_q1 = G3_p2 + G3_p3 * y2;

//     real_t G3_r = G3_q0 + G3_q1 * y4;
//     real_t G3 = G3_r + G3_p4 * y8;

//     // Computing G2
//     real_t G2_p0 = G2coeff[0] + G2coeff[1] * y;
//     real_t G2_p1 = G2coeff[2] + G2coeff[3] * y;
//     real_t G2_p2 = G2coeff[4] + G2coeff[5] * y;
//     real_t G2_p3 = G2coeff[6] + G2coeff[7] * y;
//     real_t G2_p4 = G2coeff[8]; // the lonely top term

//     real_t G2_q0 = G2_p0 + G2_p1 * y2;
//     real_t G2_q1 = G2_p2 + G2_p3 * y2;

//     real_t G2_r = G2_q0 + G2_q1 * y4;
//     real_t G2 = G2_r + G2_p4 * y8;

//     // Step 4: finish off Gs1, Gs2, Gs3
//     real_t G3_x = G3 * X; // one extra multiply
//     *Gs1 = -z * G3_x + X; // one multiply + one add
//     *Gs3 = G3_x * X2;     // one multiply
//     *Gs2 = G2 * X2;       // one multiply
// }

// Scalar Halley step function for Kepler iteration
real_t halley_step(real_t X, real_t beta, real_t r0, real_t eta0, real_t zeta0, real_t dt) {
    // #pragma HLS pipeline II=1
    #pragma HLS inline

    real_t Gs0, Gs1, Gs2, Gs3;
    stiefel_Gs03(&Gs0, &Gs1, &Gs2, &Gs3, beta, X);

    real_t f = r0 * X - dt;
    f += eta0 * Gs2;
    f += zeta0 * Gs3;

    real_t fp = r0 + eta0 * Gs1 + zeta0 * Gs2;
    real_t fpp = eta0 * Gs0 + zeta0 * Gs1;
    
    real_t denom = fp * fp * R(16.0) - R(20.0) * f * fpp;
    denom = sqrt(denom);
    denom = fp + denom;
    
    X = (X * denom - R(5.0) * f) / denom;
    
    return X;
}

// Scalar Newton step function for Kepler iteration
real_t newton_step(real_t X, real_t beta, real_t r0, real_t eta0, real_t zeta0, real_t dt)
{
// #pragma HLS pipeline II = 1
#pragma HLS inline

    real_t Gs1, Gs2, Gs3;
    stiefel_Gs13(&Gs1, &Gs2, &Gs3, beta, X);
    real_t eta0Gs1zeta0Gs2 = eta0 * Gs1 + zeta0 * Gs2;

    real_t ri = R(1.0) / (r0 + eta0Gs1zeta0Gs2);
    X = X * eta0Gs1zeta0Gs2 - eta0 * Gs2 - zeta0 * Gs3 + dt;
    X = ri * X;

    return X;
}

struct bodies_t kepler_step(struct bodies_t ss, real_t M0, real_t dt)
{
#pragma HLS inline off
#pragma HLS pipeline II = 10

    for (int i = 0; i < N_PLANETS; i++)
    {
#pragma HLS UNROLL

        real_t r2, r0, r01, v2;
        real_t beta, eta0, zeta0;
        real_t Gs1, Gs2, Gs3, eta0Gs1zeta0Gs2, ri;
        r2 = ss.x_vec[i] * ss.x_vec[i] + ss.y_vec[i] * ss.y_vec[i] + ss.z_vec[i] * ss.z_vec[i];
        r0 = sqrt(r2);
        r01 = R(1.0) / r0;

        v2 = ss.vx_vec[i] * ss.vx_vec[i] + ss.vy_vec[i] * ss.vy_vec[i] + ss.vz_vec[i] * ss.vz_vec[i];
        beta = R(2.0) * M0 * r01 - v2;

        eta0 = ss.x_vec[i] * ss.vx_vec[i] + ss.y_vec[i] * ss.vy_vec[i] + ss.z_vec[i] * ss.vz_vec[i];
        zeta0 = -beta * r0 + M0;

        // Initial guess
        real_t dtr0i = dt * r01;
        real_t X = dtr0i * eta0;
        X = -R(0.5) * X * r0 + R(1.0);
        X = dtr0i * X;

        // Iterations
        X = halley_step(X, beta, r0, eta0, zeta0, dt);
        X = halley_step(X, beta, r0, eta0, zeta0, dt);
        X = newton_step(X, beta, r0, eta0, zeta0, dt);

        stiefel_Gs13(&Gs1, &Gs2, &Gs3, beta, X);
        eta0Gs1zeta0Gs2 = eta0 * Gs1 + zeta0 * Gs2;
        ri = R(1.0) / (r0 + eta0Gs1zeta0Gs2);

        // f and g function
        real_t nf = M0 * Gs2; // negative f
        nf = nf * r01;

        real_t g = -M0 * Gs3 + dt;
        real_t nfd = M0 * Gs1; // negative fd
        nfd = nfd * r01;
        nfd = nfd * ri;

        real_t ngd = M0 * Gs2; // negative gd
        ngd = ngd * ri;

        real_t nx = -nf * ss.x_vec[i] + ss.x_vec[i];
        nx = nx + g * ss.vx_vec[i];

        real_t ny = -nf * ss.y_vec[i] + ss.y_vec[i];
        ny = ny + g * ss.vy_vec[i];

        real_t nz = -nf * ss.z_vec[i] + ss.z_vec[i];
        nz = nz + g * ss.vz_vec[i];

        ss.vx_vec[i] = -ngd * ss.vx_vec[i] + ss.vx_vec[i];
        ss.vx_vec[i] = -nfd * ss.x_vec[i] + ss.vx_vec[i];
        ss.vy_vec[i] = -ngd * ss.vy_vec[i] + ss.vy_vec[i];
        ss.vy_vec[i] = -nfd * ss.y_vec[i] + ss.vy_vec[i];
        ss.vz_vec[i] = -ngd * ss.vz_vec[i] + ss.vz_vec[i];
        ss.vz_vec[i] = -nfd * ss.z_vec[i] + ss.vz_vec[i];

        ss.x_vec[i] = nx;
        ss.y_vec[i] = ny;
        ss.z_vec[i] = nz;
    }

    return ss;
}