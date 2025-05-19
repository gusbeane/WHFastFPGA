#pragma once
#include <immintrin.h>
#include "util.h"

struct Constants {
    __m512d invfactorial512[35];
    __m512d gr_prefac;
    __m512d gr_prefac2;
    __m512d half;
    __m512d one;
    __m512d two;
    __m512d five;
    __m512d sixteen;
    __m512d twenty;
    double M0;
    __m512d _M;
    __m512i so2;
    __m512i so1;
};

extern const Constants* kConsts;

// Call this once with the desired M0/m_vec to build kConsts.
void initialize_constants(double M0, __m512d m_vec);

#ifdef PRINT_UTILITY
extern double beta_min, beta_max;
extern double X_min, X_max;
extern double halley_r0_min, halley_r0_max;
extern double halley_eta0_min, halley_eta0_max;
extern double halley_zeta0_min, halley_zeta0_max;

extern double newton_r0_min, newton_r0_max;
extern double newton_eta0_min, newton_eta0_max;
extern double newton_zeta0_min, newton_zeta0_max;
#endif // PRINT_UTILITY