#pragma once
#include <immintrin.h>
#include "util.h"

struct Constants {
    double invfactorial512[35];
    double gr_prefac[N_PLANETS];
    double gr_prefac2[N_PLANETS];
    double half;
    double one;
    double two;
    double five;
    double sixteen;
    double twenty;
    double M0;
};

extern const Constants* kConsts;

// Call this once with the desired M0/m_vec to build kConsts.
void initialize_constants(double M0, double *m_vec);

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