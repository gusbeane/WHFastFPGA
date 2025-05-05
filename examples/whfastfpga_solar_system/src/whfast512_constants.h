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
    __m512d _M;
    __m512i so2;
    __m512i so1;
};

extern const Constants* kConsts;

// Call this once with the desired M0/m_vec to build kConsts.
void initialize_constants(double M0, __m512d m_vec);