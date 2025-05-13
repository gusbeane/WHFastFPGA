#include "whfast_constants.h"
#include "whfast.h"

// Assume invfactorial is defined here or included appropriately.
static const double invfactorial[35] = {1., 1., 1./2., 1./6., 1./24., 1./120., 1./720., 1./5040.,
    1./40320., 1./362880., 1./3628800., 1./39916800., 1./479001600., 1./6227020800.,
    1./87178291200., 1./1307674368000., 1./20922789888000., 1./355687428096000.,
    1./6402373705728000., 1./121645100408832000., 1./2432902008176640000.,
    1./51090942171709440000., 1./1124000727777607680000., 1./25852016738884976640000.,
    1./620448401733239439360000., 1./15511210043330985984000000., 1./403291461126605635584000000.,
    1./10888869450418352160768000000., 1./304888344611713860501504000000.,
    1./8841761993739701954543616000000., 1./265252859812191058636308480000000.,
    1./8222838654177922817725562880000000., 1./263130836933693530167218012160000000.,
    1./8683317618811886495518194401280000000., 1./295232799039604140847618609643520000000.};

static const Constants* kConsts_local = nullptr;
const Constants* kConsts = nullptr;

static Constants buildConstants(double M0, __m512d m_vec) {
    Constants c;
    c.half = _mm512_set1_pd(0.5);
    c.one = _mm512_add_pd(c.half, c.half);
    c.two = _mm512_add_pd(c.one, c.one);
    c.five = _mm512_set1_pd(5.);
    c.sixteen = _mm512_set1_pd(16.);
    c.twenty = _mm512_set1_pd(20.);
    double M[8];
    for (int i = 0; i < 8; i++){
        M[i] = M0;
    }
    c.M0 = M0;
    c._M = _mm512_loadu_pd(M);

    c.so1 = _mm512_set_epi64(1, 2, 3, 0, 6, 7, 4, 5);
    c.so2 = _mm512_set_epi64(3, 2, 1, 0, 6, 5, 4, 7);
    for (unsigned int i = 0; i < 35; i++){
        c.invfactorial512[i] = _mm512_set1_pd(invfactorial[i]);
    }
    // GR prefactors; example setup:
    double c_val = 10065.32;
    double _gr_prefac[8] = {0};
    double _gr_prefac2[8] = {0};
    for (int p = 0; p < N_BODIES-1; p++){
        _gr_prefac[p] = 6.*M0*M0/(c_val*c_val);
        _gr_prefac2[p] = m_vec[p]/M0; // or any appropriate value from m_vec
    }
    c.gr_prefac = _mm512_loadu_pd(_gr_prefac);
    c.gr_prefac2 = _mm512_loadu_pd(_gr_prefac2);

    return c;
}

void initialize_constants(double M0, __m512d m_vec) {
    if (!kConsts_local) {
        // Build a static immutable constants once.
        static const Constants c = buildConstants(M0, m_vec);
        kConsts_local = &c;
        kConsts = kConsts_local;
    }
}