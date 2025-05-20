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

static Constants buildConstants(double M0, double *m_vec) {
    Constants c;
    c.half = 0.5;
    c.one = c.half + c.half;
    c.two = c.one + c.one;
    c.five = 5.;
    c.sixteen = 16.;
    c.twenty = 20.;
    c.M0 = M0;
    
    for (unsigned int i = 0; i < 35; i++){
        c.invfactorial512[i] = invfactorial[i];
    }

    // GR prefactors; example setup:
    double c_val = 10065.32;
    for (std::size_t p = 0; p < N_PLANETS; p++){
        c.gr_prefac[p] = 6.*M0*M0/(c_val*c_val);
        c.gr_prefac2[p] = m_vec[p]/M0; // or any appropriate value from m_vec
    }

    return c;
}

void initialize_constants(double M0, double *m_vec) {
    if (!kConsts_local) {
        // Build a static immutable constants once.
        static const Constants c = buildConstants(M0, m_vec);
        kConsts_local = &c;
        kConsts = kConsts_local;
    }
}