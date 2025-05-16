#include <cmath>
#include <ap_fixed.h>

#define N_PLANETS 8

// typedef ap_fixed<48, 10> real_t;
typedef double real_t;
#define R(x) real_t(x)

extern "C"
{
    void kepler_step(double x_vec[N_PLANETS], double y_vec[N_PLANETS], double z_vec[N_PLANETS], double vx_vec[N_PLANETS], double vy_vec[N_PLANETS], double vz_vec[N_PLANETS],double m_vec[N_PLANETS], double M0, double dt);

    void stiefel_Gs13(
        double* Gs1,
        double* Gs2,
        double* Gs3,
        double  beta,
        double  X
    );
}