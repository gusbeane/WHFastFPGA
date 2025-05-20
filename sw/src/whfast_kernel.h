#include <immintrin.h>
#include "util.h"

void whfast_kernel(double *x_vec, double *y_vec, double *z_vec,
                      double *vx_vec, double *vy_vec, double *vz_vec,
                      double *m_vec, Body *com, double dt, long Nint);

void whfast_kepler_step(double *x_vec, double *y_vec, double *z_vec,
                           double *vx_vec, double *vy_vec, double *vz_vec, double dt);

void whfast_com_step(Body *com, double dt);

void whfast_drift_step(double *x_vec, double *y_vec, double *z_vec,
                          double *vx_vec, double *vy_vec, double *vz_vec, Body *com, double dt);

void whfast_jump_step(double *x_vec, double *y_vec, double *z_vec,
                      double *vx_vec, double *vy_vec, double *vz_vec,
                      double *m_vec, double dt);

void whfast_interaction_step(double *x_vec, double *y_vec, double *z_vec,
    double *vx_vec, double *vy_vec, double *vz_vec,
    double *m_vec, double dt);

void stiefel_Gs03(double* Gs0, double* Gs1, double* Gs2, double* Gs3, double beta, double X);

void stiefel_Gs13(double* Gs1, double* Gs2, double* Gs3, double beta, double X);

void halley_step(double* X, double beta, double r0, double eta0, double zeta0, double dt);

void newton_step(double* X, double beta, double r0, double eta0, double zeta0, double dt);
