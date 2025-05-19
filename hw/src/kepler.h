#ifndef KEPLER_H
#define KEPLER_H

#include "whfast.h"

#ifdef __cplusplus
extern "C" {
#endif

    struct bodies_t kepler_step(struct bodies_t ss, double M0, double dt);

    void stiefel_Gs03(double *Gs0, double *Gs1, double *Gs2, double *Gs3,
                      double beta, double X);

    void stiefel_Gs13(double *Gs1, double *Gs2, double *Gs3, double beta, double X);

    real_t halley_step(real_t X, real_t beta, real_t r0, real_t eta0, real_t zeta0,
                     real_t dt);

    real_t newton_step(real_t X, real_t beta, real_t r0, real_t eta0, real_t zeta0, real_t dt);

#ifdef __cplusplus
}
#endif

#endif // KEPLER_H