#ifndef INTERACTION_H
#define INTERACTION_H

#include "whfast.h"

#ifdef __cplusplus
extern "C" {
#endif

    struct bodies_t interaction_step(struct bodies_t ss, real_t M0, real_t dt);

#ifdef __cplusplus
}
#endif

#endif // INTERACTION_H