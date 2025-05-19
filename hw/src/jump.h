#ifndef JUMP_H
#define JUMP_H

#include "whfast.h"

#ifdef __cplusplus
extern "C" {
#endif

    struct bodies_t jump_step(struct bodies_t ss, real_t M0, real_t dt);

#ifdef __cplusplus
}
#endif

#endif // JUMP_H