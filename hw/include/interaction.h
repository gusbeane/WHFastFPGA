#ifndef INTERACTION_H
#define INTERACTION_H

#include "whfast.h"

struct bodies_t interaction_step(struct bodies_t ss, real_t M0, real_t dt);

#endif // INTERACTION_H