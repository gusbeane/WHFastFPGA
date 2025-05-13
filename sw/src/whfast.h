#ifndef WHFAST_H
#define WHFAST_H

#include <array>
#include "util.h"

// Integrate taking time-step dt up to total time tmax; returns status
int whfast_integrate(std::array<Body, N_BODIES>& solarsystem, Body *com, double dt, long Nint);

#endif // WHFAST_H