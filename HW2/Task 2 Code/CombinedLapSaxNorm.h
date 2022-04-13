#include "Parameters.h"
float CombinedLapSaxNorm(float (&x)[XDIM][YDIM][ZDIM], float (&z)[XDIM][YDIM][ZDIM], const float (&f)[XDIM][YDIM][ZDIM],
    float (&r)[XDIM][YDIM][ZDIM],
    const float scale);