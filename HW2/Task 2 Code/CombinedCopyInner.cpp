#include "Parameters.h"

#include <algorithm>

float CombinedCopyInner(float (&r)[XDIM][YDIM][ZDIM], float (&p)[XDIM][YDIM][ZDIM])
{
    double result = 0.;

    #pragma omp parallel for reduction(+:result)
    for (int i = 1; i < XDIM-1; i++){
        for (int j = 1; j < YDIM-1; j++){
            for (int k = 1; k < ZDIM-1; k++){
                r[i][j][k] = p[i][j][k];
                result += (double) p[i][j][k] * (double) r[i][j][k];
            }
        }
    }
    return (float) result;
}




