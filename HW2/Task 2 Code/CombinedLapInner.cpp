#include "CombinedLapInner.h"

float CombinedLapInner(float (&p)[XDIM][YDIM][ZDIM], float (&z)[XDIM][YDIM][ZDIM])
{
    double result = 0.;

    #pragma omp parallel for reduction(+:result)
    for (int i = 1; i < XDIM-1; i++){
        for (int j = 1; j < YDIM-1; j++){
            for (int k = 1; k < ZDIM-1; k++){
                z[i][j][k] =
                    -6 * p[i][j][k]
                    + p[i+1][j][k]
                    + p[i-1][j][k]
                    + p[i][j+1][k]
                    + p[i][j-1][k]
                    + p[i][j][k+1]
                    + p[i][j][k-1];
                result += (double) p[i][j][k] * (double) z[i][j][k];
            }
        }
        
    }
    return (float) result;

}