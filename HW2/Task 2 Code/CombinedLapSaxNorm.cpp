#include "CombinedLapSaxNorm.h"
#include <algorithm>
float CombinedLapSaxNorm(float (&x)[XDIM][YDIM][ZDIM], float (&z)[XDIM][YDIM][ZDIM], const float (&f)[XDIM][YDIM][ZDIM],
    float (&r)[XDIM][YDIM][ZDIM],
    const float scale)
{
    float result = 0.;
    #pragma omp parallel for reduction(max:result)
    for (int i = 1; i < XDIM-1; i++){
        for (int j = 1; j < YDIM-1; j++){
            for (int k = 1; k < ZDIM-1; k++){
                z[i][j][k] =
                    -6 * x[i][j][k]
                    + x[i+1][j][k]
                    + x[i-1][j][k]
                    + x[i][j+1][k]
                    + x[i][j-1][k]
                    + x[i][j][k+1]
                    + x[i][j][k-1];
                r[i][j][k] = z[i][j][k] * scale + f[i][j][k];
                result = std::max(result, std::abs(r[i][j][k]));
            }
        }
        
    }

    return result;
    
}