#include "CombinedSaxpy.h"

// alg line 16
void CombinedSaxpy(float (&p)[XDIM][YDIM][ZDIM], float (&x)[XDIM][YDIM][ZDIM],
    float (&r)[XDIM][YDIM][ZDIM], const float alpha, const float beta)
{
    #pragma omp parallel for
    for (int i = 1; i < XDIM-1; i++){
        for (int j = 1; j < YDIM-1; j++){
            for (int k = 1; k < ZDIM-1; k++){
                x[i][j][k] = p[i][j][k] * alpha + x[i][j][k];
                p[i][j][k] = p[i][j][k] * beta + r[i][j][k];
            }
        
        }
    
    }
    
}
