#include "Laplacian.h"
#include <omp.h>
#include <stdio.h>
void ComputeLaplacian(const float (&u)[XDIM][YDIM], float (&Lu)[XDIM][YDIM])
{
    
omp_set_num_threads(4); //Since the output of the max thread possible is 4. I tried 1-4 in this line 
//printf("(num of threads %d)", omp_get_max_threads(); I wrote this line ALONE without the previous line to get the max thread that are
// able to be used and the number is 4
#pragma omp parallel for
    for (int i = 1; i < XDIM-1; i++)
    for (int j = 1; j < YDIM-1; j++)
        Lu[i][j] =
            -4 * u[i][j]
            + u[i+1][j]
            + u[i-1][j]
            + u[i][j+1]
            + u[i][j-1];

}
