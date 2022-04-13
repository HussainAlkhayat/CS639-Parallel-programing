#include "Reductions.h"
#include <mkl.h>
#include <iostream>
#include <algorithm>
using namespace std;
float Norm( float (&x)[XDIM][YDIM][ZDIM])
{
    int index = 0;
	float* tempX = **x;

	index = cblas_isamax(
		 XDIM * YDIM * ZDIM,
		 &x[0][0][0],
		 1
	 );
	float result = abs(tempX[index]);

    return result;



}

float InnerProduct(const float (&x)[XDIM][YDIM][ZDIM], const float (&y)[XDIM][YDIM][ZDIM])
{
    double result = 0.;
result = cblas_sdot (
    XDIM * YDIM * ZDIM, 
    &x[0][0][0], 
    1, 
    &y[0][0][0], 
    1);


    return (float) result;
}
