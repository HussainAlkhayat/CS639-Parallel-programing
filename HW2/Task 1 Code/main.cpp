#include "ConjugateGradientsTask1.h"
#include "Timer.h"
#include "Utilities.h"
#include <omp.h>
Timer timerLaplacian;

int main(int argc, char *argv[])
{
    using array_t = float (&) [XDIM][YDIM][ZDIM];

    float *xRaw = new float [XDIM*YDIM*ZDIM];
    float *fRaw = new float [XDIM*YDIM*ZDIM];
    float *pRaw = new float [XDIM*YDIM*ZDIM];
    float *rRaw = new float [XDIM*YDIM*ZDIM];
    float *zRaw = new float [XDIM*YDIM*ZDIM];
    
    array_t x = reinterpret_cast<array_t>(*xRaw);
    array_t f = reinterpret_cast<array_t>(*fRaw);
    array_t p = reinterpret_cast<array_t>(*pRaw);
    array_t r = reinterpret_cast<array_t>(*rRaw);
    array_t z = reinterpret_cast<array_t>(*zRaw);
    
    // Initialization
    {
        Timer timer;
        timer.Start();
        InitializeProblem(x, f);
        timer.Stop("Initialization : ");
    }

    // Call Conjugate Gradients algorithm
    omp_set_num_threads(1);
    timerLaplacian.Start();
    ConjugateGradients(x, f, p, r, z);
    timerLaplacian.Stop("Total Gradiant Time : ");
    
    return 0;
}
