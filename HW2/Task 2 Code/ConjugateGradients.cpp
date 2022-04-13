#include "Laplacian.h"
#include "Parameters.h"
#include "PointwiseOps.h"
#include "Reductions.h"
#include "Utilities.h"
#include "Timer.h"
#include "CombinedSaxpy.h"
#include "CombinedLapSaxNorm.h"
#include "CombinedCopyInner.h"
#include "CombinedLapInner.h"
#include "CombinedSaxNorm.h"

#include <iostream>

void ConjugateGradients(
    float (&x)[XDIM][YDIM][ZDIM],
    const float (&f)[XDIM][YDIM][ZDIM],
    float (&p)[XDIM][YDIM][ZDIM],
    float (&r)[XDIM][YDIM][ZDIM],
    float (&z)[XDIM][YDIM][ZDIM],
    const bool writeIterations)
{
    Timer lapTimerTotal;
    Timer saxTimerTotal;
    Timer saxTimerTotal12;
    Timer saxTimerTotal16a;

    Timer normTimerTotal;
    Timer innerTimerTotal;
    Timer innerTimerTotal13;

    Timer copyTimerTotal;


    // Algorithm : Line 2
    normTimerTotal.Start(); float nu = CombinedLapSaxNorm(x, z, f, r, -1); normTimerTotal.Pause();
    normTimerTotal.Print("CombinedLapSaxNorm on line 2: ");
    // Algorithm : Line 3
    if (nu < nuMax) return;
        
    // Algorithm : Line 4
    innerTimerTotal.Start(); float rho=CombinedCopyInner(p, r); innerTimerTotal.Pause();
    innerTimerTotal.Print("CombinedCopyInner on line 4: ");
    // Beginning of loop from Line 5
    innerTimerTotal.Reset();
    lapTimerTotal.Reset(); 
    copyTimerTotal.Reset();
    normTimerTotal.Reset();
    saxTimerTotal.Reset();
    for(int k=0;;k++)
    {
        // Algorithm : Line 6

        innerTimerTotal.Restart(); float sigma=CombinedLapInner(p, z); innerTimerTotal.Pause();
        // Algorithm : Line 7
        float alpha=rho/sigma;

        // Algorithm : Line 8
        normTimerTotal.Restart(); nu=CombinedSaxNorm(z, r, -alpha); normTimerTotal.Pause();
        // Algorithm : Lines 9-12
        if (nu < nuMax || k == kMax) {
            saxTimerTotal12.Restart(); Saxpy(p, x, x, alpha); saxTimerTotal12.Pause();
            std::cout << "Conjugate Gradients terminated after " << k << " iterations; residual norm (nu) = " << nu << std::endl;
            //if (writeIterations) WriteAsImage("x", x, k, 0, 127);
            innerTimerTotal.Print("CombinedLapInner on line 6: ");
            normTimerTotal.Print("CombinedSaxNorm on line 8: ");
            saxTimerTotal12.Print("Saxpy on line 9-12: ");
            innerTimerTotal13.Print("CombinedCopyInner on line 13: ");
            saxTimerTotal16a.Print(" CombinedSaxpy on line 16: ");
            
            return;
        }
            
        // Algorithm : Line 13
        innerTimerTotal13.Restart(); float rho_new = CombinedCopyInner(z, r); innerTimerTotal13.Pause();
        // Algorithm : Line 14
        float beta = rho_new/rho;

        // Algorithm : Line 15
        rho=rho_new;

        // Algorithm : Line 16
        saxTimerTotal16a.Restart(); CombinedSaxpy(p, x, r, alpha, beta); saxTimerTotal16a.Pause();
        
        //if (writeIterations) WriteAsImage("x", x, k, 0, 127);
    }
    



}
