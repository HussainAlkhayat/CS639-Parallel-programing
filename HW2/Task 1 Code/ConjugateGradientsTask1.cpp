#include "Laplacian.h"
#include "Parameters.h"
#include "PointwiseOps.h"
#include "Reductions.h"
#include "Utilities.h"
#include "Timer.h"


#include <iostream>
extern Timer timerLaplacian;

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
    Timer saxTimerTotal16b;

    Timer normTimerTotal;
    Timer innerTimerTotal;
    Timer innerTimerTotal13;

    Timer copyTimerTotal;


    // Algorithm : Line 2
    lapTimerTotal.Start(); ComputeLaplacian(x, z); lapTimerTotal.Pause();
    lapTimerTotal.Print("ComputeLaplacian on line 2: ");

    saxTimerTotal.Start(); Saxpy(z, f, r, -1); saxTimerTotal.Pause();
    saxTimerTotal.Print("Saxpy on line 2: ");

    normTimerTotal.Start(); float nu = Norm(r); normTimerTotal.Pause();
    normTimerTotal.Print("Norm on line 2: ");
    // Algorithm : Line 3
    if (nu < nuMax) return;
        
    // Algorithm : Line 4
    copyTimerTotal.Start(); Copy(r, p); copyTimerTotal.Pause();
    copyTimerTotal.Print("Copy on line 4: ");

    innerTimerTotal.Start(); float rho=InnerProduct(p, r); innerTimerTotal.Pause();
    innerTimerTotal.Print("InnerProduct on line 4: ");
    // Beginning of loop from Line 5
    innerTimerTotal.Reset();
    lapTimerTotal.Reset(); 
    copyTimerTotal.Reset();
    normTimerTotal.Reset();
    saxTimerTotal.Reset();
    for(int k=0;;k++)
    {
        std::cout << "Residual norm (nu) after " << k << " iterations = " << nu << std::endl;

        // Algorithm : Line 6
        lapTimerTotal.Restart(); ComputeLaplacian(p, z); lapTimerTotal.Pause();

        innerTimerTotal.Restart(); float sigma=InnerProduct(p, z); innerTimerTotal.Pause();
        // Algorithm : Line 7
        float alpha=rho/sigma;

        // Algorithm : Line 8
        saxTimerTotal.Restart(); Saxpy(z, r, r, -alpha); saxTimerTotal.Pause();

        normTimerTotal.Restart(); nu=Norm(r); normTimerTotal.Pause();
        // Algorithm : Lines 9-12
        if (nu < nuMax || k == kMax) {
            saxTimerTotal12.Restart(); Saxpy(p, x, x, alpha); saxTimerTotal12.Pause();
            std::cout << "Conjugate Gradients terminated after " << k << " iterations; residual norm (nu) = " << nu << std::endl;
            //if (writeIterations) WriteAsImage("x", x, k, 0, 127);
            lapTimerTotal.Print("ComputeLaplacian on line 6: ");
            innerTimerTotal.Print("InnerProduct on line 6: ");
            saxTimerTotal.Print("Saxpy on line 8: ");
            normTimerTotal.Print("Norm on line 8: ");
            saxTimerTotal12.Print("Saxpy on line 9-12: ");
            copyTimerTotal.Print("Copy on line 13: ");
            innerTimerTotal13.Print("InnerProduct on line 13: ");
            saxTimerTotal16a.Print("1st Saxpy on line 16: ");
            saxTimerTotal16b.Print("2nd Saxpy on line 16: ");
            return;
        }
            
        // Algorithm : Line 13
        copyTimerTotal.Restart(); Copy(r, z); copyTimerTotal.Pause();
        innerTimerTotal13.Restart(); float rho_new = InnerProduct(z, r); innerTimerTotal13.Pause();
        // Algorithm : Line 14
        float beta = rho_new/rho;

        // Algorithm : Line 15
        rho=rho_new;

        // Algorithm : Line 16
        saxTimerTotal16a.Restart(); Saxpy(p, x, x, alpha); saxTimerTotal16a.Pause();
        saxTimerTotal16b.Restart(); Saxpy(p, r, p, beta); saxTimerTotal16b.Pause();

        //if (writeIterations) WriteAsImage("x", x, k, 0, 127);
    }
    



}
