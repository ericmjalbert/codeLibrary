/*
 * -------------------------------------------------------------
 *  solver.c
 *  A program that solves an ODE.
 * -------------------------------------------------------------
 *  Author              :   Eric Jalbert
 *  Date Created        :   March 5, 2013
 *  Date Last Modified  :   August 15, 2013
 * -------------------------------------------------------------
 */
 
#include <stdio.h>      
#include <stdlib.h>     
#include <math.h>      

/* Header files with a description of contents used */
#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */

/*
 *-----------------------
 * Parameter Definitions 
 *-----------------------
 */
 
/* Boolean Definition */
#define TRUE 1
#define FALSE 0

/* Problem Constants */
#define NEQ   2                /* number of equations  */
#define RTOL  RCONST(1.0e-14)  /* scalar relative tolerance            */
#define ATOL1 RCONST(1.0e-20)  /* vector absolute tolerance components */
#define ATOL2 RCONST(1.0e-2)
#define T0    RCONST(0.0)      /* initial time           */
#define T1    RCONST(0.00001)   /* first output time      */
#define EMULT RCONST(1.0e-1)   /* multiplier factor for epsilon */ 
#define UPPER 0.1
#define LOWER 0
#define UBOUND 0.001

/* Problem Parameters */
#define ALPHA 4
#define BETA  4
#define GAMMA 4
#define MU    6
#define DELTA RCONST(1.0e-8)
#define KAPPA RCONST(0.95)

/* Looping Parameters */
#define ERROREND    RCONST(1.0e-11) /* accuracy between c values */
#define ERRORSOL    RCONST(1.0e-16)  /* accuracy between upper and lower bound *?

/* Global Variables */
realtype c;                                 /* wavespeed */
realtype epsilonRegularize = 0.99999999999;     /* regularization of D(y1) */

/*
 *-----------------------------------
 * Functions used by solver functions
 *-----------------------------------
 */
 
/*
 * D(y1), the density function
 */
realtype Diffusion(realtype y1)
{
    return DELTA * pow(y1, ALPHA) * pow(1-y1, -BETA) + epsilonRegularize * y1 * DELTA;

}

/*
 * D'(y1), the derivative of D(y1) with respect to y1
 */
realtype DiffusionPrime(realtype y1)
{
    realtype DPrime;
 
    DPrime = DELTA * pow(y1, ALPHA-1) * (ALPHA * (1-y1) + BETA * y1);
    DPrime *= pow(1-y1, -BETA-1);
 
    return DPrime + epsilonRegularize * DELTA;

}

 
/*
 * f(y1), the growth rate function
 */
realtype Growth(realtype y1)
{
    return MU * y1 * (1 - pow((y1/KAPPA), GAMMA));

}


/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/*
 * f routine. Compute function f(t,y).
 */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)

{
    realtype y1, y2, yd1, yd2;

    y1 = NV_Ith_S(y,0);
    y2 = NV_Ith_S(y,1);  

    /* This is the one where I used tau = -t to get my backwards integration */
    NV_Ith_S(ydot,0) = -Diffusion(y1) * y2;
    NV_Ith_S(ydot,1) = DiffusionPrime(y1) * y2 * y2 + c * y2 + Growth(y1);

    yd1 = NV_Ith_S(y,0);
    yd2 = NV_Ith_S(y,1);

  return(0);
}


/*
 *--------------
 * Main Program
 *--------------
 */

int main()
{
    /* Declare Variables */  
    N_Vector y;
    N_Vector abstol;
    realtype reltol;
    void * cvode_mem;
    int flag;
    realtype tout;  /* First output time */
    realtype t;
    
    /* Variables for numerical method */
    realtype waveSpeedError = 10;
    int cFound = FALSE;
    realtype cValues[100];
    realtype epsilonValues[100];
    realtype trajValues[10000];
    int i = 0;
    int j = 0;
    realtype upperBound;
    realtype lowerBound;
    int loopCount = 0;
    realtype timeChange = 1;
    realtype prevV = 100000000;
    realtype prevU = 100000000;
    
    /* Variables for writing to files */
    int goodData = FALSE;
    FILE * converg;
    FILE * phasePlane;
    FILE * travelin;
    char filename[30];
    char variLine[100];
    
    y = NULL;
    abstol = NULL;
    cvode_mem = NULL;

    /* Start loop */
    while(waveSpeedError >= ERROREND)
    {
        /* initialize problem */
        upperBound = UPPER;
        lowerBound = LOWER;
        c = (upperBound + lowerBound)/2;
        goodData = FALSE;
        
        while(cFound == FALSE)
        {
            
            /* Setting vector to initial values */
            y = N_VNew_Serial(NEQ);
            NV_Ith_S(y, 0) = (0.00001)*epsilonRegularize;
            NV_Ith_S(y, 1) = -c / (DELTA * epsilonRegularize);
            
            /* Setting absolute tolerance vector */
            abstol = N_VNew_Serial(NEQ);
            NV_Ith_S(abstol,0) = ATOL1;
            NV_Ith_S(abstol,1) = ATOL2;
            
            /* Setting relative tolerance */
            reltol = RTOL;
            
            /* Call CVodeCreate to make memory for solver and select */ 
            /* differentiation formula and use of a newton interation*/
            cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
            
            /* Initialize the memory and give the function f, */
            /* initial time T0, and the dependent variable y  */
            flag = CVodeInit(cvode_mem, f, T0, y);
            
            /* Call CVodeSVtolerances to specify the scalar relative */
            /* tolerance and vector absolute tolerances              */
            flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
            
            /* Call CVDense to specify the CVDENSE dense linear solver */
            flag = CVDense(cvode_mem, NEQ);
            
            /* Set the maximum number of steps for the solver */
            flag = CVodeSetMaxNumSteps(cvode_mem, 20000);
            
            /* Set initial time and amount of outputs to default */
            tout = T1;
            t = 0;
            loopCount = 0;
            
            /* Open file to write phase plane trajectory */
            sprintf(filename, "traj%d.dat", i);
            phasePlane = fopen(filename, "w");
            
            /*
            sprintf(filename, "travWave%02d.dat", i);
            travelin = fopen(filename, "a");
            */
            
            /* Run CVODE */
            while(1)
            {

                /* Move the vectors forward */
                flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
              
            
                /* Write data points only if they're far enough apart */
/*                if(fabs(NV_Ith_S(y,0) - prevU) >= UBOUND)*/
/*                {*/
/*                    prevU = NV_Ith_S(y,0);*/
/*                    prevV = NV_Ith_S(y,1);*/
/*                    timeChange *= .99;*/
/*                }*/
/*                else*/
/*                {*/
/*                  timeChange *= 1.01;*/
/*                }*/
                
                fprintf(phasePlane, "%d %.16f %.16f %.16f\n", loopCount, tout, NV_Ith_S(y,0), NV_Ith_S(y,1));
/*                    printf("c = %.16f | %.16f %.16f | tout = %.8f\n", c, NV_Ith_S(y,0), NV_Ith_S(y,1), tout);*/


                /* Undershot the equillibrium */
                if(NV_Ith_S(y,0) >= KAPPA)
                {
                    upperBound = c;
                    c = (upperBound + lowerBound)/2;
                    printf("c = %.16f | %.16f %.16f\n", c, NV_Ith_S(y,0), NV_Ith_S(y,1));
                    break;
                }
                
                /* Overshot the equillibrium */
                if(NV_Ith_S(y,1) >= 0)
                {
                    lowerBound = c;
                    c = (upperBound + lowerBound)/2;
                    printf("c = %.16f | %.16f %.16f\n", c, NV_Ith_S(y,0), NV_Ith_S(y,1));
                    break;
                }
                
                if (flag == CV_SUCCESS)
                {
/*                   tout += timeChange;*/
                                       
                    if(NV_Ith_S(y,0) < 0.001)
                      tout += 10;
                    else if(NV_Ith_S(y,0) < 0.01)
                      tout += 5;
                    else if(NV_Ith_S(y,0) < 0.5)
                      tout += 1;
                    else if(NV_Ith_S(y,0) < 0.8)
                      tout += 0.5;
                    else
                      tout += 0.1;
                }        
                 
                loopCount++;
            }

            /* Close the phasePlane file stream */
                  fclose(phasePlane);
                
            /*******************/
            /* Destroy problem */
            /*******************/
             
            /* Free y and abstol vectors */
            N_VDestroy_Serial(y);
            N_VDestroy_Serial(abstol);
            /* Free integrator memory */
            CVodeFree(&cvode_mem);
            /* Stop if the two bounds become sufficiently close */
                if(fabs(upperBound - lowerBound) < ERRORSOL)
                {

                    cFound = TRUE;
                    cValues[i] = c;
                    epsilonValues[i] = epsilonRegularize;
                    i++;
                                        printf("c = %.16f\n", c );
                    goodData = TRUE;

                    break;
                }
            
        }
        
        /* Now lower regularizing Epsilon */
        epsilonRegularize *= EMULT;
        
        /* Calculate waveSpeedError */    
        if(i > 1)
        {
            waveSpeedError = fabs(cValues[i-1] - cValues[i-2]);
        }
        
        /* Reset variables */
        cFound = FALSE;
        t = 0;
        
    }
    
    /* display results */
    j = i;
    while(j > 0)
    {
        printf("CValues[%d] = %.16f, Epsilon = %le\n", j-1, cValues[j-1], epsilonValues[j-1]);
        j--;
    }
    
    return 1;

}









