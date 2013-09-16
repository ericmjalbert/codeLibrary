/*
 * -------------------------------------------------------------
 *  psuedoXppaut.c
 *  A program that plots the flow lines of a phase portrait.
 * -------------------------------------------------------------
 *  Author              :   Eric Jalbert
 *  Date Created        :   May 4, 2013
 *  Date Last Modified  :   May 9, 2013
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
#define ATOL1 RCONST(1.0e-24)  /* vector absolute tolerance components */
#define ATOL2 RCONST(1.0e-2)
#define T0    RCONST(0.0)      /* initial time           */
#define T1    RCONST(0.001)      /* first output time      */
#define TMULT RCONST(1.01)

/* Problem Parameters */
#define ALPHA 4
#define BETA  4
#define GAMMA 4
#define MU    6
#define DELTA RCONST(1.0e-8)
#define KAPPA RCONST(0.95)

/* Global Variables */
#define NUMTRAJ 4
#define UBOUND 10e-4
#define VBOUND 10e2
#define EQUIL 10e-3
#define INITU 1/(float)NUMTRAJ 
#define INITV -c/DELTA/epsilonRegularize/NUMTRAJ
realtype c = 0.004436019063;                /* wavespeed */
realtype epsilonRegularize = 0.9999999;     /* regularization of D(y1) */

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
    return DELTA * pow(y1, ALPHA) / pow(1-y1, BETA) + DELTA * y1 * epsilonRegularize;

}

/*
 * D'(y1), the derivative of D(y1) with respect to y1
 */
realtype DiffusionPrime(realtype y1)
{
    realtype DPrime;
 
    DPrime = DELTA * pow(y1, ALPHA-1) * (ALPHA * (1-y1) + BETA * y1);
    DPrime *= pow(1-y1, -BETA-1);
 
    return DPrime + DELTA * epsilonRegularize;

}

 
/*
 * f(y1), the growth rate function
 */
realtype Growth(realtype y1)
{
    return MU * y1 * (1 - pow((y1/KAPPA), GAMMA));

}


/*
 * The first Vnullcline
 */
double Vcline1(double y1)
{
    if(c*c - 4* DiffusionPrime(y1) * Growth(y1) >= 0)
        return (c - sqrt(c*c - 4* DiffusionPrime(y1) * Growth(y1)))/(-2*DiffusionPrime(y1));
    else
        return 10000;
}


/*
 * The other V-Nullcline
 */
double Vcline2(double y1)
{
    if(c*c - 4* DiffusionPrime(y1) * Growth(y1) >= 0)
        return (c + sqrt(c*c - 4* DiffusionPrime(y1) * Growth(y1)))/(-2*DiffusionPrime(y1));
    else
        return 10000;
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
    realtype y1, y2;

    y1 = NV_Ith_S(y,0);
    y2 = NV_Ith_S(y,1);  

    /* This is the one where I used tau = -t to get my backwards integration */
    NV_Ith_S(ydot,0) = -Diffusion(y1) * y2;
    NV_Ith_S(ydot,1) = DiffusionPrime(y1) * y2 * y2 + c * y2 + Growth(y1);

  return(0);
}

/*
 * forwards in time function f(t,y).
 */
static int fForward(realtype t, N_Vector y, N_Vector ydot, void *user_data)

{
    realtype y1, y2;

    y1 = NV_Ith_S(y,0);
    y2 = NV_Ith_S(y,1);  

    /* This is the one where I used tau = -t to get my backwards integration */
    NV_Ith_S(ydot,0) = Diffusion(y1) * y2;
    NV_Ith_S(ydot,1) = -DiffusionPrime(y1) * y2 * y2 - c * y2 - Growth(y1);

  return(0);
}

   
/*
 *--------------------------
 * Function to write results to a graph 
 *--------------------------
 */
void GraphPhasePlane(int eloop)
{
     int i = 0;
     int j = 0;
     FILE * phasePlane;
        
        phasePlane = popen("gnuplot", "w");
    
        fprintf(phasePlane, "set terminal pdf\n");
        fprintf(phasePlane, "set xlabel \"U\"\n");
        fprintf(phasePlane, "set ylabel \"V\"\n");
        fprintf(phasePlane, "set xrange [0:1]\n");
        fprintf(phasePlane, "set yrange[%f:%f]\n", -c/DELTA/epsilonRegularize*1.1, c/DELTA/epsilonRegularize*.5);
        fprintf(phasePlane, "unset key\n");
        fprintf(phasePlane, "set arrow from 0,0 to 1,0 nohead lt 0\n");
        fprintf(phasePlane, "set output \"refPhasePlaneGrid_epsilon10e-%d\"\n", eloop);
        fprintf(phasePlane, "set multiplot\n");
        for(j = -1; j <= NUMTRAJ; j++)
        {
            for(i = 1; i <= NUMTRAJ; i++)
            {
                fprintf(phasePlane, "plot 'regPhaseGridBack%d%d_epsilon10e-%d.dat' with lines linecolor 0\n", i, j, eloop);
                fprintf(phasePlane, "plot 'regPhaseGridForw%d%d_epsilon10e-%d.dat' with lines linecolor 0\n", i, j, eloop);
            }
        }
        fprintf(phasePlane, "plot 'regPhaseTravWave_epsilon10e-%d.dat' with lines lw 3 linecolor 1 title 'Trav. Wave'\n", eloop);
        fprintf(phasePlane, "plot 'regPhaseVNull_epsilon10e-%d.dat' with lines lw 2 linecolor 3 title 'V-Nullcline'\n", eloop);
        
        pclose(phasePlane);
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
    int i = 1;
    int j = 1;
    int loopCount = 0;
    realtype prevU = 0;
    realtype prevV = 0;
    realtype timeChange;
    
    /* Variables for writing to files */
    FILE * phasePlane;
    int eloop = 1;
    char filename[100];
    char variLine[100];
    
    /* Variables for Nullclines */
    double u = 0.0001;
    double vChange1 = 0.01;
    double vChange2 = 0.01;
    double vcline1 = 0;
    double vcline2 = 0;
    double cline1Prev = 1000000;
    double cline2Prev = 1000000;   

    
for(epsilonRegularize = 0.1; epsilonRegularize >= 0.00001; epsilonRegularize*=0.1)
{

printf("calculating backwards integrations\n");

for(i = 1; i <= NUMTRAJ; i++)
{
    printf("ayj = %d, i = %d\n", j, i);
    
    for(j = -1; j <= NUMTRAJ; j++)
    {
        printf("j = %d, i = %d\n", j, i);
        y = NULL;
        abstol = NULL;
        cvode_mem = NULL;
                    
        /* Set initial time and amount of outputs to default */
        tout = T1;
        timeChange = tout/(float)i;
        t = 0;
        loopCount = 0;
            
        /* Setting vector to initial values */
        y = N_VNew_Serial(NEQ);
        NV_Ith_S(y, 0) = INITU*i - INITU / 2;
        NV_Ith_S(y, 1) = INITV*j - INITU / 2;
            
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
        flag = CVodeSetMaxNumSteps(cvode_mem, 2000);
            
        /* Open file to write phase plane trajectory */
        sprintf(filename, "regPhaseGridBack%d%d_epsilon10e-%d.dat", i,j,eloop);
        phasePlane = fopen(filename, "w");
        
        /* Write the first data point */    
        fprintf(phasePlane, "%.16f %.16f\n", NV_Ith_S(y,0), NV_Ith_S(y,1));
            
        /* Run CVODE forwards */
        while(1)
       {
            /* Move the vectors forward */
            flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
        
            /* Write data points only if they're far enough apart */
            if(fabs(NV_Ith_S(y,0) - prevU) >= UBOUND || fabs(NV_Ith_S(y,1) - prevV) >= VBOUND/epsilonRegularize)
            {
            printf("%.16f %.16f %f\n", NV_Ith_S(y,0), NV_Ith_S(y,1), t);
                fprintf(phasePlane, "%.16f %.16f %f\n", NV_Ith_S(y,0), NV_Ith_S(y,1), t);
                prevU = NV_Ith_S(y,0);
                prevV = NV_Ith_S(y,1);
                timeChange *= 0.95;
            }
            else timeChange *= 1.01;
       
            /* Trajectory outside of region */
            if(NV_Ith_S(y,0) >= 1 || NV_Ith_S(y,1) >= c/DELTA/epsilonRegularize*.5 || NV_Ith_S(y,0) <= 0 || NV_Ith_S(y,1) <= -c/DELTA/epsilonRegularize*1.1) break;
            
            /* Trajectory approchs origin equillibrium */
            if(fabs(NV_Ith_S(y,0)-0) <= EQUIL && fabs(NV_Ith_S(y,1)-0) <= EQUIL*20) break; 
  
            /* Traj approches Kappa equillibrium */
            if(fabs(NV_Ith_S(y,0) - KAPPA) <= EQUIL && fabs(NV_Ith_S(y,1) - 0) <= EQUIL*20) break; 
  
            /* Traj. approches quasi-equillibrium */
            if(fabs(NV_Ith_S(y,0) - 0) <= EQUIL && fabs(NV_Ith_S(y,1) - c/DELTA/epsilonRegularize) <= EQUIL*200) break; 
            
            /* Traj approches U = 1 */
            if(fabs(NV_Ith_S(y,0) - 1) <= EQUIL) break;
  
            if(flag == CV_SUCCESS) tout += timeChange;
        
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
    } 
}

/*
 *-------------------------------------------
 * Now Same thing but forwards intergrations
 *-------------------------------------------
 */
 printf("calculating forward trajectories\n");
for(i = 1; i <= NUMTRAJ; i++)
{
    printf("ayj = %d, i = %d\n", j, i);
    
    for(j = -1; j <= NUMTRAJ; j++)
    {
        printf("j = %d, i = %d\n", j, i);
        y = NULL;
        abstol = NULL;
        cvode_mem = NULL;
                    
        /* Set initial time and amount of outputs to default */
        tout = T1;
        timeChange = tout/(float)i;
        t = 0;
        loopCount = 0;
            
        /* Setting vector to initial values */
        y = N_VNew_Serial(NEQ);
        NV_Ith_S(y, 0) = INITU*i - INITU / 2;
        NV_Ith_S(y, 1) = INITV*j - INITU / 2;
            
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
        flag = CVodeInit(cvode_mem, fForward, T0, y);
            
        /* Call CVodeSVtolerances to specify the scalar relative */
        /* tolerance and vector absolute tolerances              */
        flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
            
        /* Call CVDense to specify the CVDENSE dense linear solver */
        flag = CVDense(cvode_mem, NEQ);
            
        /* Set the maximum number of steps for the solver */
        flag = CVodeSetMaxNumSteps(cvode_mem, 2000);
            
        /* Open file to write phase plane trajectory */
        sprintf(filename, "regPhaseGridForw%d%d_epsilon10e-%d.dat", i,j,eloop);
        phasePlane = fopen(filename, "w");
        
        /* Write the first data point */    
        fprintf(phasePlane, "%.16f %.16f\n", NV_Ith_S(y,0), NV_Ith_S(y,1));
        
        /* Run CVODE forwards */
        while(1)
       {
            /* Move the vectors forward */
            flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
        
            /* Write data points only if they're far enough apart */
            if(fabs(NV_Ith_S(y,0) - prevU) >= UBOUND || fabs(NV_Ith_S(y,1) - prevV) >= VBOUND/epsilonRegularize)
            {
            printf("%.16f %.16f %f\n", NV_Ith_S(y,0), NV_Ith_S(y,1), t);
                fprintf(phasePlane, "%.16f %.16f %f\n", NV_Ith_S(y,0), NV_Ith_S(y,1), t);
                prevU = NV_Ith_S(y,0);
                prevV = NV_Ith_S(y,1);
                timeChange *= 0.95;
            }
            else timeChange *= 1.01;
       
            /* Trajectory outside of region */
            if(NV_Ith_S(y,0) >= 1 || NV_Ith_S(y,1) >= c/DELTA/epsilonRegularize*.5 || NV_Ith_S(y,0) <= 0 || NV_Ith_S(y,1) <= -c/DELTA/epsilonRegularize*1.1) break;
            
            /* Trajectory approchs origin equillibrium */
            if(fabs(NV_Ith_S(y,0)-0) <= EQUIL && fabs(NV_Ith_S(y,1)-0) <= EQUIL*20) break; 
  
            /* Traj approches Kappa equillibrium */
            if(fabs(NV_Ith_S(y,0) - KAPPA) <= EQUIL && fabs(NV_Ith_S(y,1) - 0) <= EQUIL*20) break; 
  
            /* Traj. approches quasi-equillibrium */
            if(fabs(NV_Ith_S(y,0) - 0) <= EQUIL && fabs(NV_Ith_S(y,1) - c/DELTA/epsilonRegularize) <= EQUIL*200) break; 
            
            /* Traj approches U = 1 */
            if(fabs(NV_Ith_S(y,0) - 1) <= EQUIL) break;
  
            if(flag == CV_SUCCESS) tout += timeChange;
        
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
    } 
}

/*
 *----------------------------------------------------------------------
 * Now Graph out the travelling wave trajectory from quasi-equillibrium
 *----------------------------------------------------------------------
 */
 
 printf("Calculating Traveling wave\n");
 y = NULL;
        abstol = NULL;
        cvode_mem = NULL;
                    
        /* Set initial time and amount of outputs to default */
        tout = T1;
        timeChange = tout/(float)i;
        t = 0;
        loopCount = 0;
            
        /* Setting vector to initial values */
        y = N_VNew_Serial(NEQ);
        /*
        NV_Ith_S(y,0) = KAPPA - 0.000000001;
        NV_Ith_S(y,1) = -0.000000001;
        */
        NV_Ith_S(y, 0) = 0.0000000000000001;
        NV_Ith_S(y, 1) = -c/DELTA/epsilonRegularize;
         
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
        flag = CVodeSetMaxNumSteps(cvode_mem, 2000);
            
        /* Open file to write phase plane trajectory */
        sprintf(filename, "regPhaseTravWave_epsilon10e-%d.dat", eloop);
        phasePlane = fopen(filename, "w");
        
        /* Write the first data point */    
        fprintf(phasePlane, "%.16f %.16f\n", NV_Ith_S(y,0), NV_Ith_S(y,1));
            
        /* Run CVODE forwards */
        while(1)
       {
            /* Move the vectors forward */
            flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
        
            /* Write data points only if they're far enough apart */
            if(fabs(NV_Ith_S(y,0) - prevU) >= UBOUND || fabs(NV_Ith_S(y,1) - prevV) >= VBOUND/epsilonRegularize)
            {
                fprintf(phasePlane, "%.16f %.16f %f\n", NV_Ith_S(y,0), NV_Ith_S(y,1), t);
                printf("%.16f %.16f %f\n", NV_Ith_S(y,0), NV_Ith_S(y,1), t);
                prevU = NV_Ith_S(y,0);
                prevV = NV_Ith_S(y,1);
                timeChange *= 0.95;
            }
            else timeChange *= 1.01;
            
            /* Traj approches Kappa equillibrium */
            if(fabs(NV_Ith_S(y,0) - KAPPA) <= EQUIL && fabs(NV_Ith_S(y,1) - 0) <= EQUIL*20) break; 
            
            
                        /* Trajectory outside of region */
            if(NV_Ith_S(y,0) >= 1 || NV_Ith_S(y,1) >= c/DELTA/epsilonRegularize*.5|| NV_Ith_S(y,0) <= 0 || NV_Ith_S(y,1) <= -c/DELTA/epsilonRegularize*1.1) break;
  
            if(flag == CV_SUCCESS) tout += timeChange;
        
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
        
     
     
/*
 *-------------------------------------------------
 * Now get the regularized V-nullcline data points
 *-------------------------------------------------
 */
    printf("calculating V-nullcline\n");
    u = 0;
    cline1Prev = 10000000;
    cline2Prev = -10000000;
    sprintf(filename, "regPhaseVNull_epsilon10e-%d.dat", eloop);
    phasePlane = fopen(filename,"w");
        
    while(1)
    {
        if(u >= 1) break;
        vcline1 = Vcline1(u);
        if(vcline1 == 10000) break;
        
        if(fabs(vcline1 - cline1Prev) >= vChange1 && vcline1 != 10000)
        {
           fprintf(phasePlane, "%f %f\n", u, vcline1);
           cline1Prev = vcline1;
           vChange1 *= 1.1;
        }
        else vChange1 *= 0.9;
            
        u += 0.00001;
    }
         
    while(1)
    {
        if(u <= 0) break;
        vcline2 = Vcline2(u);
            
        if(fabs(vcline2 - cline2Prev) >= vChange2 && vcline2 != 10000)
        {
            fprintf(phasePlane, "%f %f\n", u, vcline2);      
            cline2Prev = vcline2;      
            vChange2 *= 1.1;
        }
        else vChange2 *= 0.7;

        u -= 0.00001;
    }
         
    
        fclose(phasePlane);

GraphPhasePlane(eloop);

eloop++;

}

 
    return 1;

}









