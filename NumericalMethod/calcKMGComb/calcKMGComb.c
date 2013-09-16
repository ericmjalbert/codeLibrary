/*
 * -------------------------------------------------------------
 *  calcKMGComb.c
 *  A program calculates the c_critical that was discussed in the
 *  smooth travelling wave analaysis. It is the expression max
 *  [4*d'(U)*f(U)].
 *  This program specifically solves this for multi parameter values, 
 *  namely K, Mu, and Gamma.  
 * -------------------------------------------------------------
 *  Author              :   Eric Jalbert
 *  Date Created        :   July 10, 2013
 *  Date Last Modified  :   July 10, 2013
 * -------------------------------------------------------------
 */
 
#include <stdio.h>      
#include <stdlib.h>     
#include <string.h>
#include <math.h>      

/*
 *-----------------------
 * Parameter Definitions 
 *-----------------------
 */
/* Method Parameters */
#define ERROR 1.0e-8
 
/* Problem Parameters */
#define ALPHA 4
#define BETA  4
#define DELTA 1.0e-8

/* Global Variables */
double c;                                 /* wavespeed */
double epsilonRegularize = 0.99999999999;     /* regularization of D(y1) */
double MU = 1;
double KAPPA = 0.9;
int GAMMA = 1;

/*
 *-----------------------------------
 * Functions used by solver functions
 *-----------------------------------
 */
 
/**********
 * D(y1), the density function
 */
double Diffusion(double y1)
{
    return DELTA * pow(y1, ALPHA) * pow(1-y1, -BETA) + epsilonRegularize * y1 * DELTA;

}

/**********
 * D'(y1), the derivative of D(y1) with respect to y1
 */
double DiffusionPrime(double y1)
{
    double DPrime;
 
    DPrime = DELTA * pow(y1, ALPHA-1) * (ALPHA * (1-y1) + BETA * y1);
    DPrime *= pow(1-y1, -BETA-1);
 
    return DPrime + epsilonRegularize * DELTA;
}

/**********
 * f(y1), the growth rate function
 */
double Growth(double y1)
{
    return MU * y1 * (1 - pow((y1/KAPPA), GAMMA));

}

double Function(double U)
{
  return 4 * Growth(U) * DiffusionPrime(U);
}


/*
 *--------------
 * Main Program
 *--------------
 */

int main()
{
  double upperBound, lowerBound;
  double mid, midPoint;
  char fileName[100] = "";
  FILE * grid = NULL;
  
for(GAMMA = 4; GAMMA <= 5; GAMMA += 2)
{
sprintf(fileName, "gridValues%d.dat", GAMMA);
      grid = fopen(fileName, "w");
  for(KAPPA = 0.94; KAPPA <= 0.999999; KAPPA += 0.01)
  {
    for(MU = 1; MU <= 10; MU++)
    {
      upperBound = 1;
      lowerBound = 0;
      while(fabs(upperBound - lowerBound) > ERROR)
      {
        midPoint = (upperBound +lowerBound)/2.0;
        mid = Function(midPoint);
        if(mid >= Function((midPoint + upperBound)/2.0))
          upperBound = midPoint;
        else
          lowerBound = midPoint; 
      }
      
      fprintf(grid,"%f %f %f\n",KAPPA, MU,  sqrt(Function((upperBound + lowerBound)/2.0)));
      
    }
    fprintf(grid, "\n");
  }
}

  return 0;
}

















