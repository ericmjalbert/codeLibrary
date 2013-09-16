/********************************************************************
 *  findMax.c 
 *  A program that finds the maximum value of a function in a user-inputted
 *      interval. Uses Bisection method to approximate the value up to a 
 *      certain error bound. Doesn't work for functions that have more then
 *      one locally maximum point in the interval.
 *  
 *  Author              :   Eric Jalbert
 *  Date Created        :   February 12, 2013
 *  Date Last Modified  :   February 13, 2013
 */
 
#include <stdio.h>      /*Used for user input and showing output*/
#include <stdlib.h>     /*Used for the atof command for user bound input*/
#include <math.h>       /*Used for the math commands in Function function*/


float AbsoluteVal(float x)
{
    if(x < 0)
    {
        return -x;
    }
    
    return x;
}


float Function(float x)
{
    float value;
    
    /*Parameters*/
    float alpha = 4;
    float beta  = 4;
    float delta = pow(10, -8);
    float kappa = 0.95;
    
    /*Functions*/
    float fOfx = x * (1 - pow((x/kappa), alpha));
    float dPrimeOfx = delta * pow(x, alpha-1) / pow(1-x, beta+1) * (alpha - alpha*x + beta*x);
    
    value = 4 * dPrimeOfx * fOfx;
    
    /*test function
    value = -pow(x,2);
    */
    
    return value;
}

float BisectionMethod(float epsilon, float relativeError, float lowerBound, float upperBound)
{
    float midPoint;
    float lowerMidPoint;
    float upperMidPoint;    
    
    midPoint = (lowerBound + upperBound)/2;

    if(relativeError > epsilon)
    {
        lowerMidPoint = (lowerBound + midPoint)/2;
        upperMidPoint = (midPoint + upperBound)/2;
        
        /*If lower region has a higher value*/
        if(Function(lowerMidPoint) > Function(upperMidPoint))
        {
            relativeError = AbsoluteVal((Function(lowerMidPoint) 
                - Function(midPoint)) / Function(midPoint));
        
            midPoint = 
                BisectionMethod(epsilon, relativeError, lowerBound, midPoint);
        }
        
        /*If upper region has a higher value*/
        if(Function(upperMidPoint) > Function(lowerMidPoint))
        {
            relativeError = AbsoluteVal((Function(upperMidPoint) 
                - Function(midPoint)) / Function(midPoint));        
                        
            midPoint = 
                BisectionMethod(epsilon, relativeError, midPoint, upperBound);
        }
    }
    
    return midPoint;

}


int main()
{
    float epsilon           = 0.00000001;
    float relativeError     = 1;
    
    float max;

    float lowerBound;
    float upperBound;
    char line[100];
    
    /*Start*/
    printf("\n\n\n\n\n-----------------------------------------------------\n");
    
    /*Get interval from user*/
    printf("Enter a lower bound: ");
    fgets(line, sizeof(line), stdin);
    lowerBound = atof(line);
    
    printf("Enter an upper bound: ");
    fgets(line, sizeof(line), stdin);
    upperBound = atof(line);
    
    /*Bisection Method now*/
    max = BisectionMethod(epsilon, relativeError, lowerBound, upperBound);

    printf("The max in [%f, %f] is at %f and has Function value of %f. \n",
         lowerBound, upperBound, max, sqrt(Function(max)));

    return 1;
}




