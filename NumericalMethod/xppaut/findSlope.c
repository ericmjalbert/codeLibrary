/*********************
 *  
 *  findSlope.c
 *    Just a quick program that will read my data file of points and 
 *    find the slope at each one.
 *
 *  Author: Eric Jalbert
 *  Date Created: August 1, 2013
 *  Date Modified: August 1, 2013
 */
 
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main()
{
  FILE * file = fopen("curve.dat", "r");
  FILE * output = fopen("slope.dat", "w");
  char line[100] = "";
  char * token = NULL;
  double slope=0;
  int i;
  
  while(fgets(line, sizeof(line), file) != NULL)
  {
    token = strtok(line, " ");
    token = strtok(NULL, " ");
    fprintf(output, "%f %f\n", (i-1)*.05, (atof(token)-slope)/.05);
        slope = atof(token);
    i++;
  } 
  
  
  return 0;
}





