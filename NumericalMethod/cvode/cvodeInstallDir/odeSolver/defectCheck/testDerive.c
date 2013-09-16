/*
 *  ---------------------------------
 *  testDerive.c
 *  ---------------------------------
 *  Just a program to test out the 3-point finite difference thing
 *  ---------------------------------
 */

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>

double Function(double z)
{
  return pow(z, 2); 
}

double SimpsonRule(double U0, double U1, double U2, double x0, double x2)
{
  double h = (x0 - x2)/2.0;
  
  return h/3 * (U0 + 4 * U1 + U2);
}

int main()
{
  char line[100];
  char * token = NULL;
  char ** endptr = NULL;
  FILE * output = fopen("testOutput.dat", "w");
  FILE * file = fopen("backTrans2.dat", "r");
  double z[100000];
  double U[100000];
  double uDeriv[100000];
  double uDeriv2[100000];
  double uInt[100000];
  
  int count, i;
  
  
  while(fgets(line, sizeof(line), file) != NULL)
    {
      token = strtok(line, " ");
      z[count] = strtod(token, endptr);
      count++;
    }
    
  for(i = 0; i < count; i++)
  {
    U[i] = Function(z[i]);
  }  
  
  
  /* Derive */
  
  for(i = 0; i < count; i++)
  {
    if(i == 0)
    {
      uDeriv[i] = -(z[i+1] + z[i+2] - 2*z[i])*U[i+3] / (z[i+3] - z[i+1])/(z[i+3] - z[i+2]) - (z[i+3] + z[i+2] - 2*z[i])*U[i+1] / (z[i+1] - z[i+3]) / (z[i+1] - z[i+2]) - (z[i+3] + z[i+1] - 2* z[i]) * U[i+2] / (z[i+2] - z[i+3]) / (z[i+2] - z[i+1]);
    
      uDeriv2[i] = 2 * U[i+3] / (z[i+3] - z[i+1])/(z[i+3] - z[i+2]) + 2 * U[i+1] / (z[i+1] - z[i+3]) / (z[i+1] - z[i+2]) + 2 * U[i+2] / (z[i+2] - z[i+3]) / (z[i+2] - z[i+1]);
    }
    else if (i >= count - 2)
    {
      uDeriv[i] = -(z[i-1] + z[i-2] - 2*z[i])*U[i-3] / (z[i-3] - z[i-1])/(z[i-3] - z[i-2]) - (z[i-3] + z[i-2] - 2*z[i])*U[i-1] / (z[i-1] - z[i-3]) / (z[i-1] - z[i-2]) - (z[i-3] + z[i-1] - 2* z[i]) * U[i-2] / (z[i-2] - z[i-3]) / (z[i-2] - z[i-1]);
    
      uDeriv2[i] = 2 * U[i-3] / (z[i-3] - z[i-1])/(z[i-3] - z[i-2]) + 2 * U[i-1] / (z[i-1] - z[i-3]) / (z[i-1] - z[i-2]) + 2 * U[i-2] / (z[i-2] - z[i-3]) / (z[i-2] - z[i-1]);
    }
    else
    {
      uDeriv[i] = -(z[i+1] + z[i+2] - 2*z[i])*U[i-1] / (z[i-1] - z[i+1])/(z[i-1] - z[i+2]) - (z[i-1] + z[i+2] - 2*z[i])*U[i+1] / (z[i+1] - z[i-1]) / (z[i+1] - z[i+2]) - (z[i-1] + z[i+1] - 2* z[i]) * U[i+2] / (z[i+2] - z[i-1]) / (z[i+2] - z[i+1]);
    
      uDeriv2[i] = 2 * U[i-1] / (z[i-1] - z[i+1]) / (z[i-1] - z[i+2]) + 2 * U[i+1] / (z[i+1] - z[i-1]) / (z[i+1] - z[i+2]) + 2 * U[i+2] / (z[i+2] - z[i-1]) / (z[i+2] - z[i+1]);
  
    }
  }
  
  /* Integrate */
  for(i = 0; i < count; i++)
  {
      if(i == 0)
        uInt[i] = SimpsonRule(U[i], U[i+1], U[i+2], z[i], z[i+2]);
      else if(i == count)
        uInt[i] = SimpsonRule(U[i-2], U[i-1], U[i], z[i-2], z[i]);
      else
        uInt[i] = SimpsonRule(U[i-1], U[i], U[i+1], z[i-1], z[i+1]); 
        fprintf(output, "%f %f %f %f %f\n", z[i], U[i], uDeriv[i], uDeriv2[i], uInt[i]);      
  }
  
  

  
  
  fclose(file);
  fclose(output);
  
  return 0;
}

