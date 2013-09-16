/*******************
 * alphaValsParse.c
 *
 * Parses though the alphaVals.dat file and separates each value 
 *   into two colums, the first has the K values (which are know
 *   beforehand) and the second has the alpha values.
 * 
 * Author: Eric Jalbert
 * Date Created: July 13th, 2013
 * Date Modified: July 13th, 2013
 *
 */
 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int main(){
  FILE * file = fopen("alphaValsGamma2.dat", "r");
  FILE * output = fopen("alphaValsOrderGamma2.dat", "w");
  char line[20000] = {};
  char * token = NULL;
  char * delim = " ";
  int i = 0;
  
  fgets(line, sizeof(line), file);
  
  token = strtok(line, delim);
  fprintf(output, "%f %.10f\n", 0.8+i*0.001, atof(token));
  
  for(i = 1; ; i++)
  {
    token = strtok(NULL, delim);
    if(token == NULL)
      break;
    fprintf(output, "%f %.10f\n", 0.8+i*0.001, atof(token));
  }
  
  fclose(file);
  fclose(output);
  
  return 0;
}




