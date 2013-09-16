/*******************
 * bigGNUScriptMaker.c
 *
 * Use this program to create GNUplot Scripts when there are too 
 * many lines to write it by hand.
 * 
 * Author: Eric Jalbert
 * Date Created: July 15th, 2013
 * Date Modified: July 15th, 2013
 *
 */
 
#include <stdio.h>
#include <string.h>

/* Change the scripts filename here */
#define NAME "AlphaFinder2"

int main(){
  FILE * output = NULL;
  char fileName[100] = {};
  int i;
  
  sprintf(fileName, "script%s", NAME);
  output = fopen(fileName, "w");
  
  /* small things that could probably be written by hand go here */
  fprintf(output, "set terminal postscript eps color enhanced\nset yl \"c_{/Symbol e}\"\nset xl \"c_{crit}\"\nset xr [0:]\nset yr [0:]\nset output \"cCompare4.eps\"\nunset key\n");
  
  /* Big things that need to be looped go here */
  for(i = 0; i < 200; i++)
  {
    fprintf(output, "f%d(x) = m%d * x + b%d\n", i, i, i); 
    fprintf(output, "fit f%d(x) 'GAMMA2%.3d.dat' using 6:1 via m%d, b%d\n", i, i, i, i);
  }

  fprintf(output, "set print 'alphaValsGamma2.dat'\n");
  fprintf(output, "print ");
  for(i = 0; i < 200; i++)
  {
    fprintf(output, "m%d, ", i);
  }
  fclose(output);

  return 0;
}
