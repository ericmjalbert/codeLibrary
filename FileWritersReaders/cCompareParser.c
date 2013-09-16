/*******************
 * cCompareParser.c
 *
 * Parses though the cCompare.dat files and separates each set of 
 * K values (every 20 lines) into a different file.
 * 
 * Author: Eric Jalbert
 * Date Created: July 15th, 2013
 * Date Modified: July 15th, 2013
 *
 */
 
#include <stdio.h>
#include <string.h>

int main(){
  FILE * file = fopen("cCompare2.dat", "r");
  FILE * output = NULL;
  char outputName[100] = {};
  char line[250] = {};
  int lineCount, kCount;
  
  for(kCount = 0; kCount < 200; kCount++)
  {
    sprintf(outputName, "GAMMA2%.3d.dat", kCount);
    output = fopen(outputName, "w");
    for(lineCount = 0; lineCount < 20; lineCount++)
    {
      if(fgets(line, sizeof(line), file) != NULL)
      {
        fprintf(output, "%s", line);
      }
      else
        printf("Yeah this isn't working properly...\n");
    }
  }
  fclose(file);
  fclose(output);
  
  return 0;
}




