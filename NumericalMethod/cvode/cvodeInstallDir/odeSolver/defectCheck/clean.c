/* Just cleans the wavespeed data files */

#include<stdio.h>
#include<string.h>
#include<stdlib.h>


int main()
{
  FILE * file = fopen("wavespeed.dat", "r");
  FILE * output = fopen("waveClean.dat", "w");
  char line[100] = {};
  char linecopy[100] = {};
  char * token = NULL;
  char ** endptr= NULL;
  int fileNum = 0;
  char fileName[100];  
  int i;
  
  for(i = 0; i < 8; i++)
  {
    fclose(file);
    sprintf(fileName, "backTrans%d.dat", i);
    file = fopen(fileName, "r");
      
        fclose(output);
        sprintf(fileName, "waveClean%d.dat", i);
        output = fopen(fileName, "w");      
      
    while(fgets(line, sizeof(line), file) != NULL)
    {

      if(strcmp(line, "\n") == 0 || fileNum == 0){

        fileNum++;
      }
      else
      {
        strcpy(linecopy, line);
        token = strtok(line, " ");
        
        if(strcmp(token, "nan") != 0 && strcmp(token, "-nan") != 0)
          if(strtod(token, endptr) != 0)
          {
            token = strtok(NULL, " ");
            if(strcmp(token, "nan\n") != 0 && strcmp(token, "-nan\n") != 0 && strcmp(token, "inf\n") != 0 && strcmp(token, "-inf\n") != 0)
            {
              fprintf(output, "%s", linecopy);
            }
          }
      }
    }
  }
  
fclose(file);
fclose(output);

  return 1;
}





