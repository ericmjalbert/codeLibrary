#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define KAPPA 0.95
#define ALPHA 4
#define BETA  4
#define GAMMA 4
#define MU    6
#define DELTA 1.0e-8

float D(float y1)
{
    return DELTA * pow(y1, ALPHA) / pow(1-y1, BETA) ;

}

int main()
{

    char line[2000];
    char * num;
    double temp[40000];
    int loopa[40000];
    double zetaVal[40000];
    double uVal[40000];
    double vVal[40000];
    double zVal[40000];
    int loopar[40000];
    double zetaValr[40000];
    double uValr[40000];
    double vValr[40000];
    double zValr[40000];
    int i = 0, ir = 0, j,loop;
    char fileName[100];
    FILE * file;
    
  for(loop = 7; loop < 8; loop++)
  {
    i = 0;
    ir = 0;
    /* Read in file */
    sprintf(fileName, "traj%dr.dat",loop);
    file = fopen(fileName, "r");
    while (fgets(line, sizeof(line), file) != NULL)
    {
        num = strtok(line, " ");
        loopar[ir] = atoi(num);
        num = strtok(NULL, " ");
        zetaValr[ir] = atof(num);
        num = strtok(NULL, " ");
        uValr[ir] = atof(num);
        num = strtok(NULL, " ");
        vValr[ir] = atof(num);
        ir++;
    } 
    fclose(file);
    
    sprintf(fileName, "traj%d.dat",loop);
    file = fopen(fileName, "r");
    while (fgets(line, sizeof(line), file) != NULL)
    {
        num = strtok(line, " ");
        loopa[i] = atoi(num);
        num = strtok(NULL, " ");
        zetaVal[i] = atof(num);
        num = strtok(NULL, " ");
        uVal[i] = atof(num);
        num = strtok(NULL, " ");
        vVal[i] = atof(num);
        i++;
    } 
    fclose(file);
    
    
    
    
    zValr[i] = zetaValr[i];
    for(j = ir-1; j > 0; j--)
    {
      zValr[ir-j] = zValr[ir-j+1] + 2 * (zetaValr[j] + zetaValr[j-1]) * D(uValr[j]) * D(uValr[j-1])/(D(uValr[j]) + D(uValr[j-1]));
      if(isnan(zValr[ir-j]) != 0)
        zValr[ir-j] = 0;
    }
    
    /* Reverse the order of uVal and vVal */
    for(j = 0; j < i; j++)
      temp[j] = uVal[i-j];
    for(j = 0; j < i; j++)
      uVal[j] = temp[j];
    for(j = 0; j < i; j++)
      temp[j] = vVal[i-j];
    for(j = 0; j < i; j++)
      vVal[j] = temp[j];
       
    zVal[i] = zetaVal[i];
    for(j = i-1; j > 0; j--)
    {
      zVal[i-j] = zVal[i-j+1] + 2 * (zetaVal[j] + zetaVal[j-1]) * D(uVal[j]) * D(uVal[j-1])/(D(uVal[j]) + D(uVal[j-1]));
      if(isnan(zVal[i-j]) != 0)
        zVal[i-j] = 0;
    }
    
    
    

    /* Write to file */
    sprintf(fileName, "backTrans%dReg.dat",loop);
    file = fopen(fileName, "w");
    for(j = 0; j < ir; j++)
    {
        fprintf(file, "%d %.6e %.6e %.6e %.6e\n", loopar[j], zetaValr[j], -zValr[ir-j], uValr[j], vValr[j]);
    }
    fclose(file);

  
    /* Write to file */
    sprintf(fileName, "backTrans%d.dat",loop);
    file = fopen(fileName, "w");
    for(j = 0; j < i; j++)
    {
        fprintf(file, "%d %.6e %.6e %.6e %.6e\n", loopa[j], zetaVal[j], -zVal[i-j], uVal[j], vVal[j]);
    }
    fclose(file);

  }

    return 1;

}




