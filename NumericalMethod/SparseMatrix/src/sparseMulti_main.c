/*****
 *
 *  Main program
 *  
 *  Randomly makes a sparse matrix and multiplies it with a vector
 *  
 *  Author:         Eric Jalbert
 *  Date Started:   June 20th, 2013
 *  Date Modified:  July 4th, 2013
 *  
 *****/ 

#include "sparseMulti.h"

int main()
{  
  double ** matrix = genSparseMatrix();
  Sparse * new = convertMatrix(matrix);
  int * vector = genVector();
  int * result = NULL;
  double * methodResult = malloc(sizeof(*methodResult) * SIZE);
  int i;

  printf("sparseForm:\n");
  printf("row_ptr = ");
  for(i = 0; i <= SIZE; i++)
  {
    printf("%d ", new->row_ptr[i]);
  }
  printf("\n");
  printf("val  |  col\n");
  for(i = 0; i < new->count; i++)
  {
    printf(" %.0f  |  %d\n",new->value[i], new->col[i]);
  }

  printMatrix(matrix);
  
  printf("\n      x\n\n");
  
  printVector(vector);
  
  printf("\n = \n\n");
  
  result = sparseTimesVector(new, vector);
  
  printVector(result);

  printf("\n---------------\nFrom jacobi method\n");
  
  methodResult = conjugateGradient(matrix, result);
  printDoubleVector(methodResult);
  
  freeMatrix(matrix);
  free(vector);
  free(result);
  free(methodResult);
  freeSparse(new);
  
  
  return 0;

}




