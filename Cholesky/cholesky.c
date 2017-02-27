/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cholesky.c

  Created on 11/29/97 by William J. De Meo

  Purpose: Cholesky decomposition of an n-by-n spd matrix 
           
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include "prototypes.h"
#include <math.h>

/* Subroutine cholesky:

   Arguments: 
              N  dimension of A

              A 
                 on entry: the N by N matrix to be decomposed
                 on exit: upper triangle is still A
                          lower sub-triangle is the sub-trangle
                          of the Cholesky factor L

              diag 
                 on entry: an arbitrary vector of length N
                 on exit: the diagonal of the Cholesky factor L

*/

void cholesky(long N, double *A, double *diag)
{
     long i,j,k;
     for(j=0;j<N;j++)
          diag[j] = A[N*j+j];
     for(j=0;j<N;j++)
     {
          for(k=0;k<j;k++)
               diag[j] -= A[N*k+j]*A[N*k+j];
          diag[j] = sqrt(diag[j]);
          for(i=j+1;i<N;i++)
          {
               for(k=0;k<j;k++)
                    A[N*j+i] -= A[N*k+i]*A[N*k+j];
               A[N*j+i]/=diag[j];
          }
     }
}

      
           
      
           
           
      
      
