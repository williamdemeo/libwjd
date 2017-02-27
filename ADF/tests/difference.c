/************************************************************
 * difference.c  main program for testing difference fn     *
 *                                                          *
 * Created by William DeMeo                                 *
 * on 2011.09.24                                            *
 *                                                          *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <prototypes.h>

void difference(long N, double *x);


void difference(long N, double *x) {
  double *Lx, alpha=-1.0; 
  long INC=1, nrow=N-1;
  Lx = dmalloc(nrow);
  dcopy_(&nrow, x+1, &INC, Lx, &INC);        // Lx(0:N-2) <- x(1:N-1)
  daxpy_(&nrow, &alpha, x, &INC, Lx, &INC);  // Lx(0:N-2) <- Lx(0:N-2) - x(0:N-2)
  dcopy_(&nrow, Lx, &INC, x, &INC);          // x(0:N-2) <- Lx(0:N-2) = x(1:N-1)-x(0:N-2)
  free(Lx);
}
    
  
#ifdef __TEST_DIFFERENCE__
void main()
{
  double *diffy;
  int i, j, len=0, place=1;
  long N=10, INC=1;
  //sigLevels = [0.001 (0.005:0.005:0.10) (0.125:0.025:0.20) (0.80:0.025:0.875) (0.90:0.005:0.995) 0.999];
  double y[10] = {1, 2, 3, 5, 7, 9, 10, 11, 11, 10};

  diffy = dmalloc(N);
  dcopy_(&N, y, &INC, diffy, &INC);
  printf("\ny     = ");  matprint(y, 1, N);
  printf("\ndiffy = ");  matprint(diffy, 1, N);
  difference(N, diffy);
  printf("\ndiffy = ");  matprint(diffy, 1, N-1);
  free(diffy);  free(y);
}
#endif
