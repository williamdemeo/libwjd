/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  normal_test.c

  Purpose: main program for testing routine normal()
           for generating normal(0,1) random variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "prototypes.h"
#define VERBOSE 1
#define DEBUG 1

main()
{
  long n, m, i;
  double *Z, avg, var;
  m=(long)25;

  printf("\n   normal_test.c:  Getting %ld normal(0,1) random numbers...", m);
  Z = dmalloc(m); 
  normalSequence(&m, Z);
  moment(&m, Z, &avg, &var);
  printf("\n   normal_test.c:  Average = %lf, Variance = %lf\n",avg,var);
  printf("\n   normal_test.c:  The normal random numbers are:");
  printf("\n   normal_test.c:  ");
  for(i=0;i<m;i++){
    printf("%lf ",Z[i]);
    if((i+1)%10 == 0)   printf("\n   normal_test.c:  ");
  }
  printf("\n");
  free(Z);
}

