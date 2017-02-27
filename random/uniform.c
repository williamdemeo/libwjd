/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  unif.c

  Created by William J. De Meo
  Last modified 11/30/97

  Purpose: generate uniform(0,1) random variables
           using a multiplicative congruential generator based 
           on Schrage's approximate factorization algorithm
           (see document /accounts/grad/chip/classes/s243/reports/hw2
           for details)

  Argument: 

          X (pointer to long)
             on entry:  a seed
             on exit:  A * (X mod Q) - R*[X/Q]       if >=0
                       A * (X mod Q) - R*[X/Q] + M   otherwise
            (the constants are defined below)
         
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <math.h>
#include <stdio.h>

#define A 16807 
#define M 2147483647 
#define Q 127773
#define R 12614

double unif(long *X)
{
     long k;
     double unif;

     k = (*X)/Q;                   /* [x/q]                    */
     *X = A * (*X - k*Q) - R*k;    /* x - k*q = x mod q        */
     if(*X < 0) *X += M;
     unif = (*X)*((double)1.0/M);  /* make it uniform on (0,1) */
     return(unif);
}






