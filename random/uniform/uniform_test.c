/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  uniform_test.c  
  
  Created by William J. De Meo
  Last modified on 11/30/97

  Purpose:  Test subroutine unif() for generating uniform(0,1)
            random variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <math.h>
#include "prototypes.h"

main()
{
     unsigned long n, i;
     double *u, *ave, *var;
     long *X;

     ave = dmalloc(1);
     var = dmalloc(1);
     X = lmalloc(1);
     *X = time('\0');

     printf("\nHow many uniform random variables? ");
     scanf("%d",&n);
     u = dmalloc(n);

     printf("The unif(0,1) random variables are:\n\n");
     for(i=0;i<n;i++)
     { 
          u[i] = unif(X);
          printf("%lf  ",u[i]);
          if((i+1)%6 == 0) printf("\n");
     }
     cmoment(u, n, ave, var);
     printf("\nAverage = %lf, Variance = %lf\n",*ave,*var);
}


