/*----------------------------------------------------------------------
  mconstants.c

  Created by William J. De Meo
    on 9/27/97

  This program's purpose is to find various machine parameters.
  (All numbers have been type casted in an attempt to prevent the compiler 
  from forcing higher precision.)
----------------------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utilities.h"

main()
{
  int i;
  short Short, bigS;
  long Long, bigL;
  unsigned short UShort, bigUS;
  unsigned long ULong, bigUL;
  float Float, bigF;
  double D, bigD;
  float delta, eps, negeps;
  double ddelta, deps, dnegeps;


  /* First let's check out the size of each data type on this machine */
  printf("\n      short: %d bytes.",sizeof(short));
  printf("\n        int: %d bytes.",sizeof(int));
  printf("\n       long: %d bytes.",sizeof(long));
  printf("\n      float: %d bytes.",sizeof(float));
  printf("\n     double: %d bytes.",sizeof(double));
  printf("\nlong double: %d bytes.",sizeof(long double));
printf("\n---------------------------------------------------------------------\n");

  /* Find Largest Short */
  Short = max_short();

  /* Find Largest Long */
  Long = max_long();

  /* Find Largest Unsigned Short */
  UShort = max_unsigned_short();
                                  
  /* Find Largest Unsigned Long */
  ULong = max_unsigned_long();;

  /* Find Largest Single Precision Float */
  Float = max_single_float();

  /* Find Largest Double Precision Float */
  D = max_double_float();

  /*** MACHINE EPSILON ***/
  /* Find Smallest Single Precision Float */
  /* Two methods produce different results */
  /* The results are related by:  eps = 2 * negeps  */

  printf("\n\nProof that computed epsilons are not zero:");
  printf("\n----------------------------------------------------------------------\n");
  printf("\nSingle Precision");

  printf("\nFirst method: test that addition to 1 is greater than 1\n");
  for( delta=1.0F ; (1.0F + delta) > 1.0F ;  ) {
    eps=delta;
    delta /= 1.1F;		
  }
  if( !((1.0F + eps) > 1.0F) )
    fprintf(stderr,"\nERROR: 1 + eps = %g is not greater than 1\n",
	    (double)(1.0F + eps));
  else{
    for (i = 0; i < sizeof(eps); i++)
      printf("%02x", (unsigned) ((unsigned char *)&eps)[i]);
    fprintf(stdout,"\nsingle eps: %.50f", (double)eps);
    fprintf(stdout,"\n  1 + eps = %.50f\n", (1.0 + 10000*eps));
  }

  printf("\nSecond method: test that subtraction from 1 is less than 1");
  for(delta= 1.0F; (1.0F - delta) < 1.0F;  ){
    negeps=delta;
    delta /= 1.1F;
  }
  if( !((1.0F - negeps) < 1.0F) )
    fprintf(stderr,"\nERROR: 1 - negeps = %g is not less than 1\n",
	    (double)(1.0F - negeps));
  else{
    fprintf(stdout,"\nsingle negeps: %.50f", (double)negeps);
    fprintf(stdout,"\n  1 - negeps = %.50f\n", (1.0 - 10000*negeps));
  }

  printf("\n\nDouble precision");

  printf("\nFirst method: test that addition to 1 is greater than 1\n");
  for(ddelta = 1.0; (1.0 + ddelta) > 1.0; ) {
    deps=ddelta;
    ddelta /= 1.1;		
  }
  if( !((1.0 + deps) > 1.0) )
    fprintf(stderr,"\nERROR: 1 + deps = %g is not greater than 1\n",1.0+deps);
  else{
    for (i = 0; i < sizeof(deps); i++)
      printf("%02x", (unsigned) ((unsigned char *)&deps)[i]);
    fprintf(stdout,"\ndouble eps: %.50f", deps);
    fprintf(stdout,"\n 1 + deps = %.50f",1.0+10000*deps);
  }

  printf("\nSecond method: test that subtraction from 1 is less than 1");
  ddelta = 1.0;
  while( (1.0 - ddelta) < 1.0 ) {
    dnegeps=ddelta;
    ddelta /= 1.1;		
  }
  if( !((1.0 - dnegeps) < 1.0) )
    fprintf(stderr,"\nERROR: 1 - deps = %g is not less than 1\n",1.0-dnegeps);
  else
    fprintf(stdout,"\ndouble negeps: %.50f", dnegeps);
    fprintf(stdout,"\n 1 - dnegeps = %.50f\n", 1.0 - 10000*dnegeps);

  printf("\n----------------------------------------------------------------------\n");
  printf("\n         largest short and long:  %d, %d",bigS, bigL);
  printf("\nlargest unsinged short and long:  %u, %u",bigUS, bigUL);
  printf("\n       largest float and double:  %g, %g",bigF, bigD);
  printf("\n        float mach eps, neg-eps:  %g, %g",eps, negeps);
/*	 *(unsigned *)&eps, *(unsigned *)&negeps); */
  printf("\n       double mach eps, neg-eps:  %g, %g",deps,dnegeps);
  printf("\n-------------------------------------------------------------------\n");
  printf("\nValues of Variables After Overflow");
  printf("\n----------------------------------------------------------------------\n");
  printf("                     short and long:  %d, %d",Short,Long);
  printf("\n unsinged short and unsigned long:  %u, %u",UShort, ULong);
  printf("\nsingle and double precision float:  %g, %g",Float, D);
  printf("\n----------------------------------------------------------------------\n");
  printf("\nValues of Variables After Underflow");
  printf("\n----------------------------------------------------------------------\n");
  printf("  single precision machine eps:  %g",(float)delta);
  printf("\ndouble precision machine eps:  %g",ddelta);
  printf("\n----------------------------------------------------------------------\n");

}










