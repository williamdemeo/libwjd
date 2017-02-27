/*----------------------------------------------------------------------
  mconstants.c

  Created by William J. De Meo
    on 9/27/97

  This program's purpose is to find various machine parameters.
  (All numbers have been type casted in an attempt to prevent the compiler 
  from forcing higher precision.)
----------------------------------------------------------------------*/
#include <stdlib.h>

main()
{

  short S, bigS;
  long L, bigL;
  unsigned short US, bigUS;
  unsigned long UL, bigUL;
  float F, bigF;
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
  S = (short) 1;
  do{
    bigS = S;
    S *= (short) 2;
  }while(bigS == (S / (short) 2)); 
  /* If S*2 too big for machine, then  (S*2)/2 != S
   * and loop exits.  S stores the spurious (oversized) value.  */

  /* Find Largest Long */
  L = (long) 1;
  do{
    bigL = L;
    L *= (long) 2;
  }while(bigL == (L / (long) 2)); 


  /* Find Largest Unsigned Short */
  US = (unsigned short) 1;
  do{
    bigUS = US;
    US *= (unsigned short) 2;
  }while(bigUS == (US / (unsigned short) 2)); 
                                  
  /* Find Largest Unsigned Long */
  UL = (unsigned long) 1;
  do{
    bigUL = UL;
    UL *= (unsigned long) 2;
  }while(bigUL == (UL / (unsigned long) 2)); 

  /* Find Largest Single Precision Float */
  F = (float) 1;
  do{
    bigF = F;
    F *= (float) 2;
  }while(bigF == (F / (float) 2)); 

  /* Find Largest Double Precision Float */
  D = (double) 1;
  do{
    bigD = D;
    D *= (double) 2;
  }while(bigD == (D / (double) 2)); 


  /*** MACHINE EPSILON ***/
  /* Find Smallest Single Precision Float */
  /* Two methods produce different results */
  /* The results are related by:  eps = 2 * negeps  */

  printf("\n\nProof that computed epsilons are not zero:");
  printf("\n----------------------------------------------------------------------\n");
  printf("\nSingle precision (30 decimal display):");
  printf("\nFirst method: test that addition to 1 is greater than 1");
  for( delta=(float)1.0 ; (float)1.0 + delta > (float) 1.0 ;  ) {
    eps=delta;
    delta /= (float) 1.1;		
  }
  if(!((float)1.0 + eps) > (float)1.0)
  printf("\n   1 + eps = %.30f",(float)(1.0 + eps));
  /* Second method: test that subtraction from 1 is less than 1 */

  for(delta= (float) 1.0;(float)((float) 1.0 - delta) < (float) 1.0;){
    negeps=delta;
    delta /= (float) 1.1;
  }

  /* Find Smallest Double Precision Float (agian by two methods) */
  ddelta = (double) 1.0;
  while( (double)1.0 - ddelta < (double)1 ) {
    dnegeps=ddelta;
    ddelta /= (double) 1.1;		
  }

  for(ddelta = (double) 1.0; (double) 1.0 + ddelta > (double) 1.0; ) {
    deps=ddelta;
    ddelta /= (double) 1.1;		
  }

  printf("\n   1 + eps = %.30f",(float)(1.0 + 2*eps));
  printf("\n1 - negeps = %.30f\n", (float)(1.0 - 2*negeps));
  printf("\n\nDouble precision (25 decimal display):");
  printf("\n   1 + eps = %.30f",(double)1.0 + 2*deps);
  printf("\n1 - negeps = %.30f\n", (double)1.0 - 2*dnegeps);
  printf("\n----------------------------------------------------------------------\n");
  printf("\n         largest short and long:  %d, %d",bigS, bigL);
  printf("\nlargest unsinged short and long:  %u, %u",bigUS, bigUL);
  printf("\n       largest float and double:  %g, %g",bigF, bigD);
  printf("\n        float mach eps, neg-eps:  %.30f, %.30f",eps, negeps);
/*	 *(unsigned *)&eps, *(unsigned *)&negeps); */
  printf("\n       double mach eps, neg-eps:  %g, %g",deps,dnegeps);
  printf("\n-------------------------------------------------------------------\n");
  printf("\nValues of Variables After Overflow");
  printf("\n----------------------------------------------------------------------\n");
  printf("                     short and long:  %d, %d",S,L);
  printf("\n unsinged short and unsigned long:  %u, %u",US, UL);
  printf("\nsingle and double precision float:  %g, %g",F, D);
  printf("\n----------------------------------------------------------------------\n");
  printf("\nValues of Variables After Underflow");
  printf("\n----------------------------------------------------------------------\n");
  printf("  single precision machine eps:  %g",(float)delta);
  printf("\ndouble precision machine eps:  %g",ddelta);
  printf("\n----------------------------------------------------------------------\n");

}










