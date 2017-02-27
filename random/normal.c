/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   normal.c

     Subroutines for constructing normal(0,1) random numbers 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include "prototypes.h"
#define __VERBOSE 1
#define __DEBUG 1

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function:   void normalSequence(const long *m, double *Z)
   INPUTS
      m   the number of N(0,1) random numbers desired.
      Z   on entry:  an arbitrary vector (for storing m doubles).
          on exit:  a vector of m normal(0,1) random numbers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void normalSequence(const long *M, double *Z){
  double *u, *g, *gtemp;
  long i, j, n, I= (long)0;
  long X, m=*M;

  double avg, var;  

  //X = lmalloc((long)1);
  n = 3*m; /* first generate three times as many uniforms */

  // SEED:
  X = (long) time( NULL );
  //srand ( (unsigned)time ( NULL ) );   // (if using built in uni rand num gen)

  u = dmalloc(n);
  for(i=0;i<n;i++) u[i] = unif(&X);
  // If using built in uni rand num gen:
  //for(i=0;i<n;i++) u[i] = (double) (rand() / (RAND_MAX + 1.0));

  // Check basic stats of the uniform(0,1) sequence
  moment(&n, u, &avg, &var);
  assert(avg > 0.0);  assert(avg < 1.0); assert(var > 0.001); 
  assert(var < .2); // Var(U) = 1/12 = 0.083333...

#ifdef __VERBOSE
  printf("\n   normal.c:  UAverage = %lf, UVariance = %lf", avg, var);
#endif
  g = dmalloc(m+1);  //In case m is odd, get one extra normal rv.
  gtemp = dmalloc(2);

#ifdef __DEBUG
  for (i=0; i<m+1; i++) g[i] = 11111111;  // helps quickly identify errors
  for (i=0; i<2; i++) gtemp[i] = 11111111;
#endif

  i=0; I=0;
  while(i<m){

    normal(&n, u, &I, gtemp);  // get two more normal(0,1) random numbers
    
    if (I>-1){
      g[i] = gtemp[0]; 
      g[i+1] = gtemp[1];
      i = i+2;
      I = I+2;
    }
    else { /* didn't get enough normals -- need new uniforms (VERY unlikely) */
      i=0; I=0;
      printf("\n   normal.c:  >>>>>>>>  Generating more unif(0,1) variables...\n\n");
      for(j=0;j<n;j++) u[j] = unif(&X);
      //for(j=0;i<n;j++) u[j] = (double) (rand() / (RAND_MAX + 1.0));
    }
  }

  moment(&m, g, &avg, &var);
  assert(avg<20.0); assert(avg>-20.0); assert(var<5.0); assert(var>0.0001);
#ifdef __VERBOSE
  printf("\n   normal.c:  NAverage = %lf, NVariance = %lf", avg, var);
  printf("\n   normal.c:  The normal(0,1) random sample is:");
  printf("\n   normal.c:  ");
  j=0;
  for(i=0;i<m;i++){
    j++;
    printf("%lf ",g[i]);
    if(j%10 == 0)  printf("\n   normal.c:  ");
  }
  printf("\n");
#endif

  for(i=0;i<m;i++) Z[i] = g[i];

  free(u); free(g); free(gtemp);
}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function normal(const long *n, const double *u, long *I, double *x)
  Subroutine for constructing normal(0,1) random numbers 

     INPUTS:
          n  the length of u
          u  a vector of uniform random variables
	  i  the index of the first element of u used to generate the normal rv
          x  on entry:  an arbitrary vector of length at least 2
             on exit:  first two elements are two normal random variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

void normal(const long *n, const double *u, long *I, double *x)
{
  double s;
  assert(*I>-1);  assert(*I<(*n+1));
  s = (2*u[*I] - 1)*(2*u[*I] - 1) + (2*u[*I+1] - 1)*(2*u[*I+1] - 1);
  while(s >= 1){
    *I = *I+2;
    if(*I > *n-2) {
      *I=-1;
      return;
    }
    s = (2*u[*I] - 1)*(2*u[*I] - 1) + (2*u[*I+1] - 1)*(2*u[*I+1] - 1);
  } 

  x[0] = (2*u[*I] - 1)*sqrt(-2*log(s)/s); 
  x[1] = (2*u[*I+1] - 1)*sqrt(-2*log(s)/s);

}

