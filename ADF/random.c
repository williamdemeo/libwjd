/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  random.c

  Purpose: main program for testing routine normal()
           for generating normal(0,1) random variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "prototypes.h"

double max(double a, double b, double c);
double min(double a, double b, double c);
long I;
#define E 2.71828182845905

main()
{
     long n, m, i,j;
     double *u, *uAlt, *g, *gtemp, *ave, *var, minimum = 100, maximum = -100;
     long *X;
     FILE *ofp;
     I=(long)0;
     ave = dmalloc((long)1);
     var = dmalloc((long)1);
     gtemp = dmalloc((long)2);
     X = lmalloc((long)1);
     *X = (long) time( NULL );
     
     printf("\n\nlog(E) = %lf\n", log(E));
     printf("\nX = %ld", *X);
     //printf("\nHow many normal random variables? ");
     //scanf("%ld",&m);
     m=25;
     n = 3*m; /* generate three times as many uniforms */

     //ofp = fopen("norm.out", "a");

     srand ( (unsigned)time ( NULL ) );

     u = dmalloc(n);
     for(i=0;i<n;i++) u[i] = unif(X);
     //for(i=0;i<n;i++) u[i] = (double) (rand() / (RAND_MAX + 1.0));


     cmoment(u, n, ave, var);
     //fprintf(ofp,"\n\nUAverage = %lf, UVariance = %lf\n",*ave,*var);
     printf("\n\nUAverage = %lf, UVariance = %lf\n",*ave,*var);

     //fprintf(ofp,"\nThe normal(0,1) random variables are:\n\n");
     //printf("\nThe normal(0,1) random variables are:\n\n");
     g = dmalloc(m+1);  //In case m is odd, get one extra normal rv.
     for (i=0; i<m+1; i++) g[i] = 11111111;  // This helps identify errors faster.

     i=0;
     while(i<m){
       if(normal(u, n, gtemp)==1){
	 g[i] = gtemp[0]; 
	 g[i+1] = gtemp[1];
	 maximum = max(maximum,g[i], g[i+1]);
	 minimum = min(minimum,g[i], g[i+1]);
	 i = i+2;
	 ////fprintf(ofp,"%lf  %lf  ",g[i-1], g[i]); 
	 //if((i+1)%6 == 0) 
	 ////fprintf(ofp,"\n"); 
	 //printf("\n"); 
       }
       else { /* didn't get enough normals -- need new uniforms */
	 I=0;
	 printf("\n\n >>>>>>>>  Generating different unif(0,1) variables...\n\n");
	 for(j=0;j<n;j++) u[j] = unif(X);
	 //for(j=0;i<n;j++) u[j] = (double) (rand() / (RAND_MAX + 1.0));
       }
     }
     cmoment(g, m, ave, var);
     //fprintf(ofp,"\n\nAverage = %lf, Variance = %lf, min = %lf, max = %lf \n",*ave,*var, minimum, maximum);
     printf("\n\nAverage = %lf, Variance = %lf, min = %lf, max = %lf \n",*ave,*var, minimum, maximum);
     printf("\n\n The normal random variables are:\n");
     for(i=0;i<m;i++){
       printf("%lf ",g[i]);
       if((i+1)%10 == 0) printf("\n"); 
     }
     printf("\n");

     //fclose(ofp);
}

double max(double a, double b, double c)
{
     double max;
     if(a >= b) max = a;
     else max = b;

     if(max >= c) return(max);
     else return(c);
}

double min(double a, double b, double c)
{
     double min;
     if(a <= b) min = a;
     else min = b;

     if(min <= c) return(min);
     else return(c);
}
