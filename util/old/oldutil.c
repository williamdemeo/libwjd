/*======================================================================
  memalloc.c

  Created by William J. De Meo
    on 10/01/97

  Purpose:  define functions for dynamic memory allocation
======================================================================*/
#include <stdio.h>
#include <stdlib.h>
/* #include <prototypes.h> */
/*
  readdata()
    purpose: Place all data values of a file into an nx1 array
    arguments:
               ifp = pointer to the input data file
               n = total number of observations in the whole file
               byrow = indicator: T -> data stored by row
                                  F -> data stored by column
*/

void prn_info(char *prog_name)
{
  printf("\n%s%s%s\n\n",
         "Usage:  ",prog_name," infile  outfile");
}

double *readdata(FILE *ifp, unsigned long n, int byrow)
{
  double *x;
  char d;
  
unsigned long i;
  
    x = dmalloc(n);
    for(i=0;i<n;i++)x[i]=0;
    
  if(byrow == 0){
    while((d = getc(ifp)) != EOF)
          {
               while(d == ' ' || d == '\n'){
                    d = getc(ifp); /* skip white space and newlines */
               }
          *x++ += (double)d;
          }
  }
  else if(byrow == 1){
  }

else{
printf("\nUsage: readdata(): invalid third argument.");
exit(1);
                }
  return(x);
}


/*----------------------------------------------------------------------
function: dmalloc() 
purpose: allocates memory for arrays containing doubles 
arguments: unsigned long N = number of elements in the array 
----------------------------------------------------------------------*/

double *dmalloc(unsigned long N)
{
     double *x;

     if ((x=(double*) malloc((unsigned)(N * sizeof(double)))) == NULL)
     {
	  printf("Error: Couldn't get %ul bytes of memory \n", 
		 (unsigned)(N*sizeof(double)));
	  exit(1);
     }
     return(x);
}


/*----------------------------------------------------------------------
function: ddmalloc() 
purpose: allocates memory for 2d arrays containing doubles 
arguments: unsigned long M = number of rows, 
           unsigned long N = number of columns,
remarks: we can refer to element A[i][j] as A[j*M+i]
----------------------------------------------------------------------*/

double **ddmalloc(unsigned long M, unsigned long N)
{
     double *At;
     double **A;
     int i;

     /* Allocate all space needed for the array at once */
     /* so memory is contiguously allocated */
     if((A=(double**) malloc((unsigned)(M * sizeof(double*)))) == NULL) 
     {
	  printf("Error: Couldn't get %ul bytes of memory \n",
		 (unsigned)(M * sizeof(double*)));
	  exit(1);
     }

     /* Then give all the space to row 1 */
     A[0]= dmalloc(M*N);
     At = A[0];

     /* Then divide it up among rows 2 thru M */
     for(i=1;i<M;i++) A[i] = (At += N);

     return(A);
}


