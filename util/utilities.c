/*======================================================================
  utilities.c

  Created by William J. De Meo
    on 10/01/97
    revised: 11/13/97

  Purpose:  define utility functions
----------------------------------------------------------------------
Functions:
1) dmalloc: allocates memory for arrays containing doubles 
2) ddmalloc: allocates memory for 2d arrays containing doubles 
3) cmalloc: allocates memory for pointers to characters 
4) lmalloc: allocates memory for pointers to long
5) nextpow2(N): returns the first P such that 2^P >= abs(N)
6) dmin: returns minimum of two doubles
7) dmax: returns maximum of two doubles
8) umin: returns minimum of two usigned longs
9) lmin: returns minimum of two longs
10) d3max: returns maximum of three doubles
11) d3min: returns minimum of three doubles
======================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*(1)-------------------------------------------------------------------
function: dmalloc() 
purpose: allocates memory for arrays containing doubles 
arguments: long N = number of elements in the array 
remarks: can be used to store matrices by column.  
         refer to element A[i][j] as A[j*M+i]
----------------------------------------------------------------------*/
double *dmalloc(long N)
{
     double *x;

     if ((x=(double*) malloc((long)(N * sizeof(double)))) == NULL)
     {
       printf("Error: Couldn't get %ld bytes of memory \n", (long)(N*sizeof(double)));
       exit(1);
     }
     return(x);
}


/*(2)-------------------------------------------------------------------
function: ddmalloc() 
purpose: allocates memory for 2d arrays containing doubles 
arguments: long M = number of rows, 
           long N = number of columns,
----------------------------------------------------------------------*/

double **ddmalloc(long M, long N)
{
     double *At;
     double **A;
     int i;

     /* Allocate all space needed for the array at once */
     /* so memory is contiguously allocated */
     if((A=(double**) malloc((long)(M * sizeof(double*)))) == NULL) 
     {
	  printf("Error: Couldn't get %ld bytes of memory \n",(long)(M * sizeof(double*)));
	  exit(1);
     }

     /* Then give all the space to row 1 */
     A[0]= dmalloc(M*N);
     At = A[0];

     /* Then divide it up among rows 2 thru M */
     for(i=1;i<M;i++) A[i] = (At += N);

     return(A);
}

/*(3)-------------------------------------------------------------------
function: cmalloc() 
purpose: allocates memory for pointers to characters 
arguments: long N = number of characters,
----------------------------------------------------------------------*/
char *cmalloc(long N)
{
  char *x;
  
  if((x=(char*) malloc((long)(N * sizeof(char)))) == NULL)
     {
	  printf("Error: Couldn't get %ld bytes of memory \n", (long)(N * sizeof(char*)));
	  exit(1);
     }
  return(x);
}

/*(4)-------------------------------------------------------------------
function: lmalloc() 
purpose: allocates memory for pointers to long
arguments: long N = number of longs
----------------------------------------------------------------------*/
long *lmalloc(long N)
{
  long *x;
  
  if((x=(long*) malloc((long)(N * sizeof(long)))) == NULL)
     {
	  printf("Error: Couldn't get %ld bytes of memory \n", (long)(N * sizeof(long*)));
	  exit(1);
     }
  return(x);
}

/*(5)--nextpow2----------------------------------------------------------
function: nextpow2()
purpose: for given argument N, returns P such that 2^P >= abs(N)
----------------------------------------------------------------------*/
/*
inline int nextpow2( const double& N )
{
  double absN, logN;
  double log2 = log((double)2.0);

  absN = fabs(N);

  logN = log(absN);

  return (int)ceil(logN/log2);
}
*/
/*(6)-------------------------------------------------------------------
function: dmin() 
purpose: returns minimum of two doubles
arguments: double a, double b are the doubles to be compared
----------------------------------------------------------------------*/
double dmin(double a, double b)  
{
  return((a < b) ? a : b);
}
/*(7)-------------------------------------------------------------------
function: dmax() 
purpose: returns maximum of two doubles
arguments: double a, double b are the doubles to be compared
----------------------------------------------------------------------*/
double dmax(double a, double b)  
{
  return((a > b) ? a : b);
}

/*(8)-------------------------------------------------------------------
function: umin() 
purpose: returns minimum of two unsigned longs
arguments: unsigned long a, unsigned long b are the numbers to be compared
----------------------------------------------------------------------*/
unsigned long umin(unsigned long a, unsigned long b)
{
  return((a < b) ? a : b);
}

/*(9)-------------------------------------------------------------------
function: lmin() 
purpose: returns minimum of two longs
arguments: long a, long b are the numbers to be compared
----------------------------------------------------------------------*/
long lmin(long a, long b)
{
  return((a < b) ? a : b);
}

/*(10)-------------------------------------------------------------------
function: d3min() 
purpose: returns maximum of three doubles
arguments: double a, double b, double c are the numbers to be compared
----------------------------------------------------------------------*/
double d3max(double a, double b, double c)
{
     double max;
     max = ((a > b) ? a : b);

     return((max > c) ? max : c);
}

/*(11)-------------------------------------------------------------------
function: d3min() 
purpose: returns minimum of three doubles
arguments: double a, double b, double c are the numbers to be compared
----------------------------------------------------------------------*/
double d3min(double a, double b, double c)
{
     double min;
     min = ((a < b) ? a : b);

     return((min < c) ? min : c);
}

