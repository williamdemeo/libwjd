/************************************************************
 * regress.c  main program for performing regression        *
 *                                                          *
 * Created by William J. De Meo                             *
 * on 11/28/97                                              *
 *                                                          *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "prototypes.h"
#include "clapack.h"
#define MAX_NAME 100

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 BLAS Subroutines prototypes */

/*
void dcopy_(long *N, double *X, long *INCX, double *Y, long *INCY); // y <- x
double ddot_(long *N, double *X, long *INCX, double *Y, long *INCY); //  returns x * y 
void dtrsv_(char *UPLO, char *TRANSA, char *DIAG, long *N, double *A, 
            long *LDA, double *Y, long *INCY);  // y <- inv(A)*y

// y <- (alpha)Ax + (beta)y   (or A^t if TRANSA='T')
void dgemv_(char *TRANSA, long *M, long *N, double *alpha, double *A, long *LDA, 
            double *x, long *INCX, double *beta, double *y, long *INCY);

// C <- (alpha)AB + (beta)C 
void dgemm_(char *TRANSA, char *TRANSB, long *M, long *N, long *K, double *alpha, 
            double *A, long *LDA, double *B, long *LDB, double *beta, double *C,long *LDC); 

// B <- alpha*inv(A)*B
void dtrsm_(char *SIDE, char *UPLO, char *TRANSA, char *DIAG, long *M, long *N,
            double *alpha, double *A, long *LDA, double *B, long *LDB);
*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void reg(long M, long N, double *QR, double *leadu, double *E, 
         double *y, double *B, double *cov, double *se, double *e, double *sigma);

void read_name(char *);

main()
{
     char *filename;
     double *x, *y, *leadu, *E, *B, *cov, *se, *e, *sigma;
     long i, j, nrow, ncol, mindim;
  
     filename = cmalloc(MAX_NAME);

/* matrix must be of the form [X,y] where first column of X is 
   a vector of 1's if an intercept term is desired */
     printf("\n%s\n%s","Enter file name containing the matrix [X,y] ",
            "(with intercept column included if desired): ");
     read_name(filename);
     printf("\nEnter the number of observations: ");
     scanf("%ld",&nrow);
     printf("\nYou entered %ld \n",nrow);
     printf("\nEnter the number of parameters (including intercept): ");
     scanf("%ld",&ncol);
     printf("\nYou entered %ld \n",ncol);
     if(nrow <= ncol) 
          printf("\n\nWARNING: #obs = %ld <= %ld = #parameters\n\n",nrow,ncol);
     mindim = lmin(nrow,ncol); /* mindim is the smaller dimension */

     x = dmalloc(nrow*(ncol+1));
     
     y = x+(nrow*ncol); /* y is assigned the address of last col of x */

     leadu = dmalloc(mindim);
     E = dmalloc(ncol*ncol);
     B = dmalloc(ncol);
     cov = dmalloc(ncol*ncol);
     se = dmalloc(ncol);
     e = dmalloc(nrow);
     sigma = dmalloc((long)1);
     
     //matlabread(x, nrow, ncol, filename);
     matread(x, nrow, ncol+1, filename); 
     matprint(x, nrow, ncol+1);
     /*matrix is stored contiguously column-wise */

     qrpivot(nrow,ncol,x,E,leadu); /* only send first ncol columns of x */

     printf("\nRough estimate of smallest singular value of X:");
     printf("\nR(%ld,%ld) = %lf",ncol,ncol,x[nrow*(ncol-1)+(ncol-1)]);

     reg(nrow,ncol,x,leadu,E,y,B,cov,se,e,sigma);

     printf("\n\nMSE = %lf\n",*sigma);
     printf("\nCOEFFICIENT \t SE \n");
     for (i = 0; i < ncol; i++)
          printf("%4.5lf \t %4.5lf\n", B[i],se[i]);
     printf("\n\nOBS \t RESIDUAL\n");
     for (i = 0; i < nrow; i++)
          printf("%ld \t %4.5lf\n", i+1,e[i]);
}
