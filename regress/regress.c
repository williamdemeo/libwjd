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


/* Subroutine reg()
   Arguments:
           
        M number of rows of X
 
        N number of columns of X (expect N < M)

        QR the matrix resulting from applying qrpivot() to X

        leadu
              on entry: the vector of leading u's resulting from qrpivot()
              on exit: the vector of coefficient estimates B, where y = XB

        E the permutation matrix resulting from qrpivot()

        y  a vector (length M) of "observables" (the rhs in XB = y)

        B     on entry: an arbitrary length N vector
              on exit: the coefficient estimates
                     
        cov   on entry: an arbitrary NxN matrix
              on exit: the covariance matrix

        se    on entry: an arbitrary length N vector
              on exit: the s.e.'s of the coefficient estimates

        e     on entry: an arbitrary length M vector
              on exit: the vector of residuals:  e = y - XB
           
        sigma   on exit: the mse =  y^te / (M-N)
           */
        
void reg(long M, long N, double *QR, double *leadu, double *E, 
         double *y, double *B, double *cov, double *se, double *e, double *sigma)
{

     long i,j,mindim;
     double a, *Qy, *invR, *EiR, *coef;

     /* BLAS arguments */
     long INC=(long)1;
     double alpha = (double)1, beta = (double)0;
     char UPLO, NOTRANS, TRANS, DIAG, SIDE; 
     UPLO='U'; NOTRANS = 'N'; TRANS = 'T'; DIAG = 'N'; SIDE='L';
  
     dcopy_(&M, y, &INC, e, &INC);     /* e <- y */

     mindim = lmin(M-1,N);    /* expect mindim = N */

     /* Apply P(n)...P(1) to e to get e <- (Q_1 Q_2)^t Y*/
     for(j=0;j<mindim;j++)
     {
          a = leadu[j]*e[j];        /* initialize  a = u(1)e(1) */
          for(i=j+1;i<M;i++)
               a += QR[M*j+i]*e[i]; /* a = u^t e */
          a *= (double)(-2);
          e[j] += a * leadu[j];     /* e(1) <- e(1) - 2 u(1)u^te */
          for(i=j+1;i<M;i++)
               e[i] +=  a* QR[M*j+i]; /* e <- e + (-2)uu^t e */
     }

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  COMPUTE COEFFICIENTS
  
  BLAS 3 method:        (currently used method) */

     /* compute inv(R) */
     invR = dmalloc(N*N);     /* workspace */
     for(j=0;j<N;j++)         /* begin with identity matrix */
     {
          for(i=0;i<N;i++)
               invR[N*j+i]=(double)0;
          invR[N*j+j]=(double)1;
     }
     /* invR <- alpha*inv(R)*invR = alpha*inv(R)*eye  */
     dtrsm_(&SIDE, &UPLO, &NOTRANS, &DIAG, &N, &N, &alpha, QR, &M,invR,&N);

     /* compute the E*inv(R) matrix */
     EiR = dmalloc(N*N);      
     /*  EiR <- (alpha)E*invR + (beta)EiR  */
     dgemm_(&NOTRANS, &NOTRANS, &N, &N, &N, &alpha, E, &N, 
            invR, &N, &beta, EiR, &N); 
     free(invR);

     /* compute the coefficients */
     dgemv_(&NOTRANS, &N, &N, &alpha, EiR, &N, e, &INC, &beta, B, &INC);
     /* B <- (alpha)EiR*e + (beta)B    (beta = 0)  only references 1st N elements of e */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
     
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  COMPUTE COEFFICIENTS

  Alternative method (BLAS 2):       */

     if(0) /* not currently used */
     {
          coef = dmalloc(N);  /* workspace */
          dcopy_(&N, e, &INC, coef, &INC); /* coef <- e(1:N) */
          dtrsv_(&UPLO, &NOTRANS, &DIAG, &N, QR, &M, coef, &INC);  /* coef <- Inv(R)*coef */  
          dgemv_(&NOTRANS, &N, &N, &alpha, E, &N, coef, &INC, &beta, B, &INC);
          /* B <- (alpha)E*coef + (beta)B    (beta = 0)  */
          free(coef);
     }

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
     
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  COMPUTE RESIDUALS
  */
     for(i=0;i<N;i++) e[i]=(double)0; /* annihilate first N elements of e */

     /* Apply P(1)...P(n) to e  to get e <- Q2 Q2^t Y*/
     for(j=(mindim-1);j>=0;j--)
     {
          a = leadu[j]*e[j];        /* initialize  a = u(1)e(1) */
          for(i=j+1;i<M;i++)
               a += QR[M*j+i]*e[i]; /* a = u^t e */
          a *= (double)(-2);
          e[j] += a * leadu[j];     /* e(1) <- e(1) - 2 u(1)u^te */
          for(i=j+1;i<M;i++)
               e[i] +=  a* QR[M*j+i]; /* e <- e + (-2)uu^t e */
     }
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  COMPUTE MSE
  */
     *sigma = ddot_(&M, y, &INC, e, &INC);
     *sigma /= (M - N);

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  COMPUTE COVARIANCE MATRIX and SE's

  cov <- (alpha)EiR*(EiR)' + (beta)cov  */
     dgemm_(&NOTRANS, &TRANS, &N, &N, &N, &alpha, EiR, &N, 
            EiR, &N, &beta, cov, &N); 

     for(j=0;j<N;j++)
          se[j] = sqrt((*sigma)*cov[N*j+j]);
}

void read_name(char *name)
{
     int c, i = 0;
  
     while ((c = getchar()) != EOF && c != ' ' && c != '\n')
          name[i++] = c;
     name[i] = '\0';
}

  
