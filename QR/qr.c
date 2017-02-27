/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  qr.c

  Created on 20110819 by williamdemeo@gmail.com based on file House.c.
  Last modified on 20110819.

  Purpose: QR decomposition of an m-by-n matrix using
           Householder reflections
           
  Further Details:  This implementation uses BLAS 2
                    (matrix-vector mult. and rank 1 updates)

  Dependencies:  Requires subroutines found in the libraries:
                 libf77blas, libatlas, libgfortran, libm
                 the first two are distributed with Atlas
                 The linking order is as follows:
                 -lf77blas -lgfortran -latlas -lm

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include "prototypes.h"
#include "clapack.h"
#include <stdlib.h>
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
BLAS Subroutine prototypes                                            
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*double dnrm2_(long *N, double *x, long *INC);
void dcopy_(long *N, double *X, long *INCX, double *Y, long *INCY);
void dgemv_(char *TRANSA, long *M, long *N, double *alpha, double *A, long *LDA, 
  double *x, long *INCX, double *beta, double *y, long *INCY);

void dger_(long *M,long *N,double *alpha,double *x,long *INCX,double *y,long *INCY,
  double *A,long *LDA);
//void dswap_(const int N, double *X, const int incX, double *Y, const int incY);
void dswap_(long *N, double *X, long *incX, double *Y, long *incY);
*/
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/* Subroutine qr:

   Arguments: 
              M  (long) number of rows of A

              N  (long) number of columns of A

              A  (pointer to double) 
                 on entry: the M by N matrix to be decomposed
                 on exit: upper-right-triangle = R
                 column i of lower trapezoid = a(i+1:M,i) = u(1:m,i)
                 where P(i) = I - 2 u(1:m,i)u^t(1:m,i) is the ith
                 Householder transformation (of dimension (M-i)x(M-i))

              leadu  (pointer to doulbe)
                     on entry: an arbitrary vector of length min(M-1,N)
                     on exit: the leading entries of the Householder vectors u(i) 
                              i.e. u(i) = (leadu(i), a(i+1:M,i)) i=1,...,,min(M-1,N)

Note that Q is obtained from augmenting the Householder tranformations
to be of proper dimensions, and then multiplying:
If P'(i) denotes augmented P(i),

Q = P'(1) P'(2) P'(3) ... =
            
|        |  |1|  0   |  |1   |     |
|  P(1)  |  |--------|  |   1|  0  |
|        |  |0| P(2) |  |----------|  ...  
|        |  | |      |  |  0 | P(3)|   

But this is left to the calling function and is not performed in qr().
*/

void qr(long M, long N, double *A, double *leadu)
{
  char T = 'T';  
  long i, nrow, ncol, mindim;
  double unit, zero, alpha=(double)-2;
  double *u, *y;  /* used for temporary work space */
  long INC=1;      /* INC is used to represent storage 
                     spacing between elements */
  unit = (double)1; zero = (double)0;

  u = dmalloc(M);  /* work space */
  y = dmalloc(N);

  mindim = lmin(M-1,N);

  for(i=0;i<mindim;i++)
    {
      nrow=M-i;
      ncol=N-i;

      House(nrow, A+(M*i+i), u);

      /* y <- (A^t)u  (is working)*/ 
      dgemv_(&T, &nrow, &ncol, &unit, A+(M*i)+i, 
	      &M,u,&INC,&zero,y,&INC);
        
      /* Rank 1 update: A <- A + (-2)uy^t   i.e. A - 2uu^tA */
      dger_(&nrow,&ncol, &alpha, u, &INC, y, &INC, A+(M*i)+i, &M);
      leadu[i] = u[0];

      /* store u(2:nrow) in A(i+1:M,i) */
      nrow--; 
      dcopy_(&nrow, u+1, &INC, A+(M*i)+i+1,&INC);

    }
        free(u); free(y);
}

/* Subroutine qrpivot:

   Arguments: same as qr() with one exception:

              E (pointer to double) 
                on entry: an N by N matrix of zero's
                on exit: the permutation matrix
                The final decomposition is AE = QR
*/

void qrpivot(long M, long N, double *A, double *E, double *leadu)
{
     char T = 'T';  
     long i,j, nrow, ncol, mindim, perm=0;
     double unit, zero, alpha=(double)-2, maxnorm, norm;
     double *u, *y;  /* used for temporary work space */
     long INC=1;      /* INC is used to represent storage 
                         spacing between elements */
     unit = (double)1; zero = (double)0;

     u = dmalloc(M);  /* work space */
     y = dmalloc(N);
     /* Start permutation matrix as the identity */
     for(j=0;j<N;j++)
       for(i=0;i<N;i++)
         E[N*j+i] = (double)0;
     for(i=0;i<N;i++)
       E[N*i+i] = (double)1;

     mindim = lmin(M-1,N);

     for(i=0;i<mindim;i++)
     {
          nrow=M-i;
          ncol=N-i;

          /* column pivot */
          maxnorm=0;
          for(j=i;j<N;j++)
          {
               norm = dnrm2_(&nrow,A+(M*j+i),&INC);
               if(norm>maxnorm)
               {
                    perm=j;
                    maxnorm=norm;
               }
          }
          
          if(perm>i)
          { /* If the i'th column was not the largest in norm, permute cols
               of A, and note it by swapping cols of pivot matrix*/
	    dswap_(&M,A+(M*i),&INC,A+(M*perm),&INC);	    /*               dswap_(&M,A+(M*i),&INC,A+(M*perm),&INC);*/
	    dswap_(&N,E+(N*i),&INC,E+(N*perm),&INC);
	    /*	       dswap_(&N,E+(N*i),&INC,E+(N*perm),&INC);*/
          }

          House(nrow, A+(M*i+i), u);

          /* y <- (A^t)u  */ 
          dgemv_(&T, &nrow, &ncol, &unit, A+(M*i)+i, 
		 &M,u,&INC,&zero,y,&INC);
        
          /* Rank 1 update: A <- A + (-2)uy^t   i.e. A - 2uu^tA */
          dger_(&nrow,&ncol, &alpha, u, &INC, y, &INC, A+(M*i)+i, &M);
          leadu[i] = u[0];

          /* store u(2:nrow) in A(i+1:M,i) */
            nrow--; 
	    dcopy_(&nrow, u+1, &INC, A+(M*i)+i+1,&INC);

     }
          free(u); free(y);
}


void House(long N, double *a, double *u)
{
  double sign = (double)1, norm;
  long inc=1;
  long i;
    
  if(a[0] < (double)0) sign = -1;

  /* u <- A(i:M,i) */
  dcopy_(&N, a, &inc, u, &inc);

  norm = dnrm2_(&N, a, &inc);       /* L2 norm of a = A(i:M,i) */
  u[0] += sign*norm;      
  norm = dnrm2_(&N, u, &inc);        /* L2 norm of new u */
  /* could also try the relation: norm = sqrt(2*(norma*norma + sign*u[1]*norma)) */

  for(i=0;i<N;i++)  u[i] /= norm;
  /* consider skipping this normalization and putting the norms in the coefficient:
     alpha = 2/(unorm * unorm) */
}

  




