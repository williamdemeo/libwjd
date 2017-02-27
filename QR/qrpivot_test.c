/************************************************************
 * qrpivot_test.c  main program for testing QR              *
 * Orthogonalization routine qrpivot()                      *
 *                                                          *
 * Created 2001.11.01 by William J. DeMeo.                  *
 * Modified 2011.08.19.                                     *
 *                                                          *
 * Note: differences between qrpivot_test.c and qr_test.c   *
 *       are marked with the comment "qrpivot_test.c"       *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "prototypes.h"
#include "timer.h"
#define MAX_NAME 100

void qrpivot(long M, long N, double *A, double *E, double *leadu); 
void read_name(char *);

main()
{
     char *filename;
     double *x, *leadu, *E;
     long i, j, nrow, ncol, mindim;
  
     filename = cmalloc(MAX_NAME);

     printf("\nEnter file name containing the matrix: ");
     read_name(filename);
     printf("\nEnter the number of rows: ");
     scanf("%lu",&nrow);
     printf("\nEnter the number of columns: ");
     scanf("%lu",&ncol);

     mindim = lmin(nrow,ncol); /* mindim is the smaller dimension */

     x = dmalloc(nrow*ncol);
     leadu = dmalloc(mindim);
     E = dmalloc(ncol*ncol);
  
     matlabread(x, nrow, ncol, filename); 
     /*matrix is stored contiguously column-wise */

     TIME0  /* timing qr routine */
     /* Test qrpivot:   */
     qrpivot(nrow,ncol,x,E,leadu);                /* qrpivot_test.c */
     TIME1("time for C qr factorization (pivot)");

     printf("\nThe orthogonalization produced: \n");
     matprint(x,nrow,ncol);
     printf("\nWith permutation matrix: \n");     /* qrpivot_test.c */
     matprint(E,ncol,ncol);                       /* qrpivot_test.c */
     printf("\nand leading u's:\n");
     for(i=0;i<mindim;i++)
          printf("%4.5lf \t", leadu[i]);
     printf("\n");
}
  
void read_name(char *name)
{
     int c, i = 0;
  
     while ((c = getchar()) != EOF && c != ' ' && c != '\n')
          name[i++] = c;
     name[i] = '\0';
}

  
