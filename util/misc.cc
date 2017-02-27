/**@name misc.cc
 * miscellaneous utility routines (c++ version)
 * @see misc.h, misc.c
 * @author  William DeMeo
 * \URL[william.demeo@verizon.net]{mailto:william.demeo@verizon.net}
 */
//@{
/* miscellaneous useful functions */

/* (1)   matprint -- print a matrix of doubles                       */
/* (2)   lmatprint -- print a matrix of longs                        */
/* (3)   matcopy -- copy a matrix of doubles                         */
/* (4.a) matread -- read an n*p matrix of doubles from a file        */
/* (4.b) matlabread -- read a matrix stored column-wise in a file    */
/* (5)   matwrite -- write an n*p matrix of doubles to a file        */
/* (6)   check  -- checks file opened properly                       */
/* (7)   invert_ut -- invert an upper triagular matrix               */
/* (8)   x_t_y -- matrix mult of X Transpose * Y                     */
/* (9)   x_t_x -- matrix mult of X Transpose * X                     */
/* (10)  factorial -- compute n!                                     */
/* (11)  binomial -- compute binomial coeficient "n choose x"        */
/* (12)  kron -- find kronecker product of two vectors               */
/* (13)  save_string -- return copy of char* arg  (Stroustrup p.128) */

#include "misc.h"

// set USING_EXCEPTIONS to 0 for old method of handling errors
#define USING_EXCEPTIONS  1

/* (1) */
/** prints a matrix of doubles
 *  @param x  holds the matrix data
 *  @param nrow  number of rows in x (leading dimension)
 *  @param ncol  number of columns in x (trailing dimension)
 */
void    matprint(double *x, long nrow, long ncol)
{
	long    i, j;

	for (i = 0; i < nrow; i++)
	{
		for (j = 0; j < ncol; j++)
			printf("%4.5lf \t", x[nrow * j + i]);
		printf("\n");
	}
}

/* (2) */
/** prints a matrix of longs
 *  @param x  holds the matrix data
 *  @param nrow  number of rows in x (leading dimension)
 *  @param ncol  number of columns in x (trailing dimension)
 */
void    lmatprint(const long *x, long nrow, long ncol)
{
	long    i, j;

	for (i = 0; i < nrow; i++)
	{
		for (j = 0; j < ncol; j++)
			printf("%ld \t", x[nrow * j + i]);
		printf("\n");
	}
}

/* (3) */
/** copy a matrix of doubles
 *  @param from  holds the matrix to be copied (the source)
 *  @param to  the destination
 *  @param nrow  number of rows (leading dimension)
 *  @param ncol  number of columns (trailing dimension)
 */
void    matcopy(const double *from, double *to, long nrow, long ncol)
{
	long    i;

	for (i = 0; i < (nrow * ncol); i++)
		to[i] = from[i];
}

/* (4.a) */
/** read a matrix and store it column-wise
 *  @param x  where the matrix will be stored
 *  @param nrow  number of rows (leading dimension)
 *  @param ncol  number of columns (trailing dimension)
 *  @param file_to_read  string holding the name of the file to read
 */
void matread(double *x, long nrow, long ncol, const char *file_to_read) throw(char *)
{
  void    check(FILE *);
  long    c, i;
  FILE   *infile;
  
  infile = fopen(file_to_read, "r");

  if(USING_EXCEPTIONS && infile==NULL) 
    throw "matread: could not open file for writing";
  else
    check(infile);
  
  for (i = 0; i < nrow; i++) {
    for (c = 0; c < ncol; c++)
      fscanf(infile, "%lf", x + nrow * c + i);
  }
  fclose(infile);
}

/* (4.b) */
/** read a matrix stored column-wise and store it column-wise
 *  @param x  where the matrix will be stored
 *  @param nrow  number of rows (leading dimension)
 *  @param ncol  number of columns (trailing dimension)
 *  @param file_to_read  string holding the name of the file to read
 */
void matlabread(double *x, long nrow, long ncol, char *file_to_read) throw(char *)
{
  void    check(FILE *);
  long    i;
  FILE   *infile;

  infile = fopen(file_to_read, "r");

  if(USING_EXCEPTIONS && infile==NULL) 
    throw "matread: could not open file for writing";
  else
    check(infile);

  for (i = 0; i < (nrow*ncol); i++)
    fscanf(infile, "%lf", x + i);
  fclose(infile);
}

/* (5) */
/** write an nrow x ncol matrix of doubles to a file
 *  @param x  holds the matrix data
 *  @param nrow  number of rows (leading dimension)
 *  @param ncol  number of columns (trailing dimension)
 *  @param file_to_write  string holding the name of the file to write
 */
void matwrite(const double *x, long nrow, long ncol, const char *file_to_write) throw(char * )
{
    void    check(FILE *);
    long    c, i;
    FILE   *outfile;

    outfile = fopen(file_to_write, "w");

    if(USING_EXCEPTIONS && outfile == NULL) 
      throw "matwrite: could not open file for writing";
    else
      check(outfile);

    for (i = 0; i < nrow; i++)
    {
        for(c = 0; c < ncol; c++)
            fprintf(outfile, "%4.5lf   ", x[nrow * c + i]);
        fprintf(outfile, "\n");
    }
    fclose(outfile);
}    

/* (6) */
/** check if we can open a given file
 *  @param a  file pointer
 */
void    check(FILE * a)
{
  if (a == NULL) {   
    fprintf(stderr,"null pointer exception\n");
    exit(1);  
  }
}

/* (7) */
/** invert an upper triagular matrix 
 *  @param t  p x p matrix to be inverted
 *  @param u  p x p array which, on exit, holds the inverted matrix
 *  @param p  is the number of rows of the (square) matrices
 */
void    invert_ut(double *t, double *u, long p)
{				/* inverse is u */
  long    i, j, k;
  double  temp;

  for (j = 0; j < p; j++) {
    /* check to see if upper triangular t is singular */
    if (t[j * p + j] < 1.e-8){
      fprintf(stderr,"singular matrix passed to 'invert_ut'\n");
      exit(2);
    }
  }
  for (j = p - 1; j >= 0; j--){
    u[j * p + j] = 1 / t[j * p + j];
    for (k = j - 1; k >= 0; k--){
      temp = 0;
      for (i = k + 1; i <= j; i++)
	temp += t[i * p + k] * u[j * p + i];
      u[j * p + k] = -1. * temp / t[k * p + k];
    }
  }
}

/* (8) */
/** matrix mult of transpose(X) * Y.
 *  x is n by p, and y is n by q matrix , stored by columns: ret. via z
 *  @param x  an n x p array
 *  @param y  an n x q array
 *  @param z  on exit, a p x q array holding x'*y
 */
void    x_t_y(double *x, double *y, double *z, long n, long p, long q)
{
	long    i, j, k;
	double  temp;

	for (i = 0; i < p; i++)
	{
		for (j = 0; j < q; j++)
		{
			temp = 0;
			for (k = 0; k < n; k++)
				temp += x[n * i + k] * y[n * j + k];
			z[p * j + i] = temp;
		}
	}
}

/* (9) */
/** matrix mult of transpose(X) * X.
 *  x is n by p, stored by columns, x transpose x returned via z
 *  @param x  an n x p array
 *  @param z  on exit, a p x p array holding x'*y
 */
void    x_t_x(double *x, double *z, long n, long p)
{
	long    i, j, k;
	double  temp;

	for (i = 0; i < p; i++)
	{
		for (j = i; j < p; j++)
		{
			temp = 0;
			for (k = 0; k < n; k++)
				temp += x[n * i + k] * x[n * j + k];
			z[p * j + i] = z[p * i + j] = temp;
		}
	}
}

/* (10) */
/** computes n facorial */
long     factorial(long n)
{
    long t, k;
    
    t = 1;
    for(k=2;k<=n; t *= k++);
    return(t);                                                           
}

/* (11) */
/** computes n choose x */
long     binomial(long n, long x)
{
    long factorial(long n);
    if(x>n)
    {
        printf("\nERROR: Parameter out of range in 'binomial(n,x)'\n");
        exit(2);
    }
    return(factorial(n)/(factorial(x) * factorial(n-x)));
}

/* (12) */
/** Find kronecker product of m by 1 vector a and n by 1 vector b.
 *  Resulting m*n by 1 product is passed back via c 
 *  @param a  an m x 1 vector
 *  @param m  length of a
 *  @param b  an n x 1 vector
 *  @param n  length of b
 *  @param c  an m x n array which, on exit, holds the kronecker product
 */
void kron(long *a, long m, long *b, long n, long *c)
{
  long i, j; 

  for(j=0;j<m;j++)
    for(i=0;i<n;i++)
      c[j*n + i] = a[j]*b[i];
}

/* (13) */
/** return a copy of the given string 
 * @param p  source string to be copied 
 * @return s  copy of p  */
char* save_string(const char* p)
{
  char* s = new char[strlen(p)+1];
  strcpy(s,p);  // from <string.h>
  return s;
}

//@}
