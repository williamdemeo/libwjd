/* misc.h
 * header file for misc.c and misc.cc programs, which contain 
 * miscellaneous utility functions (in C and C++ respectively).
 */
#ifndef _MISC_H_
#define _MISC_H_
/**@name misc
 * miscellaneous utility routines
 * @author  William DeMeo
 * \URL[william.demeo@verizon.net]{mailto:william.demeo@verizon.net}
 */
//@{
#include <stdio.h>
#include <cstdlib>
#include <string.h>

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

/* (1) */
/** prints a matrix of doubles
 *  @param x  holds the matrix data
 *  @param nrow  number of rows in x (leading dimension)
 *  @param ncol  number of columns in x (trailing dimension)
 */
void matprint(const double *x, long nrow, long ncol);

/* (2) */
/** prints a matrix of longs
 *  @param x  holds the matrix data
 *  @param nrow  number of rows in x (leading dimension)
 *  @param ncol  number of columns in x (trailing dimension)
 */
void lmatprint(const long *x, long nrow, long ncol);

/* (3) */
/** copy a matrix of doubles
 *  @param from  holds the matrix to be copied (the source)
 *  @param to  the destination
 *  @param nrow  number of rows (leading dimension)
 *  @param ncol  number of columns (trailing dimension)
 */
void matcopy(const double *from, double *to, long nrow, long ncol);

/* (4.a) */
/** read a matrix and store it column-wise
 *  @param x  where the matrix will be stored
 *  @param nrow  number of rows (leading dimension)
 *  @param ncol  number of columns (trailing dimension)
 *  @param file_to_read  string holding the name of the file to read
 */
void matread(double *x, long nrow, long ncol, const char *file_to_read) throw(char *);

/* (4.b) */
/** read a matrix stored column-wise and store it column-wise
 *  @param x  where the matrix will be stored
 *  @param nrow  number of rows (leading dimension)
 *  @param ncol  number of columns (trailing dimension)
 *  @param file_to_read  string holding the name of the file to read
 */
void matlabread(double *x, long nrow, long ncol, const char *file_to_read) throw(char *);

/* (5) */
/** write an nrow x ncol matrix of doubles to a file
 *  @param x  holds the matrix data
 *  @param nrow  number of rows (leading dimension)
 *  @param ncol  number of columns (trailing dimension)
 *  @param file_to_write  string holding the name of the file to write
 */
void matwrite(const double *x, long nrow, long ncol, const char *file_to_write) throw(char *);

/* (6) */
/** check if we can open a given file
 *  @param a  file pointer
 */
void check(FILE * a);

/* (7) */
/** invert an upper triagular matrix 
 *  @param t  p x p matrix to be inverted
 *  @param u  p x p array which, on exit, holds the inverted matrix
 *  @param p  is the number of rows of the (square) matrices
 */
void    invert_ut(double *t, double *u, long p);

/* (8) */
/** matrix mult of transpose(X) * Y.
 *  x is n by p, and y is n by q matrix , stored by columns: ret. via z
 *  @param x  an n x p array
 *  @param y  an n x q array
 *  @param z  on exit, a p x q array holding x'*y
 */
void    x_t_y(double *x, double *y, double *z, long n, long p, long q);

/* (9) */
/** matrix mult of transpose(X) * X.
 *  x is n by p, stored by columns, x transpose x returned via z
 *  @param x  an n x p array
 *  @param z  on exit, a p x p array holding x'*y
 */
void    x_t_x(double *x, double *z, long n, long p);

/* (10) */
/** computes n facorial */
long     factorial(long n);

/* (11) */
/** computes n choose x */
long     binomial(long n, long x);

/* (12) */
/** Find kronecker product of m by 1 vector a and n by 1 vector b.
 *  Resulting m*n by 1 product is passed back via c 
 *  @param a  an m x 1 vector
 *  @param m  length of a
 *  @param b  an n x 1 vector
 *  @param n  length of b
 *  @param c  an m x n array which, on exit, holds the kronecker product
 */
void kron(long *a, long m, long *b, long n, long *c);

/* (13) */
/** return a copy of the given string 
 * @param p  source string to be copied 
 * @return s  copy of p  */
char* save_string(const char* p);
//@}
#endif
