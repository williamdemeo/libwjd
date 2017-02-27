/* utilities.h 
 * header file for utilities.c and utilities.cc programs, which contain
 * utility functions (in C and C++ respectively).
 */
#ifndef _UTILITIES_H_
#define _UTILITIES_H_
/**@name utilities
 * various utility functions
 * @author  William DeMeo
 * \URL[williamdemeo@gmail.com]{mailto:williamdemeo@gmail.com}
 */
//@{
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* IF TEMPLATES ARE AVAILABLE, INCLUDE THE FOLLOWING LINE */
#define TEMPLATES

/*----------------------------------------------------------------------
Functions (some available only if TEMPLATES variable is defined)
0) tmalloc: allocates memory for arrays of arbitrary type
1) dmalloc: allocates memory for arrays containing doubles 
2) ddmalloc: allocates memory for 2d arrays containing doubles 
3) cmalloc: allocates memory for pointers to characters 
4) lmalloc: allocates memory for pointers to long
5) nextpow2(N): returns the first P such that 2^P >= abs(N)
6) round: returns a rounded version of the scalar argument
7) min: returns minimum of two arguments of arbitrary type
8) max: returns maximum of two arguments of arbitrary type
9) dmin: returns minimum of two doubles
10) dmax: returns maximum of two doubles
11) umin: returns minimum of two usigned longs
12) umax: returns maximum of two usigned longs
13) lmin: returns minimum of two longs
14) lmax: returns maximum of two longs
15) d3min: returns minimum of three doubles
16) d3max: returns maximum of three doubles
17) macheps: returns the smallest double eps such that eps + 1 > 1
18) max_short: returns the largest representable short
19) max_unsigned_short: returns the largest representable unsigned short
20) max_long: returns the largest representable long
21) max_unsigned_long: returns the largest representable unsigned long
22) max_single_float: returns the largest representable float
23) max_double_float: returns the largest representable double
======================================================================*/

#ifdef TEMPLATES
#include <typeinfo> /* for typeid() */

/** allocates memory for arrays of arbitrary type. <BR>
 *  can be used to store matrices by row or column.  <BR>
 *  refer to element A[i][j] as A[i*TDA+j] or A[j*LDA+i] <BR>
 *  where LDA and TDA are the leading and trailing dimensions of A, resp.
 *  @param N  number of elements in the array 
 */
template <class Type> Type *tmalloc(long N);
#endif

/** allocates memory for arrays containing doubles. <BR>
 *  can be used to store matrices by row or column.  <BR>
 *  refer to element A[i][j] as A[i*TDA+j] or A[j*LDA+i] <BR>
 *  where LDA and TDA are the leading and trailing dimensions of A, resp.
 *  @param N  number of elements in the array 
 */
double *dmalloc(long N);

/** allocates memory for 2d arrays containing doubles 
 *  @param  LD  number of rows (leading dimension)
 *  @param  TD  number of columns (trailing dimension)
 */
double **ddmalloc(long LD, long TD);

/** allocates memory for pointers to characters 
 *  @param N  number of characters
 */
char *cmalloc(long N);

/** allocates memory for pointers to long
 *  @param N  number of longs
 */
long *lmalloc(long N);

/** for given argument N, returns P such that 2^P >= abs(N)
 */
int nextpow2( const double& N );

#ifdef TEMPLATES
/** for given argument N, returns P such that 2^P >= abs(N)
 */
template <class Type> inline int nextpow2( const Type& N );

/** returns a rounded version of a scalar argument
 * @param t1  the number to round
 */
template <class Type> inline const int& round( const Type& t1);

/** returns minimum of two arguments of arbitrary type
 *  @param t1  number to be compared with second arg
 *  @param t2  number to be compared with first arg
 */
template <class Type> inline const Type& min( const Type& t1, const Type& t2);

/** returns maximum of two arguments of arbitrary type
 *  @param t1  number to be compared with second arg
 *  @param t2  number to be compared with first arg
 */
template <class Type> inline const Type& max( const Type& t1, const Type& t2);
#endif

/** returns minimum of two doubles
 *  @param a  double to be compared with second arg
 *  @param b  double to be compared with first arg
 */
inline const double& dmin(const double& a, const double& b);

/* returns maximum of two doubles
 * @param a  double to be compared with second arg
 * @param b  double to be compared with first arg
 */
inline const double& dmax(const double& a, const double& b);

/** returns minimum of two unsigned longs
 *  @param a  unsigned long to be compared with second arg
 *  @param b  unsigned long to be compared with first arg
 */
inline const unsigned long& umin(const unsigned long& a, const unsigned long& b);

/** returns maximum of two unsigned longs
 *  @param a  unsigned long to be compared with second arg
 *  @param b  unsigned long to be compared with first arg
 */
inline const unsigned long& umax(const unsigned long& a, const unsigned long& b);

/** returns minimum of two longs
 *  @param a  long to be compared with second arg
 *  @param b  long to be compared with first arg
 */
inline const long& lmin(const long& a, const long& b);

/** returns maximum of two longs
 *  @param a  long to be compared with second arg
 *  @param b  long to be compared with first arg
 */
inline const long& lmax(const long& a, const long& b);

/** returns minimum of three doubles
 * @param a  double to be compared with other args
 * @param b  double to be compared with other args
 * @param c  double to be compared with other args
 */
inline const double& d3min(const double& a, const double& b, const double& c);

/** returns maximum of three doubles
 * @param a  double to be compared with other args
 * @param b  double to be compared with other args
 * @param c  double to be compared with other args
 */
inline const double& d3max(const double& a, const double& b, const double& c);

/** returns the smallest double eps such that eps + 1 > 1 */
double macheps();

/** return largest short */
short max_short();

/** return largest unsigned short */
unsigned short max_unsigned_short();

/** return largest long */
long max_long();

/** return largest unsigned long */
unsigned long max_unsigned_long();

/** return largest single precision float */
float max_single_float();

/** return largest double precision float */
double max_double_float();


//@}
#endif
