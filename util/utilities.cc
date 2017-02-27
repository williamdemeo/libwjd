/**@name utilities.cc
 * various utility functions
 * @author  William DeMeo
 * \URL[william.demeo@verizon.net]{mailto:william.demeo@verizon.net}
 */
//@{

/*----------------------------------------------------------------------
Functions:
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

/* If templates are available, 
   include the following line either here or in the file utilities.h  */
// #define TEMPLATES

#include "utilities.h"
#include<iostream>
using namespace std;
#ifdef TEMPLATES
/*(0)--tmalloc------------------------------------------------------*/
/** allocates memory for arrays of arbitrary type. <BR>
 *  can be used to store matrices by row or column.  <BR>
 *  refer to element A[i][j] as A[i*TDA+j] or A[j*LDA+i] <BR>
 *  where LDA and TDA are the leading and trailing dimensions of A, resp.
 *  @param N  number of elements in the array 
 */
template <class Type>
Type *tmalloc(long N)
{
  Type *x;
  x=(Type*) malloc((long)(N * sizeof(Type)));
  if (x == NULL) {
    fprintf(stderr,"couldn't allocate %l bytes of memory \n",(long)(N*sizeof(Type)));
    exit(1);
  }
  return(x);
}
#endif  

/*(1)--dmalloc------------------------------------------------------*/
/** allocates memory for arrays containing doubles. <BR>
 *  can be used to store matrices by row or column.  <BR>
 *  refer to element A[i][j] as A[i*TDA+j] or A[j*LDA+i] <BR>
 *  where LDA and TDA are the leading and trailing dimensions of A, resp.
 *  @param N  number of elements in the array 
 */
double *dmalloc(long N)
{
  double *x;
  x=(double*) malloc((long)(N * sizeof(double)));
  if (x == NULL) {
    //    fprintf(stderr,"couldn't allocate %l bytes of memory \n",(long)(N*sizeof(double)));
    cout << "couldn't allocate " << (long)(N*sizeof(double)) << " bytes of memory." << endl;
    exit(1);
  }
  return(x);
}

/*(2)--ddmalloc-----------------------------------------------------*/
/** allocates memory for 2d arrays containing doubles 
 *  @param  LD  number of rows (leading dimension)
 *  @param  TD  number of columns (trailing dimension)
 */
double **ddmalloc(long LD, long TD)
{
  double *At;
  double **A;
  int i;

  /* Allocate all space needed for the array at once */
  /* so memory is contiguously allocated */
  A=(double**) malloc((long)(LD*sizeof(double*)));
  if( A == NULL) {
    fprintf(stderr, "can't allocate %ld bytes of memory \n",(long)(LD*sizeof(double*)));
    exit(1);
  }

  /* Then give all the space to row 1 */
  A[0]= dmalloc(LD*TD);
  At = A[0];

  /* Then divide it up among rows 2 thru LD */
  for(i=1;i<LD;i++) A[i] = (At += TD);

  return(A);
}

/*(3)--cmalloc-----------------------------------------------------*/
/** allocates memory for pointers to characters 
 *  @param N  number of characters
 */
char *cmalloc(long N)
{
  char *x;
  x=(char*) malloc((long)(N * sizeof(char)));
  if(x == NULL){
    fprintf(stderr,"can't allocate %ld bytes of memory \n",(long)(N * sizeof(char*)));
    exit(1);
  }
  return(x);
}

/*(4)--lmalloc------------------------------------------------------*/
/** allocates memory for pointers to long
 *  @param N  number of longs
 */
long *lmalloc(long N)
{
  long *x;
  x=(long*) malloc((long)(N * sizeof(long)));
  if(x == NULL){
    fprintf(stderr,"can't allocate %ld bytes of memory \n",(long)(N * sizeof(long*)));
    exit(1);
  }
  return(x);
}

/*(5)--nextpow2----------------------------------------------------------*/
/** for given argument N, returns P such that 2^P >= abs(N)
 */
int nextpow2( const double& N )
{
  double absN, logN;
  double log2 = log((double)2.0);

  absN = fabs(N);
  logN = log(absN);

  return (int)ceil(logN/log2);
}

#ifdef TEMPLATES
/*(5b)--nextpow2----------------------------------------------------------*/
/** for given argument N, returns P such that 2^P >= abs(N)
 */
template <class Type>
inline int nextpow2( const Type& N )
{
  double absN, logN;
  double log2 = log((double)2.0);

  if(typeid(int) == typeid(N))
    absN = (double) abs(N);
  else if(typeid(double) == typeid(N))
    absN = fabs(N);

  logN = log(absN);

  return (int)ceil(logN/log2);
}

/*(6)--round---------------------------------------------------------*/
/** returns a rounded version of a scalar argument
 * @param t1  the number to round
 */
template <class Type>
inline const int& round( const Type& t1)
{
  return ( ( t1 < 0 ) ? (int)((t1)-0.5) : (int)((t1)+0.5) );
}

/*(7)--min----------------------------------------------------------*/
/** returns minimum of two arguments of arbitrary type
 *  @param t1  number to be compared with second arg
 *  @param t2  number to be compared with first arg
 */
template <class Type>
inline const Type& min( const Type& t1, const Type& t2)
{
  return ( ( t1 < t2 ) ? t1 : t2 );
}

/*(8)--max----------------------------------------------------------*/
/** returns maximum of two arguments of arbitrary type
 *  @param t1  number to be compared with second arg
 *  @param t2  number to be compared with first arg
 */
template <class Type>
inline const Type& max( const Type& t1, const Type& t2)
{
  return ( ( t1 > t2 ) ? t1 : t2 );
}
#endif

//#ifndef TEMPLATES
/*(9)--dmin-----------------------------------------------------------*/
/** returns minimum of two doubles
 *  @param a  double to be compared with second arg
 *  @param b  double to be compared with first arg
 */
inline const double& dmin(const double& a, const double& b)  
{
  return((a < b) ? a : b);
}

/*(10)--dmax----------------------------------------------------------*/
/* returns maximum of two doubles
 * @param a  double to be compared with second arg
 * @param b  double to be compared with first arg
 */
inline const double& dmax(const double& a, const double& b)  
{
  return((a > b) ? a : b);
}

/*(11)--umin-----------------------------------------------------------*/
/** returns minimum of two unsigned longs
 *  @param a  unsigned long to be compared with second arg
 *  @param b  unsigned long to be compared with first arg
 */
inline const unsigned long& umin(const unsigned long& a, const unsigned long& b)
{
  return((a < b) ? a : b);
}

/*(12)--umax-----------------------------------------------------------*/
/** returns maximum of two unsigned longs
 *  @param a  unsigned long to be compared with second arg
 *  @param b  unsigned long to be compared with first arg
 */
inline const unsigned long& umax(const unsigned long& a, const unsigned long& b)
{
  return((a > b) ? a : b);
}

/*(13)--lmin-----------------------------------------------------------*/
/** returns minimum of two longs
 *  @param a  long to be compared with second arg
 *  @param b  long to be compared with first arg
 */
// inline const long& lmin(const long& a, const long& b) { return((a < b) ? a : b);}
long lmin(long a, long b)
{
  return((a < b) ? a : b);
}

/*(14)--lmax-----------------------------------------------------------*/
/** returns maximum of two longs
 *  @param a  long to be compared with second arg
 *  @param b  long to be compared with first arg
 */
inline const long& lmax(const long& a, const long& b)
{
  return((a > b) ? a : b);
}

/*(15)--d3min----------------------------------------------------------*/
/** returns minimum of three doubles
 * @param a  double to be compared with other args
 * @param b  double to be compared with other args
 * @param c  double to be compared with other args
 */
inline const double& d3min(const double& a, const double& b, const double& c)
{
     const double min = ((a < b) ? a : b);

     return((min < c) ? min : c);
}

/*(16)--d3max---------------------------------------------------------*/
/** returns maximum of three doubles
 * @param a  double to be compared with other args
 * @param b  double to be compared with other args
 * @param c  double to be compared with other args
 */
inline const double& d3max(const double& a, const double& b, const double& c)
{
     const double max = ((a > b) ? a : b);

     return((max > c) ? max : c);
}

/*(17)--macheps--------------------------------------------------------*/
/** returns the smallest double eps such that eps + 1 > 1 */
double macheps()
{
  double eps, delta;
  
  for(delta = 1.0; (1.0 + delta) > 1.0; ) {
    eps=delta;
    delta /= 1.1;		
  }
  return eps;
}

/** return largest short */
short max_short()
{
  short Short, bigS;
  Short = (short) 1;
  do{
    bigS = Short;
    Short *= (short) 2;
  }while(bigS == (Short / (short) 2)); 
  /* If S*2 too big for machine, then  (S*2)/2 != S
   * and loop exits.  S stores the spurious (oversized) value.  */
  return bigS;
}

/** return largest unsigned short */
unsigned short max_unsigned_short()
{
  unsigned short UShort, bigUS;
  UShort = (unsigned short) 1;
  do{
    bigUS = UShort;
    UShort *= (unsigned short) 2;
  }while(bigUS == (UShort / (unsigned short) 2)); 

  return bigUS;
}                       

/** return largest long */
long max_long()
{
  long Long, bigL;
  Long = 1L;
  do{
    bigL = Long;
    Long *= 2L;
  }while(bigL == (Long / 2L)); 

  return bigL;
}

/** return largest unsigned long */
unsigned long max_unsigned_long()
{
  unsigned long ULong, bigUL;
  ULong = 1UL;
  do{
    bigUL = ULong;
    ULong *= 2UL;
  }while(bigUL == (ULong / 2UL)); 

  return bigUL;
}

/** return largest single precision float */
float max_single_float()
{
  float Float, bigF;
  Float = 1.0F;
  do{
    bigF = Float;
    Float *= 2.0F;
  }while(bigF == (Float / 2.0F)); 

  return bigF;
}

/** return largest double precision float */
double max_double_float()
{
  double D, bigD;
  D = (double)1.0;
  do{
    bigD = D;
    D *= (double)2.0;
  }while(bigD == (D / (double)2.0)); 

  return bigD;
}

//@}
