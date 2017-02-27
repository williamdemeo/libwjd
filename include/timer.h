#include <sys/time.h>
#include <sys/resource.h> 

struct rusage R_usage_0;
static double R_time_0;

/*
   These macros use the system routine "getrusage" to determine how
   long a particular segment of code takes to run.  (See man getrusage
   for more information).  TIME0 initializes a variable
   to hold the current time, and TIME1 prints out a message giving 
   the number of seconds since the call to TIME0.  For example, to
   print out the time it takes to find the sum of n numbers, use:
 
        TIME0
        sum = 0.;
        for(i=0;i<n;i++)
           sum += x[i];
        TIME1("time to sum n numbers")
 
   Notice that semicolons are not required.
*/

#define TIME0      getrusage(0,&R_usage_0);\
  R_time_0=(double)R_usage_0.ru_utime.tv_sec + \
	   (double)R_usage_0.ru_utime.tv_usec / 1.e6;

#define TIME1(msg) getrusage(0,&R_usage_0);\
          printf("%s : %7.4lf sec\n",msg,\
          ((double)R_usage_0.ru_utime.tv_sec + \
	   (double)R_usage_0.ru_utime.tv_usec / 1.e6) - R_time_0);
