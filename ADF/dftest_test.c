/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  dftest_test.c

  Purpose: main program for testing routine dftest()
           for performing Dickey-Fuller test for unit root.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <math.h>
#include <stdio.h>
// #include <iostream>
// #include <cstdlib>
#include <stdlib.h>
#include <time.h>
#include "prototypes.h"
#include "dftest.h"

//using namespace std;

//void normalSequence(const long &m, double *Z);

int main()
{
  long n, m, i;
  int testH, model;
  double *y, *Z, avg, var, pVal; //, pVal=-1.0;
  DFTest_ModelType Model = DFTest_ARD;
  double alpha=0.05;  // significance level for the test. Must be between 0.001 and 0.999.
  if (Model = DFTest_AR){ model = 0; }
  if (Model = DFTest_ARD){ model = 1; }
  m= (long)100;

  printf("\n   dftest_test.c:  Getting %ld normal(0,1) random numbers...", m);
  Z = dmalloc(m); 
  normalSequence(m, Z);
  moment(m, Z, avg, var);
  printf("\n   dftest_test.c:  Average = %lf, Variance = %lf\n",avg,var);

  printf("\n   dftest_test.c:  Generating non-stationary sequence...\n");
  y = dmalloc(m); 
  y[0] = Z[0];
  for (i=1;i<m;++i){
    y[i] = y[i-1]+Z[i];
  }

  dftest(m, y, Model, alpha, testH, pVal);

  printf("\n   dftest_test.c:  ");
  printf("RESULTS\n");
  printf("\n   dftest_test.c:  testH = %d,  pVal = %lf\n", testH , pVal);
  
  free(Z);
}

