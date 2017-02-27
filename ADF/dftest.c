/************************************************************
 * dftest.c  The Dicky-Fuller test in its simplest form.    *
 *                                                          *
 * Created by William DeMeo                                 *
 * on 2011.09.25                                            *
 *                                                          *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "f2c.h"
#include "clapack.h"
#include "prototypes.h"
#include "dftest.h"
#define DEBUG_dftest 1

void reg(long M, long N, double *QR, double *leadu, double *E, 
         double *y, double *B, double *cov, double *se, double *e, double *sigma);


void dftest(long &N, double *y, const DFTest_ModelType &Model, double &alpha, int &testH, double &pVal){

  double *diffy, *x, *leadu, *E, *B, *cov, *se, *e, *sigma, *cvtable, testStat;
  long i, j, nrow, ncol, mindim, INC=1;  // n=N since blas doesn't like const argument.
  int model;

  if (Model==DFTest_AR){  model=0; }
  else if(Model==DFTest_ARD){ model=1; }
  else{   
    fprintf(stderr, "   dftest.c: unkonwn model type. Aborting... \n");
    exit(EXIT_FAILURE); 
  }

  nrow = N-1;     // Effective sample size 
  ncol = 1+model;  // if model==1, first column is a vector of 1's (for drift term)
  assert(ncol<nrow);
  mindim = ncol;

  diffy = dmalloc(N);
  dcopy_(&N, y, &INC, diffy, &INC);  // diffy <- y
  difference(N, diffy);   // first-differenced series diffy(0) = y1-y0, diffy(1) = y2-y1, etc

  x = dmalloc(nrow*(ncol+1));
  if (model==1){  // first column of x is a vector of 1's
    for (i=0; i<nrow; ++i) x[i] = 1.0;
  }
  dcopy_(&nrow, y, &INC, x+(nrow*model), &INC);         // x(:,1) =  y(0:N-2)
  dcopy_(&nrow, diffy, &INC, x+(nrow*ncol), &INC); // x(:,2) = diffy(0:N-2) (right-hand-side)
  /* The design matrix is now:      1    y0    y1-y0
                                    1    y1    y2-y1
                                    1    y2    y3-y2
				    ...                 */
#ifdef DEBUG_dftest
  char *outfile="datafile";
  printf("\n X = \n");
  matprint(x,nrow, ncol+1);
  matwrite(x, nrow, ncol+1, outfile);
#endif

  leadu = dmalloc(mindim);
  E = dmalloc(ncol*ncol);
  B = dmalloc(ncol);
  cov = dmalloc(ncol*ncol);
  se = dmalloc(ncol);
  e = dmalloc(nrow);
  sigma = dmalloc((long)1);
     
  // Perform the regression:
  qrpivot(nrow,ncol,x,E,leadu); /* only send first ncol columns of x */

  printf("\nRough estimate of smallest singular value of X:");
  printf("\nR(%ld,%ld) = %lf",ncol,ncol,x[nrow*(ncol-1)+(ncol-1)]);

  reg(nrow,ncol,x,leadu,E,diffy,B,cov,se,e,sigma);
  
  // Compute the statistic:
  testStat = (B[model]-1)/se[model];
  // Later make this a function to accommodate the various kinds of tests we could do.
  // void getStat(i,testT,testLags,testModel,testType,testReg,sigLevels,sampSizes,CVTable,needPValue);
  
#ifdef DEBUG_dftest
  printf("\n\nMSE = %lf\n",*sigma);
  printf("\nCOEFFICIENT \t SE \n");
  for (i = 0; i < ncol; i++)  printf("%4.8lf \t %4.8lf\n", B[i],se[i]);
  //printf("\n\nOBS \t RESIDUAL\n");  for (i = 0; i < nrow; i++)  printf("%ld \t %4.5lf\n", i+1,e[i]);
  printf("\n\nCov\n");  matprint(cov, ncol, ncol);
  printf("\n\nTest statistic: %lf\n", testStat);
#endif
    
  // Test the statistic and get p-value:
  runTest(alpha, testStat, nrow, Model, testH, pVal);

  
}    

/* difference(const long& N, double *x) -- return the first-differenced series: D(x) = x - Lag(x)
 * N  length of the vector x
 * x  on entry, a vector of length N
 *    on exit, a vector of length N-1 containing first differenced series:
 *    {x[0], x[1]-x[0], x[2]-x[1], ..., x[N-1]-x[N-1]}
 */
void difference(const long& N, double *x) {
  double *Lx, alpha=-1.0; 
  long INC=1, n = N-1;
  Lx = dmalloc(n);
  dcopy_(&n, x+1, &INC, Lx, &INC);        // Lx(0:N-2) <- x(1:N-1)
  daxpy_(&n, &alpha, x, &INC, Lx, &INC);  // Lx(0:N-2) <- Lx(0:N-2) - x(0:N-2)
  dcopy_(&n, Lx, &INC, x, &INC);          // x(0:N-2) <- Lx(0:N-2) = x(1:N-1)-x(0:N-2)
  free(Lx);
}
    
// Test the statistic
void runTest(const double& alpha, const double& testStat, const long& sampleSize, const  DFTest_ModelType& Model, int& testH, double& pVal){
  int row, col;  
  double testCValue=-999; // initialize to -999 so we don't accidentally reject the null if testCValue doesn't get asigned a value.

  getSampleSizeIndex(sampleSize, row);  // row is the row of the CVTable we need
  getSigIndex(alpha, col);              // col is the col of the CVTable we need
  switch(Model){
  case DFTest_AR:
    // CVTable's have size 15x50=DFTest_CVTableNumRows*DFTest_CVTableNumCols...
    testCValue = DFTest_CVTable_AR[row*DFTest_CVTableNumCols + col]; // ...stored in ROW-major format!
    getPValue(testStat, DFTest_CVTable_AR+(row*DFTest_CVTableNumCols), pVal);
    break;
  case DFTest_ARD:
    testCValue = DFTest_CVTable_ARD[row*DFTest_CVTableNumCols + col];
    getPValue(testStat, DFTest_CVTable_ARD+(row*DFTest_CVTableNumCols), pVal);
    break;
  default:
    // This error message should show what type of model was passed in.
    fprintf(stderr, "   dftest.c: unkonwn model type. Aborting... \n");
    exit(EXIT_FAILURE); 
    // How do you do that with an enum type?
    break;
  }
  
  if (testStat < testCValue)
    testH = 1;  // reject the null in favor of the alternative
  else 
    testH = 0;  // fail to reject the null
}

int getPValue(const double& testStat, const double *table, double& pVal){
  int j;
  double fd;
  if (testStat<=table[0]){
    pVal = DFTest_sigLevel[0];
    if (testStat<table[0]){
      printf("   WARNING: (dftest.c) Test statistic, %lf, is outside of table values (p-value only approximated).", testStat);
    }
    return(1);
  }
  j=0;
  while(testStat>table[j] & j<DFTest_CVTableNumCols){
    j++;
  }
  if (j==DFTest_CVTableNumCols){
    printf("   WARNING: (dftest.c) Test statistic, %lf, is outside of table values (p-value only approximated).", testStat);
    pVal = DFTest_sigLevel[j-1];
    return(1);
  }
  // At this point testStat must be somewhere between the values table[j-1] and table[j].
  // Check that:
  if (testStat < table[j-1] || testStat > table[j]){ 
    fprintf(stderr, "   ERROR: (dftest.c) unexpected test statistic... t = %lf\n", testStat);
    fprintf(stderr, "   ERROR: (dftest.c) p-value not set.    Aborting... \n");
    exit(EXIT_FAILURE); 
  }

  // Compute the "fractional distance" for the purpose of interpolating sig values:
  fd = (testStat - table[j-1])/(table[j] - table[j-1]);  
  // When fd closer to 1 (0), testStat is closer to table[j] (table[j-1]).

  // Linear interpolation (to get approximate p-value).
  pVal = (1-fd)*DFTest_sigLevel[j-1] + fd*DFTest_sigLevel[j];
  // To get the most conservative p-value estimate, we could just take pVal = table[j]. 
  // In that case, if we reported a p-value of 0.020 = table[j], it would mean the 
  // true p-value is somewhere in the interval (0.015, 0.020) = (table[j-1], table[j]).
  // This seems overly conservative, and interpolation is probably better.
  return(1);
}
      
// getSigIndex -- get index of the element in DFTest_sigLevel array (see dftest.h) 
//                that is the greatest value not above alpha.
int getSigIndex(const double& alpha, int& i){
  int j;
  i = -1;

#ifdef DEBUG_dftest
  assert(alpha>0.0009); assert(alpha<0.9999);
#endif

  if (alpha==DFTest_sigLevel[0]){ 
    i=0;
    return(1);
  }
  if (DFTest_sigLevel[0]<alpha){
    j=0;
    while (DFTest_sigLevel[j]<alpha){
      j++;
      assert(j<DFTest_sigLevelLength);
    }
    i=j-1;  // DFTest_sigLevel[i] should be the greatest sigLevel that is no greater than alpha.
            // This is the significance level at which we will perform the test.
            // So, i is the column of the CVTable we will use.
    // Check that it's what we expect:
    if (alpha < DFTest_sigLevel[i] || alpha > DFTest_sigLevel[i+1]) {
      fprintf(stderr, "   WARNING: (dftest.c) couldn't find appropriate significance level for test....\n");
      fprintf(stderr, "   WARNING: (dftest.c) ...sig index may not be set appropriately.\n");
    return(0);
    } else if (DFTest_sigLevel[i]<alpha) { // If DFTest_sigLevel[i] is not exactly equal to alpha, issue a warning.
      printf("\n    WARNING: (in dftest.c) Testing at significance level %lf < alpha = %lf... hope that's okay!", DFTest_sigLevel[i], alpha);
      return(1);
    } else { 
      return(1); 
    }
  }
  // If you made it to this point, something's wrong.
  fprintf(stderr, "   ERROR: (dftest.c) couldn't find appropriate significance level...   Aborting... \n");
  exit(EXIT_FAILURE); 
}



// getSampleSizeIndex -- get the index of the element in DFTest_CVsampleSize array (see dftest.h) 
//                       that is the greatest value not above sampleSize.
int getSampleSizeIndex(const int& sampleSize, int &i){
  int j;
  i = -1;
  assert(sampleSize>9);
  if (sampleSize==DFTest_CVSampleSize[0]){   // sample size is 10
    i=0;
    return(1);
  }
  if (DFTest_CVSampleSize[0]<sampleSize){
    j=0;
    while (DFTest_CVSampleSize[j]<sampleSize){
      j++;
      assert(j<DFTest_CVSampleSizeLength);
    }
    i=j-1;  // DFTest_CVSampleSize[i] should be the greatest sample size that is no greater than sampleSize.
            // So, i is the row of the CVTable we will use.
    // Check that it's what we expect:
    if (sampleSize < DFTest_CVSampleSize[i] || sampleSize > DFTest_CVSampleSize[i+1]) {
      fprintf(stderr, "   WARNING: (dftest.c) couldn't find appropriate sample size for test....\n");
      fprintf(stderr, "   WARNING: (dftest.c) ...sample size index may not be set appropriately.\n");
    return(0);
    } else { 
      return(1); 
    }
  }
  // If you made it to this point, something's wrong.
  fprintf(stderr, "   ERROR: (dftest.c) couldn't find appropriate sample size.  Aborting... \n");
  exit(EXIT_FAILURE); 

}

    
