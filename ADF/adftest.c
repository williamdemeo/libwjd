/*
               o When the value is 'AR', the null model

                   y(t) = y(t-1) + b1*(1-L)y(t-1)
                                 + b2*(1-L)y(t-2)
                                 + ... 
                                 + bp*(1-L)y(t-p)
                                 + e(t)                   

                 is tested against the alternative model

                   y(t) = a*y(t-1) + b1*(1-L)y(t-1)
                                   + b2*(1-L)y(t-2)
                                   + ... 
                                   + bp*(1-L)y(t-p)
                                   + e(t)

                 with AR(1) coefficient a < 1.

               o When the value is 't1', a standard t statistic

                   t1 = (a-l)/se

                 is computed from OLS estimates of the AR(1) coefficient a
                 and its standard error se in the alternative model. The
                 test assesses the significance of the restriction a = 1.
%  adf - Structure of regression statistics from the OLS estimation of 
%       coefficients in the alternative model. It has the following fields:
%
%       num         Length of the input series y, with NaNs removed
%       size        Effective sample size, adjusted for lags, difference*
%       names       Regression coefficient names			
%       coeff       Estimated coefficient values
%       se          Estimated coefficient standard errors
%       Cov         Estimated coefficient covariance matrix
%       tStats      t statistics of coefficients and p-values
%       FStat       F statistic and p-value
%       yMu         Mean of y, adjusted for lags, difference*
%       ySigma      Standard deviation of y, adjusted for lags, difference*
%       yHat        Fitted values of y, adjusted for lags, difference*
%       res         Regression residuals
%       DWStat      Durbin-Watson statistic
%       SSR         Regression sum of squares
%       SSE         Error sum of squares
%       SST         Total sum of squares
%       MSE         Mean squared error
%       RMSE        Standard error of the regression
%       RSq         R^2 statistic
%       aRSq        Adjusted R^2 statistic
%       LL          Loglikelihood of data under Gaussian innovations
%       AIC         Akaike information criterion
%       BIC         Bayesian (Schwarz) information criterion
%       HQC         Hannan-Quinn information criterion
%
%       *Lagging and differencing a time series reduce the sample size.
%       Absent any presample values, if y(t) is defined for t = 1:N, then
%       the lagged series y(t-k) is defined for t = k+1:N. Differencing
%       reduces the time base to k+2:N. With p lagged differences, the
%       common time base is p+2:N and the effective sample size is N-(p+1).
% Notes:
%
%   o A suitable value for 'lags' must be determined in order to draw valid
%     inferences from the test. One method is to begin with a maximum lag,
%     such as the one recommended by Schwert [7], and then test down by
%     assessing the significance of the coefficient of the largest lagged
%     change in y, bp. The usual t statistic is appropriate, as reported in
%     the reg output structure. Another method is to combine a measure of
%     fit, such as SSR, with information criteria such as AIC, BIC, and
%     HQC. These statistics are also reported in the reg output structure.
%     Ng and Perron [6] provide further guidelines.
%
%   o The value of 'model' is determined by the growth characteristics of
%     the time series being tested, and should be chosen with a specific
%     testing strategy in mind. As discussed in Elder & Kennedy [4],
%     including too many regressors results in lost power, while including
%     too few biases the test in favor of the null. In general, if a series
%     is growing, the 'TS' model provides a reasonable trend-stationary
%     alternative to a unit-root process with drift. If a series is not
%     growing, 'AR' and 'ARD' models provide reasonable stationary
%     alternatives to a unit-root process without drift. The 'ARD'
%     alternative has mean c/(1-a); the 'AR' alternative has mean 0.
%
%   o Dickey-Fuller statistics follow nonstandard distributions under the
%     null, even asymptotically. Critical values for a range of sample
%     sizes and significance levels have been tabulated using Monte Carlo
%     simulations of the null model with Gaussian innovations and five
%     million replications per sample size. For small samples, values are
%     valid only for Gaussian innovations; for large samples, values are
%     also valid for non-Gaussian innovations. Critical values and p-values
%     are interpolated from the tables. Tables for tests of type 't1' and
%     't2' are identical to those for PPTEST.
%
% Example:
%
%   % Test GDP data for a unit root using a trend-stationary alternative
%   % with 0, 1, and 2 lagged differences:
%
%   load Data_GDP
%   y = log(Data);
%   h = adftest(y,'model','TS','lags',0:2)
%
%   % The test fails to reject the unit-root null with each alternative.

*/

/*
INPUTS
   N        length(y) == total number of observations.
   y        vector of observations (no missing values).
   lags     max number of lagged y's to regress against.
   reject   on return, reject=1 if we can reject the null at alpha sig level (i.e. pVal < alpha)
   pVal     on return, pVal=probability of observing the computed test stat under the null. [1]

NOTES
   [1]      In statistical hypothesis testing, the p-value is the probability of obtaining 
	    a test statistic at least as extreme as the one that was actually observed, assuming 
	    that the null hypothesis is true. One often "rejects the null hypothesis" when the p-value 
	    is less than the significance level α (Greek alpha), which is often 0.05 or 0.01. 
	    When the null hypothesis is rejected, the result is said to be statistically significant. 
*/

void adftest(long *N, double *y, int lags, int *reject, double *pVal);
void adftest(long *N, double *y, int lags, int *reject, double *pVal){ 
  double alpha=0.05;  // significance level for the test. Must be between 0.001 and 0.999.

  double *x, *y, *leadu, *E, *B, *cov, *se, *e, *sigma;
  long i, j, nrow, ncol, mindim;

  T = N-(lags+1);   // Effective sample size 

  // Perform the regression:
  //testReg = runReg(i,y,testT,testLags,testModel,needRegOut);
  mindim = lmin(nrow,ncol); /* mindim is the smaller dimension */
  x = dmalloc(nrow*(ncol+1));
  y = x+(nrow*ncol); /* y is assigned the address of last col of x */

     leadu = dmalloc(mindim);
     E = dmalloc(ncol*ncol);
     B = dmalloc(ncol);
     cov = dmalloc(ncol*ncol);
     se = dmalloc(ncol);
     e = dmalloc(nrow);
     sigma = dmalloc((long)1);
     
     //matlabread(x, nrow, ncol, filename);
     matread(x, nrow, ncol+1, filename); 
     matprint(x, nrow, ncol+1);
     /*matrix is stored contiguously column-wise */

     qrpivot(nrow,ncol,x,E,leadu); /* only send first ncol columns of x */

     printf("\nRough estimate of smallest singular value of X:");
     printf("\nR(%ld,%ld) = %lf",ncol,ncol,x[nrow*(ncol-1)+(ncol-1)]);

     reg(nrow,ncol,x,leadu,E,y,B,cov,se,e,sigma);

     printf("\n\nMSE = %lf\n",*sigma);
     printf("\nCOEFFICIENT \t SE \n");
     for (i = 0; i < ncol; i++)
          printf("%4.5lf \t %4.5lf\n", B[i],se[i]);
     printf("\n\nOBS \t RESIDUAL\n");
     for (i = 0; i < nrow; i++)
          printf("%ld \t %4.5lf\n", i+1,e[i]);
    
  // Set row/column grid for tables of critical values:
  int sampleSizes[15] = {10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 300, 500, 1000, 10000};
  int minT = min(sampSizes);              // Minimum effective sample size
  int maxT = max(max(sampSizes),max(T));  // Maximum effective sample size
  int sampSizes[end] = maxT;              // Force maxT into table

  double *sigLevels, *ptr;
  int i, j, len=0, place=1;
  //sigLevels = [0.001 (0.005:0.005:0.10) (0.125:0.025:0.20) (0.80:0.025:0.875) (0.90:0.005:0.995) 0.999];
  double sigs[4][3] = {{0.005,0.005,0.10}, {0.125,0.025,0.20}, {0.80,0.025,0.875}, {0.90,0.005,0.995}};

  sigLevels = dmalloc(50);
  sigLevels[0] = 0.001;
  ptr = sigLevels+1;
  for (i=0; i<4; ++i){
    len = partition(ptr,sigs[i][0],sigs[i][2],sigs[i][1]);
    assert(len > 1);
    ptr += len;
  }
  sigLevels[49]=0.999;

  double minAlpha = min(sigLevels); // Minimum significance level
  double maxAlpha = max(sigLevels); // Maximum significance level

  // Check if any tests are outside of tables:
  assert(T >= minT);  
  assert(alpha >= minAlpha);  
  assert(alpha <= maxAlpha);
  // Significance levels must be between %5.3f and %5.3f in table of critical vals.',minAlpha,maxAlpha)
  // Get appropriate table of critical values:
  //CVTable = getCVTable(testModel,testType);
    
  // Compute the statistic:
  //[testStat,testPValue] = getStat(i,testT,testLags,testModel,testType,testReg,sigLevels,sampSizes,CVTable,needPValue);
    
  // Test the statistic:
  //[testCValue,testH] = runTest(sigLevels,sampSizes,testT,testType,testAlpha,testStat,CVTable);
    
