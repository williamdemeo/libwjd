function [h,pValue,stat,cValue,reg] = kpsstest(y,varargin)
%KPSSTEST KPSS test for stationarity
%
% Syntax:
%
%   [h,pValue,stat,cValue,reg] = kpsstest(y)
%   [h,pValue,stat,cValue,reg] = kpsstest(y,param1,val1,param2,val2,...)
%
% Description:
%
%   The test of Kwiatkowski, Phillips, Schmidt and Shin (KPSS) assesses the
%   null hypothesis that a univariate time series y is trend stationary
%   against the alternative that it is a nonstationary unit-root process.
%   The test uses the structural model
%
%       y(t) = c(t) + d*t + u1(t),
%
%       c(t) = c(t-1) + u2(t),
%
%   where u1(t) is a stationary process and u2(t) ~ iid(0,sigma^2). The
%   null hypothesis is that sigma^2 = 0, so that the random walk term c(t)
%   becomes a constant intercept. The alternative is that sigma^2 > 0,
%   which introduces the unit root in the random walk.
%
% Input Arguments:
%
%   y - Vector of time-series data. The last element is the most recent
%       observation. NaNs indicating missing values are removed.
%
% Optional Input Parameter Name/Value Pairs:
%
%   NAME        VALUE
%
%   'lags'      Scalar or vector of nonnegative integers indicating the
%               number of autocovariance lags to include in the Newey-West
%               estimator of the long-run variance. The default value is 0.
%
%   'trend'     Scalar or vector of Boolean values indicating whether or
%               not to include the deterministic trend term d*t in the
%               model. The default value is true.
%
%   'alpha'     Scalar or vector of nominal significance levels for the
%               tests. Values must be between 0.01 and 0.10. The default
%               value is 0.05.
%
%   Scalar parameter values are expanded to the length of any vector value
%   (the number of tests). Vector values must have equal length. If any
%   value is a row vector, all outputs are row vectors.
%
% Output Arguments:
%
%   h - Vector of Boolean decisions for the tests, with length equal to the
%       number of tests. Values of h equal to 1 indicate rejection of the
%       trend-stationary null in favor of the unit-root alternative. Values
%       of h equal to 0 indicate a failure to reject the null.
%
%   pValue - Vector of p-values of the test statistics, with length equal
%       to the number of tests. Values are right-tail probabilities.
%
%   stat - Vector of test statistics, with length equal to the number of
%       tests. Statistics are computed using an OLS regression of y on an
%       intercept and, if 'trend' is true, a trend. Residuals e from the
%       regression are used to form the test statistic
%
%           stat = sum(S(t)^2)/((s^2)*(T^2))
%
%       where S(t) = e1 + ... + et, s^2 is the Newey-West estimator of the
%       long-run variance, and T is the sample size.
%
%   cValue - Vector of critical values for the tests, with length equal to
%       the number of tests. Values are for right-tail probabilities.
%
%   reg - Structure of statistics from the OLS regression. The number of
%       records is equal to the number of tests. Each record has the
%       following fields:
%
%       num         Length of the input series y, with NaNs removed
%       size        Effective sample size (same as num)
%       names       Regression coefficient names			
%       coeff       Estimated coefficient values
%       se          Estimated coefficient standard errors
%       Cov         Estimated coefficient covariance matrix
%       tStats      t statistics of coefficients and p-values
%       FStat       F statistic and p-value
%       yMu         Mean of y
%       ySigma      Standard deviation of y
%       yHat        Fitted values of y
%       res         Regression residuals
%       autoCov     Estimated residual autocovariances
%       NWEst       Newey-West estimator
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
% Notes:
%
%   o A suitable value for 'lags' must be determined in order to draw valid
%     inferences from the test. One method is to begin with a small number
%     of lags and then evaluate the sensitivity of the results by adding
%     more lags. For consistency of the Newey-West estimator, the number of
%     lags must go to infinity as the sample size increases. [2] suggests
%     that a number of lags on the order of sqrt(T), where T is the sample
%     size, is often satisfactory under both the null and the alternative.
%
%   o The value of 'trend' is determined by the growth characteristics of
%     the time series being tested, and should be chosen with a specific
%     testing strategy in mind. If a series is growing, including a trend
%     term (setting 'trend' to true) provides a reasonable comparison of a
%     trend-stationary null and a unit-root process with drift. If a series
%     does not exhibit long-term growth characteristics, including a trend
%     term is inappropriate.
%
%   o Test statistics follow nonstandard distributions under the null, even
%     asymptotically. Asymptotic critical values for a standard set of
%     significance levels between 0.01 and 0.1, for models with and without
%     a trend, have been tabulated in [2] using Monte Carlo simulations.
%     Critical values and p-values reported by KPSSTEST are interpolated
%     from the tables.
%
% Example:
%
%   % Reproduce the first row of the second half of Table 5 in [2]:
%
%   load Data_NelsonPlosser
%   y = log(Dataset.GNPR);
%   [~,~,stat] = kpsstest(y,'lags',0:8,'trend',true)
%
% References:
%
%   [1] Hamilton, J. D. Time Series Analysis. Princeton, NJ: Princeton
%       University Press, 1994.
%
%   [2] Kwiatkowski, D., P. C. B. Phillips, P. Schmidt and Y. Shin.
%       "Testing the Null Hypothesis of Stationarity against the
%       Alternative of a Unit Root." Journal of Econometrics. Vol. 54,
%       1992, pp. 159-178.
%
%   [3] Newey, W. K., and K. D. West. "A Simple Positive Semidefinite,
%       Heteroskedasticity and Autocorrelation Consistent Covariance
%       Matrix." Econometrica. Vol. 55, 1987, pp. 703-708.
%  
% See also LMCTEST, PPTEST, ADFTEST, VRATIOTEST. 

% Copyright 2009-2010 The MathWorks, Inc.
% $Revision: 1.1.6.8 $ $Date: 2010/10/08 16:41:20 $

% Parse inputs and set defaults:

parseObj = inputParser;
parseObj.addRequired('y',@yCheck);
parseObj.addParamValue('lags',0,@lagsCheck);
parseObj.addParamValue('trend',true,@trendCheck);
parseObj.addParamValue('alpha',0.05,@alphaCheck);

parseObj.parse(y,varargin{:});

y = parseObj.Results.y;
lags = parseObj.Results.lags;
trend = parseObj.Results.trend;
alpha = parseObj.Results.alpha;

% Check parameter values for commensurate lengths, expand scalars and
% single strings, and convert all variables to columns:

[numTests,rowOutput,lags,trend,alpha] = sizeCheck(lags,trend,alpha);
y = y(:);

% Adjust y for missing values:

y(isnan(y)) = []; % Remove missing values
T = length(y);    % Effective sample size

if T < 1
    
    error(message('econ:kpsstest:EffectiveSampleSizeLessThanOne'))
      
end

if any(trend) && (T < 2)
    
	error(message('econ:kpsstest:OnePointTrend'))
      
end

% Preallocate output variables:

switch nargout
    
    case {0,1}
        
        needPValue = false;
        needRegOut = false;
        
        h = false(numTests,1);
        
    case 2
        
        needPValue = true;
        needRegOut = false;
        
        h = false(numTests,1);
        pValue = NaN(numTests,1);
        
    case 3
        
        needPValue = true;
        needRegOut = false;
        
        h = false(numTests,1);
        pValue = NaN(numTests,1);
        stat = pValue;
        
    case 4
        
        needPValue = true;
        needRegOut = false;
        
        h = false(numTests,1);
        pValue = NaN(numTests,1);
        stat = pValue;
        cValue = pValue;
        
    case 5
        
        needPValue = true;
        needRegOut = true;
        
        h = false(numTests,1);
        pValue = NaN(numTests,1);
        stat = pValue;
        cValue = pValue;
        regFields = {'num','size','names','coeff','se','Cov','tStats',...
                     'FStat','yMu','ySigma','yHat','res','autoCov',...
                     'NWEst','DWStat','SSR','SSE','SST','MSE','RMSE',...
                     'RSq','aRSq','LL','AIC','BIC','HQC'};
        reg = cell2struct(cell(length(regFields),numTests),regFields,1);

end

% Initialize loop variables:

lastLags = [];
lastTrend = [];

% Run the tests:

for i = 1:numTests
    
    testLags = lags(i);
    testTrend = trend(i);
    testAlpha = alpha(i);
    
    % Check to see if statistic is unchanged from last test:
    
    sameStat = isequal(testLags,lastLags) && isequal(testTrend,lastTrend);
    
    if ~sameStat % Recompute the statistic only if it changes
    
        % OLS regression of y on p + q*t:
    
        testReg = runReg(y,testLags,testTrend,needRegOut);
    
        % Get appropriate table of critical values from [2]:

        sigLevels = [0.010 0.025 0.050 0.100];

        if testTrend
        
            % Alpha    0.010  0.025  0.050  0.100    
            % -------------------------------------
            CVTable = [0.216  0.176  0.146  0.119];       
        
        else
        
            % Alpha    0.010  0.025  0.050  0.100    
            % -------------------------------------
            CVTable = [0.739  0.574  0.463  0.347]; 
        
        end
    
        % Compute the statistic:
    
        [testStat,testPValue] = getStat(i,T,testReg,sigLevels,CVTable,needPValue);
    
    end
    
    % Test the statistic:
    
    [testCValue,testH] = runTest(sigLevels,testAlpha,testStat,CVTable);
    
    % Add the test results to the outputs:
    
    switch nargout
    
        case {0,1}
        
            h(i) = testH;
        
        case 2
        
            h(i) = testH;
            pValue(i) = testPValue;
        
        case 3
        
            h(i) = testH;
            pValue(i) = testPValue;
            stat(i) = testStat;
        
        case 4
        
            h(i) = testH;
            pValue(i) = testPValue;
            stat(i) = testStat;
            cValue(i) = testCValue;
        
        case 5
        
            h(i) = testH;
            pValue(i) = testPValue;
            stat(i) = testStat;
            cValue(i) = testCValue;
            reg(i) = testReg;

    end
    
    % Save values to check for changes, next loop:
    
    lastLags = testLags;
    lastTrend = testTrend;

end

% Display outputs as row vectors if any parameter value is a row vector:

if rowOutput
    
    switch nargout
        
        case {0,1}

            h = h';

        case 2
            
            h = h';
            pValue = pValue';
            
        case 3
            
            h = h';
            pValue = pValue';
            stat = stat';
            
        case 4
            
            h = h';
            pValue = pValue';
            stat = stat';
            cValue = cValue';
            
        case 5
            
            h = h';
            pValue = pValue';
            stat = stat';
            cValue = cValue';
            reg = reg';
        
    end
    
end

%-------------------------------------------------------------------------
% Check input y
function OK = yCheck(y)
            
    if isempty(y)
        
        error(message('econ:kpsstest:DataUnspecified'))
          
    elseif ~isnumeric(y)
        
        error(message('econ:kpsstest:DataNonNumeric'))
          
    elseif ~isvector(y)
        
        error(message('econ:kpsstest:DataNonVector'))
          
    else
        
        OK = true;
        
    end

%-------------------------------------------------------------------------
% Check value of 'lags' parameter
function OK = lagsCheck(lags)
    
    if ~isnumeric(lags)
        
        error(message('econ:kpsstest:LagsNonNumeric'))
          
    elseif ~isvector(lags)
        
        error(message('econ:kpsstest:LagsNonVector'))
          
    elseif any(mod(lags,1) ~= 0) || any(lags < 0)
        
        error(message('econ:kpsstest:LagsNonNonnegativeInteger'))
          
    else
        
        OK = true;
        
    end

%-------------------------------------------------------------------------
% Check value of 'trend' parameter
function OK = trendCheck(trend)
    
    if ~islogical(trend)
        
        error(message('econ:kpsstest:TrendNonBoolean'))
          
    elseif ~isvector(trend)
        
        error(message('econ:kpsstest:TrendNonVector'))
          
    else
        
        OK = true;
        
    end

%-------------------------------------------------------------------------
% Check value of 'alpha' parameter
function OK = alphaCheck(alpha)

minAlpha = 0.01;
maxAlpha = 0.1;
    
    if ~isnumeric(alpha)
        
        error(message('econ:kpsstest:AlphaNonNumeric'))
          
    elseif ~isvector(alpha)
        
        error(message('econ:kpsstest:AlphaNonVector'))
          
    elseif any(alpha < minAlpha) || any(alpha > maxAlpha)
        
        error('econ:kpsstest:AlphaOutOfRange',...
              'Significance levels must be between the minimum %5.3f \nand the maximum %5.3f in the table of critical values.',minAlpha,maxAlpha)
          
    else
        
        OK = true;
        
    end
 
%-------------------------------------------------------------------------
% Check parameter values for commensurate lengths, expand scalars, and
% convert all variables to columns
function [numTests,rowOutput,varargout] = sizeCheck(varargin)

% Initialize outputs:

numTests = 1;
rowOutput = false;

% Determine vector lengths, number of tests, row output flag:

for i = 1:nargin
        
    ivar = varargin{i};
    iname = inputname(i);
    
    paramLength.(iname) = length(ivar);
    numTests = max(numTests,paramLength.(iname));
    
    if ~isscalar(ivar)
        rowOutput = rowOutput || (size(ivar,1) == 1);
    end    
    
end

% Check for commensurate vector lengths:

for i = 1:(nargin-1)
    iname = inputname(i);
    for j = (i+1):nargin
        jname = inputname(j);
        if (paramLength.(iname) > 1) && (paramLength.(jname) > 1) ...
            && (paramLength.(iname) ~= paramLength.(jname))
        
            error(message('econ:kpsstest:ParameterSizeMismatch', iname, jname))
              
        end        
    end
end

% Expand scalars:

for i = 1:nargin
    
    ivar = varargin{i};
    if paramLength.(inputname(i)) == 1
        varargout{i} = ivar(ones(numTests,1));
    else
        varargout{i} = ivar(:);  % Column output
    end
    
end

%-------------------------------------------------------------------------
% Perform the test regression
function testReg = runReg(y,testLags,testTrend,needRegOut)
                                        
T = length(y);

% Set up the regression:

if testTrend
    
    X = [ones(T,1),(1:T)'];
    names = {'p';'q'};
    
else
    
    X = ones(T,1);
    names = {'p'};
    
end

% Run the regression:

[Q,R] = qr(X,0);
coeff = R\(Q'*y);
yHat = X*coeff;
res = y-yHat;

% Compute the Newey-West estimator. Follows Hamilton [3] p. 514.

	% Estimated residual autocovariances:
    
	gamma = zeros(testLags+1,1);
	for j = 0:testLags
        gamma(j+1) = (1/T)*res(j+1:end)'*res(1:end-j);
	end

	% Newey-West estimator:
    
	S = 0;
	for j = 1:testLags
        S = S+(1-j/(testLags+1))*gamma(j+1);
	end
	lambdaSq = gamma(1)+2*S;

% Write results to the regression record:

if needRegOut % Compute all statistics

    yBar = mean(y);
    regRes = yHat-yBar;
    SSR = regRes'*regRes;
    numParams = length(coeff);
    dfr = numParams-1;
    dft = T-1;
    diffRes = diff(res);
    SSE = res'*res;
    dfe = T-numParams;
    MSE = SSE/dfe;
    S = R\eye(numParams);
    Cov = S*S'*MSE;
 
    testReg.num = T;        
    testReg.size = T;
    testReg.names = names;
    testReg.coeff = coeff;
    testReg.se = sqrt(diag(Cov));
    testReg.Cov = Cov;
    testReg.tStats.t = coeff./testReg.se;
    testReg.tStats.pVal = 2*(tcdf(-abs(testReg.tStats.t),dfe));    
    testReg.FStat.F = (SSR/dfr)/(SSE/dfe);
    testReg.FStat.pVal = 1-fcdf(testReg.FStat.F,dfr,dfe);
    testReg.yMu = yBar;
    testReg.ySigma = std(y);
    testReg.yHat = yHat;
    testReg.res = res;
    testReg.autoCov = gamma;
    testReg.NWEst = lambdaSq;
    testReg.DWStat = (diffRes'*diffRes)/SSE;
    testReg.SSR = SSR;
    testReg.SSE = SSE;
    testReg.SST = SSR+SSE;
    testReg.MSE = MSE;
    testReg.RMSE = sqrt(MSE);
    testReg.RSq = 1-SSE/testReg.SST;
    testReg.aRSq = 1-(SSE/testReg.SST)*(dft/dfe);
    testReg.LL = -normlike([0,sqrt(MSE)],res);
    testReg.AIC = 2*numParams-2*testReg.LL;
    testReg.BIC = numParams*log(T)-2*testReg.LL;
    testReg.HQC = 2*numParams*log(log(T))-2*testReg.LL;
                 
else % Compute only the statistics needed to run the test 
    
    testReg.res = res;
    testReg.NWEst = lambdaSq;
    
end

%-------------------------------------------------------------------------
% Compute the statistic:
function [testStat,testPValue] = getStat(i,T,testReg,sigLevels,CVTable,needPValue)

eSum = cumsum(testReg.res);
testStat = (eSum'*eSum)/(testReg.NWEst*T^2);

if needPValue

    if testStat < CVTable(end)
    
        testPValue = sigLevels(end);
        warning('econ:kpsstest:StatTooSmall',...
                'Test statistic #%d below tabulated critical values: \nmaximum p-value = %5.3f reported.',i,testPValue)
    
    elseif testStat > CVTable(1)
    
        testPValue = sigLevels(1);
        warning('econ:kpsstest:StatTooBig',...
                'Test statistic #%d above tabulated critical values: \nminimum p-value = %5.3f reported.',i,testPValue)
    
    else
    
        testPValue = interp1(CVTable,sigLevels,testStat,'linear');
    
    end
    
else
    
    testPValue = NaN;

end

%-------------------------------------------------------------------------
% Test the statistic:
function [testCValue,testH] = runTest(sigLevels,testAlpha,testStat,CVTable)

testCValue = interp1(sigLevels,CVTable,testAlpha,'linear');
testH = (testStat > testCValue); % Right-tailed test