function [h,pValue,stat,cValue,reg1,reg2] = egcitest(Y,varargin)
%EGCITEST Engle-Granger cointegration test
%
% Syntax:
%
%   [h,pValue,stat,cValue,reg1,reg2] = egcitest(Y)
%   [h,pValue,stat,cValue,reg1,reg2] = egcitest(Y,param,val,...)
%
% Description:
%
%   Engle-Granger tests assess the null hypothesis of no cointegration
%   among the time series in Y. The test regresses Y(:,1) on Y(:,2:end),
%   then tests the residuals for a unit root.
%
% Input Arguments:
%
%   Y - numObs-by-numDims matrix representing numObs observations of a
%       numDims-dimensional time series y(t), with the last observation the
%       most recent. Observations containing NaN values are removed.
%
% Optional Input Parameter Name/Value Pairs:
%
%   NAME        VALUE
%
%   'creg'      String or cell vector of strings indicating the form of the
%               cointegrating regression, where  y1 = Y(:,1) is regressed
%               on Y2 = Y(:,2:end) and optional deterministic terms in X:
%
%                   y1 = X*a + Y2*b + e
%
%               Values are 'nc' (no constant or trend in X), 'c' (constant
%               but no trend in X), 'ct' (constant and linear trend in X),
%               or 'ctt' (constant, linear trend, and quadratic trend in
%               X). The default value is 'c'.
%
%   'cvec'      Vector or cell vector of vectors containing coefficients
%               [a;b] to be held fixed in the cointegrating regression. The
%               length of a is 0, 1, 2 or 3, depending on 'creg', with
%               coefficient order: constant, linear trend, quadratic trend.
%               The length of b is numDims-1. It is assumed that the
%               coefficient of y1 = Y(:,1) has been normalized to 1. NaN
%               values indicate coefficients to be estimated. If 'cvec' is
%               completely specified (no NaN values), no cointegrating
%               regression is performed. The default value is a completely
%               unspecified cointegrating vector (all NaN values).
%
%   'rreg'      String or cell vector of strings indicating the form of the
%               residual regression. Values are 'ADF', for an augmented
%               Dickey-Fuller test, or 'PP', for a Phillips-Perron test.
%               Test statistics are computed by calling ADFTEST or PPTEST
%               with the 'model' parameter set to 'AR', assuming data have
%               been demeaned or detrended, as necessary, by the choice of
%               'creg'. The default value is 'ADF'.
%
%   'lags'      Scalar or vector of nonnegative integers indicating the
%               number of lags used in the residual regression. The meaning
%               of the parameter depends on the value of 'rreg' (see the
%               documentation for the 'lags' parameter in ADFTEST and
%               PPTEST). The default value is 0.
%
%   'test'      String or cell vector of strings indicating the type of
%               test statistic computed from the residual regression.
%               Values are 't1' (a "tau test") or 't2' (a "z test"). The
%               form of the statistic depends on the value of 'rreg' (see
%               the documentation for the 'test' parameter in ADFTEST and
%               PPTEST). The default value is 't1'.
%
%   'alpha'     Scalar or vector of nominal significance levels for the
%               tests. Values must be between 0.001 and 0.999. The default
%               value is 0.05.
%
%   Single-element parameter values are expanded to the length of any
%   vector value (the number of tests). Vector values must have equal
%   length. If any value is a row vector, all outputs are row vectors.
%
% Output Arguments:
%
%   h - Vector of Boolean decisions for the tests, with length equal to the
%       number of tests. Values of h equal to 1 (true) indicate rejection
%       of the null in favor of the alternative of cointegration. Values of
%       h equal to 0 (false) indicate a failure to reject the null.
%
%   pValue - Vector of p-values of the test statistics, with length equal
%       to the number of tests. p-values are left-tail probabilities.
%
%   stat - Vector of test statistics, with length equal to the number of
%       tests. Statistic depends on the values of 'rreg' and 'test' (see
%       the documentation for ADFTEST and PPTEST).
%
%   cValue - Vector of critical values for the tests, with length equal to
%       the number of tests. Values are for left-tail probabilities. Since
%       residuals are estimated rather than observed, critical values are
%       different from those used in ADFTEST or PPTEST (unless the
%       cointegrating vector is completely specified by 'cvec'). EGCITEST
%       loads tables of critical values from the file Data_EGCITest.mat,
%       then linearly interpolates test values from the tables. Critical
%       values in the tables were computed using methods described in [3].
%
%   reg1 - Structure of regression statistics from the cointegrating
%       regression. 
%
%   reg2 - Structure of regression statistics from the residual regression.
%   
%   The number of records in both reg1 and reg2 is equal to the number of
%   tests. Each record has the following fields:
%
%       num         Length of the regression response y, with NaNs removed
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
%
% Notes:
%
%   o A suitable value for 'lags' must be determined in order to draw valid
%     inferences from the test. See notes on the 'lags' parameter in the
%     documentation for ADFTEST and PPTEST.
%
%   o Samples with less than ~20-40 observations (depending on the
%     dimension of the data) may yield unreliable critical values. See [3].
%
%   o If cointegration is inferred, residuals from the reg1 output can be
%     used as data for the error-correction term in a VEC representation of
%     y(t). Estimation of autoregressive terms can then be performed with
%     VGXVARX, treating the residual series as exogenous. See [1].
%
% Example:
%
%   % Data on the term structure of interest rates in Canada:
% 
%   load Data_Canada
%   Y = Data(:,3:end);
%   names = series(3:end);
%   plot(dates,Y)
%   legend(names,'location','NW')
%   grid on
% 
%   % Test for cointegration (and reproduce row 1 of Table II in [3]):
% 
%   [h,pValue,stat,cValue,reg] = egcitest(Y,'test',{'t1','t2'});
% 
%   % Plot the estimated cointegrating relation y1-Y2*b-X*a:
% 
%   a = reg(2).coeff(1);
%   b = reg(2).coeff(2:3);
%   plot(dates,Y*[1;-b]-a)
%   grid on
%
% References:
%
%   [1] Engle, R. F. and C. W. J. Granger. "Co-Integration and Error-
%       Correction: Representation, Estimation, and Testing." Econometrica.
%       v. 55, 1987, pp. 251-276.
%  
%   [2] Hamilton, J. D. Time Series Analysis. Princeton, NJ: Princeton
%       University Press, 1994.
%
%   [3] MacKinnon, J. G. "Numerical Distribution Functions for Unit Root
%       and Cointegration Tests." Journal of Applied Econometrics. v. 11,
%       1996, pp. 601-618.
%
% See also JCITEST, ADFTEST, PPTEST, VECTOVAR.

% Copyright 2011 The MathWorks, Inc.
% $Revision: 1.1.6.8 $ $Date: 2010/11/08 02:21:33 $

% Parse inputs and set defaults:

parseObj = inputParser;
parseObj.addRequired('Y',@YCheck);
parseObj.addParamValue('creg','c',@cregCheck);
parseObj.addParamValue('cvec',NaN,@cvecCheck);
parseObj.addParamValue('rreg','ADF',@rregCheck);
parseObj.addParamValue('lags',0,@lagsCheck);
parseObj.addParamValue('test','t1',@testCheck);
parseObj.addParamValue('alpha',0.05,@alphaCheck);

parseObj.parse(Y,varargin{:});

Y = parseObj.Results.Y;
creg = lower(parseObj.Results.creg);
cvec = parseObj.Results.cvec;
rreg = upper(parseObj.Results.rreg);
lags = parseObj.Results.lags;
test = lower(parseObj.Results.test);
alpha = parseObj.Results.alpha;

% Check parameter values for commensurate lengths, expand single-element
% values, and convert all variables to columns:

[numTests,rowOutput,creg,cvec,rreg,lags,test,alpha] = sizeCheck(creg,cvec,rreg,lags,test,alpha);

% Remove rows of Y with missing values:

Y(any(isnan(Y),2),:) = [];
[numObs,numDims] = size(Y);

% Set row/column grid for tables of critical values:

sampSizes = [10 15 20 25 30 40 50 75 100 150 200 300 500 1000 10000 Inf];
minT = sampSizes(1);
maxT = sampSizes(end-1); % Report asymptotic values if numObs > maxT

sigLevels = [0.001 (0.005:0.005:0.10) (0.125:0.025:0.20) ...
                   (0.80:0.025:0.875) (0.90:0.005:0.995) 0.999];
minAlpha = sigLevels(1);
maxAlpha = sigLevels(end);

maxDims = 12; % Table rows 1, ..., maxDims

% Check if the sample size is too small:

if numObs < minT
    
	error(message('econ:egcitest:SampleSizeLessThanTabulatedValues', minT))
      
end

% Check if any alpha is outside of the tables:

if any(alpha < minAlpha) || any(alpha > maxAlpha)
    
    error('econ:egcitest:AlphaOutOfRange',...
          'Significance levels must be between the minimum %5.3f \nand the maximum %5.3f in the tables of critical values.',minAlpha,maxAlpha)
      
end

% Check if the number of dimensions is too large:

if numDims > maxDims
    
	error(message('econ:egcitest:numDimsExceedsTabulatedValues', maxDims))
      
end

% Preallocate output variables:

switch nargout
    
    case {0,1}
        
        needPValue = false;
        needCValue = false;
        needReg1Out = false;
        needReg2Out = false;
        
        h = false(numTests,1);
        
    case 2
        
        needPValue = true;
        needCValue = false;
        needReg1Out = false;
        needReg2Out = false;
        
        h = false(numTests,1);
        pValue = NaN(numTests,1);
        
    case 3
        
        needPValue = true;
        needCValue = false;
        needReg1Out = false;
        needReg2Out = false;
        
        h = false(numTests,1);
        pValue = NaN(numTests,1);
        stat = pValue;
        
    case 4
        
        needPValue = true;
        needCValue = true;
        needReg1Out = false;
        needReg2Out = false;
        
        h = false(numTests,1);
        pValue = NaN(numTests,1);
        stat = pValue;
        cValue = pValue;
        
    case 5
        
        needPValue = true;
        needCValue = true;
        needReg1Out = true;
        needReg2Out = false;
        
        h = false(numTests,1);
        pValue = NaN(numTests,1);
        stat = pValue;
        cValue = pValue;
        reg1Fields = {'num','size','names','coeff','se','Cov','tStats',...
                      'FStat','yMu','ySigma','yHat','res','DWStat','SSR',...
                      'SSE','SST','MSE','RMSE','RSq','aRSq','LL','AIC',...
                      'BIC','HQC'};
        reg1 = cell2struct(cell(length(reg1Fields),numTests),reg1Fields,1);
        
    case 6
        
        needPValue = true;
        needCValue = true;
        needReg1Out = true;
        needReg2Out = true;
        
        h = false(numTests,1);
        pValue = NaN(numTests,1);
        stat = pValue;
        cValue = pValue;
        reg1Fields = {'num','size','names','coeff','se','Cov','tStats',...
                      'FStat','yMu','ySigma','yHat','res','DWStat','SSR',...
                      'SSE','SST','MSE','RMSE','RSq','aRSq','LL','AIC',...
                      'BIC','HQC'};
        reg1 = cell2struct(cell(length(reg1Fields),numTests),reg1Fields,1);
        reg2Fields = {'num','size','names','coeff','se','Cov','tStats',...
                      'FStat','yMu','ySigma','yHat','res','autoCov',...
                      'NWEst','DWStat','SSR','SSE','SST','MSE','RMSE',...
                      'RSq','aRSq','LL','AIC','BIC','HQC'};
        reg2 = cell2struct(cell(length(reg2Fields),numTests),reg2Fields,1);

end

% Load critical values:

load Data_EGCITest EGCV

% Initialize loop variables:

lastCreg = [];
lastCvec = [];
lastRreg = [];
lastLags = [];
lastType = [];

% Run the tests:

for testNum = 1:numTests
    
    testCreg = creg{testNum};
    testCvec = cvec{testNum};
    testRreg = rreg{testNum};
    testLags = lags(testNum);
    testType = test{testNum};
    testAlpha = alpha(testNum);
    
    % Get cvec coefficients:
    
    [a,b] = getCvecCoeffs;
       
    % Check to see if statistic is unchanged from last test:
    
    sameStat = isequal(testCreg,lastCreg) && ...
               isequal(testCvec,lastCvec) && ...
               isequal(testRreg,lastRreg) && ...
               isequal(testLags,lastLags) && ...
               isequal(testType,lastType);
    
    if ~sameStat % Recompute the statistic only if it changes
    
        % Perform the cointegrating regression:
    
        testReg1 = runCReg;
        res = testReg1.res;
        
        % Get appropriate table of critical values:
        
        tableDim = 1 + sum(isnan(b)); % 1 + number of estimated 
                                      % stochastic coefficients
        CVTable = getCVTable;
        
        % Compute the statistic:
        
        if needReg2Out
            
            [testStat,testReg2] = runRReg;
            
        else
            
            testStat = runRReg;
            
        end
        
        % Compute the p-value, if requested:
        
        if needPValue
            
            testPValue = getPValue;
            
        end
    
    end
    
    % Test the statistic:
    
    [testCValue,testH] = runTest;
    
    % Add the test results to the outputs:
    
    switch nargout
    
        case {0,1}
        
            h(testNum) = testH;
        
        case 2
        
            h(testNum) = testH;
            pValue(testNum) = testPValue;
        
        case 3
        
            h(testNum) = testH;
            pValue(testNum) = testPValue;
            stat(testNum) = testStat;
        
        case 4
        
            h(testNum) = testH;
            pValue(testNum) = testPValue;
            stat(testNum) = testStat;
            cValue(testNum) = testCValue;
        
        case 5
        
            h(testNum) = testH;
            pValue(testNum) = testPValue;
            stat(testNum) = testStat;
            cValue(testNum) = testCValue;
            reg1(testNum) = testReg1;
        
        case 6
        
            h(testNum) = testH;
            pValue(testNum) = testPValue;
            stat(testNum) = testStat;
            cValue(testNum) = testCValue;
            reg1(testNum) = testReg1;
            reg2(testNum) = testReg2;
            
    end
    
    % Save values to check for changes, next loop:
    
    lastCreg = testCreg;
    lastCvec = testCvec;
    lastRreg = testRreg;
    lastLags = testLags;
    lastType = testType;

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
            reg1 = reg1';
                        
        case 6
            
            h = h';
            pValue = pValue';
            stat = stat';
            cValue = cValue';
            reg1 = reg1';
            reg2 = reg2';
        
    end
    
end

%-------------------------------------------------------------------------
% Check input Y
function OK = YCheck(Y)
            
    if isempty(Y)
        
        error(message('econ:egcitest:DataUnspecified'))
          
    elseif isvector(Y)
        
        error(message('econ:egcitest:DataIsVector'))
          
                    
    elseif ~isnumeric(Y)
        
        error(message('econ:egcitest:DataNonNumeric'))
          
    else
        
        OK = true;
        
    end
    
end % YCHECK
    
%-------------------------------------------------------------------------
% Check value of 'creg' parameter
function OK = cregCheck(creg)
    
    if ~isvector(creg)
        
        error(message('econ:egcitest:CRegNonVector'))
          
    elseif isnumeric(creg) || (iscell(creg) && any(cellfun(@isnumeric,creg)))
        
        error(message('econ:egcitest:CRegNumeric'))
          
    elseif ~all(ismember(lower(creg),{'nc','c','ct','ctt'}))
        
        error(message('econ:egcitest:CRegInvalid'))
          
    else
        
        OK = true;
        
    end
    
end % CREGCHECK

%-------------------------------------------------------------------------
% Check value of 'cvec' parameter
function OK = cvecCheck(cvec)
    
    if ~isvector(cvec) || (iscell(cvec) && any(~cellfun(@isvector,cvec)))
        
        error(message('econ:egcitest:CVecNonVector'))
          
    elseif (~iscell(cvec) && ~isnumeric(cvec)) || (iscell(cvec) && ~all(cellfun(@isnumeric,cvec)))
        
        error(message('econ:egcitest:CVecNonNumeric'))      
          
    else
        
        OK = true;
        
    end
    
end % CVECCHECK
    
%-------------------------------------------------------------------------
% Check value of 'rreg' parameter
function OK = rregCheck(rreg)
    
    if ~isvector(rreg)
        
        error(message('econ:egcitest:RRegNonVector'))
          
    elseif isnumeric(rreg) || (iscell(rreg) && any(cellfun(@isnumeric,rreg)))
        
        error(message('econ:egcitest:RRegNumeric'))      
          
    elseif ~all(ismember(upper(rreg),{'ADF','PP'}))
        
        error(message('econ:egcitest:RRegInvalid'))
          
    else
        
        OK = true;
        
    end
    
end % RREGCHECK
    
%-------------------------------------------------------------------------
% Check value of 'lags' parameter
function OK = lagsCheck(lags)
    
    if ~isvector(lags)
        
        error(message('econ:egcitest:LagsNonVector'))
          
    elseif ~isnumeric(lags)
        
        error(message('econ:egcitest:LagsNonNumeric'))
          
    elseif any(mod(lags,1) ~= 0) || any(lags < 0)
        
        error(message('econ:egcitest:LagsOutOfRange'))
          
    else
        
        OK = true;
        
    end
    
end % LAGSCHECK

%-------------------------------------------------------------------------
% Check value of 'test' parameter
function OK = testCheck(test)
    
    if ~isvector(test)
        
        error(message('econ:egcitest:TypeOfTestNonVector'))
          
    elseif isnumeric(test) || (iscell(test) && any(cellfun(@isnumeric,test)))
        
        error(message('econ:egcitest:TypeOfTestNumeric'))
    
    elseif ~all(ismember(lower(test),{'t1','t2'}))
        
        error(message('econ:egcitest:TypeOfTestInvalid'))
          
    else
        
        OK = true;
        
    end
    
end % TESTCHECK

%-------------------------------------------------------------------------
% Check value of 'alpha' parameter
function OK = alphaCheck(alpha)
    
    if ~isvector(alpha)
        
        error(message('econ:egcitest:AlphaNonVector'))
          
    elseif ~isnumeric(alpha)
        
        error(message('econ:egcitest:AlphaNonNumeric'))
          
    else
        
        OK = true;
        
    end
    
end % ALPHACHECK
    
%-------------------------------------------------------------------------
% Check parameter values for commensurate lengths, expand single-element
% values, and convert all variables to columns
function [numTests,rowOutput,varargout] = sizeCheck(varargin)

% Initialize outputs:

numTests = 1;
rowOutput = false;

% Determine vector lengths, number of tests, row output flag:

for i = 1:nargin
    
    ivar = varargin{i};
    iname = inputname(i);
    
    if (isnumeric(ivar) && ~strcmp(iname,'cvec')) || iscell(ivar)
        paramLength.(iname) = length(ivar);
        if ~isscalar(ivar)
        	rowOutput = rowOutput || (size(ivar,1) == 1);
        end    
    else        
        paramLength.(iname) = 1;   % Single string or cvec
        varargin{i} = varargin(i); % Convert to cell        
    end
    
    numTests = max(numTests,paramLength.(iname));
    
end

% Check for commensurate vector lengths:

for i = 1:(nargin-1)
    iname = inputname(i);
    for j = (i+1):nargin
        jname = inputname(j);
        if (paramLength.(iname) > 1) && (paramLength.(jname) > 1) ...
            && (paramLength.(iname) ~= paramLength.(jname))
        
            error(message('econ:egcitest:ParameterSizeMismatch', iname, jname))
              
        end        
    end
end

% Expand scalars and single strings:

for i = 1:nargin
    
    ivar = varargin{i};
    if paramLength.(inputname(i)) == 1
        varargout{i} = ivar(ones(numTests,1)); % Expand to column output
    else
        varargout{i} = ivar(:);  % Column output
    end
    
end

end % SIZECHECK

%-------------------------------------------------------------------------
% Get cvec coefficients
function [a,b] = getCvecCoeffs
        
% Allow any all-NaN cvec to proceed to the estimation step;
% check other user-specified cvecs for proper dimension:
    
switch testCreg

    case 'nc'

        if all(isnan(testCvec))

           testCvec = NaN(numDims-1,1); % Ensure proper dimension

        end

        if length(testCvec) ~= (numDims-1)

            error(message('econ:egcitest:NCCvecWrongSize'))

        end

        a = []; % User-specified deterministic coefficients
        b = testCvec; % User-specified stochastic coefficients

    case 'c'

        if all(isnan(testCvec))

           testCvec = NaN(numDims,1); % Ensure proper dimension

        end

        if length(testCvec) ~= numDims

            error(message('econ:egcitest:CCvecWrongSize'))

        end

        a = testCvec(1); % User-specified deterministic coefficients
        b = testCvec(2:end); % User-specified stochastic coefficients

    case 'ct'

        if all(isnan(testCvec))

           testCvec = NaN(numDims+1,1); % Ensure proper dimension

        end

        if length(testCvec) ~= (numDims+1)

            error(message('econ:egcitest:CTCvecWrongSize'))

        end

        a = testCvec(1:2); % User-specified deterministic coefficients
        b = testCvec(3:end); % User-specified stochastic coefficients
        
	case 'ctt'

        if all(isnan(testCvec))

           testCvec = NaN(numDims+2,1); % Ensure proper dimension

        end

        if length(testCvec) ~= (numDims+2)

            error(message('econ:egcitest:CTTCvecWrongSize'))

        end

        a = testCvec(1:3); % User-specified deterministic coefficients
        b = testCvec(4:end); % User-specified stochastic coefficients

end
    
end % GETCVECCOEFFS

%-------------------------------------------------------------------------
% Perform the cointegrating regression
function testReg1 = runCReg

% Set up the regression:

y1 = Y(:,1);
Y2 = Y(:,2:end);
a = a(:); % Convert to column
b = b(:); % Convert to column
coeff = [a;b];
specified = ~isnan(coeff); % Indices of user-specified coefficients

switch testCreg
    
    case 'nc' % Regression: y1 = Y2*b + e
        
        % Design matrix:
        
        X = Y2;
        
        % Coefficient names:
        
        names = strcat({'b'},num2str((1:numDims-1)','%-d'));
        
    case 'c' % Regression: y1 = a + Y2*b + e
        
        % Design matrix:
        
        X = [ones(numObs,1),Y2];
        
        % Coefficient names:
        
        names = cat(1,'a',strcat({'b'},num2str((1:numDims-1)','%-d')));
        
    case 'ct' % Regression: y1 = a1 + a2*t + Y2*b + e
        
        % Design matrix:
        
        X = [ones(numObs,1),(1:numObs)',Y2];
        
        % Coefficient names:
        
        names = cat(1,'a1','a2',strcat({'b'},num2str((1:numDims-1)','%-d')));
        
	case 'ctt' % Regression: y1 = a1 + a2*t + a3*t^2 + Y2*b + e
        
        % Design matrix:
        
        X = [ones(numObs,1),(1:numObs)',((1:numObs)').^2,Y2];
        
        % Coefficient names:
        
        names = cat(1,'a1','a2','a3',strcat({'b'},num2str((1:numDims-1)','%-d')));
        
end

% Adjust the regression to accommodate user-specified coefficients:
        
y1_adjusted = y1-X(:,specified)*coeff(specified,1);
X_adjusted = X(:,~specified);

if size(X_adjusted,1) < size(X_adjusted,2)
    
    error(message('econ:egcitest:RankDeficientDesignMatrix'));
      
end

% Run the regression:

[Q,R] = qr(X_adjusted,0);
estimated_coeff = R\(Q'*y1_adjusted);
coeff(~specified) = estimated_coeff;
numEstParams = length(estimated_coeff);
numParams = length(coeff);
yHat = X*coeff;
res = y1-yHat;
SSE = res'*res;
dfe = numObs-numEstParams;
MSE = SSE/dfe;
S = R\eye(numEstParams);
EstCov = S*S'*MSE;
Cov = zeros(numParams);
Cov(~specified,~specified) = EstCov;

% Write results to the regression record:

if needReg1Out % Compute all statistics

    yBar = mean(y1);
    regRes = yHat-yBar;
    SSR = regRes'*regRes;
    dfr = numEstParams-1;
    dft = numObs-1;
    diffRes = diff(res);
    
    testReg1.num = numObs;        
    testReg1.size = numObs;
    testReg1.names = names;
    testReg1.coeff = coeff;
    testReg1.se = sqrt(diag(Cov));
    testReg1.Cov = Cov;
    testReg1.tStats.t = coeff./testReg1.se;
    testReg1.tStats.pVal = 2*(tcdf(-abs(testReg1.tStats.t),dfe));    
    testReg1.FStat.F = (SSR/dfr)/(SSE/dfe);
    testReg1.FStat.pVal = 1-fcdf(testReg1.FStat.F,dfr,dfe);
    testReg1.yMu = yBar;
    testReg1.ySigma = std(y1);
    testReg1.yHat = yHat;
    testReg1.res = res;
    testReg1.DWStat = (diffRes'*diffRes)/SSE;
    testReg1.SSR = SSR;
    testReg1.SSE = SSE;
    testReg1.SST = SSR+SSE;
    testReg1.MSE = MSE;
    testReg1.RMSE = sqrt(MSE);
    testReg1.RSq = 1-SSE/testReg1.SST;
    testReg1.aRSq = 1-(SSE/testReg1.SST)*(dft/dfe);
    testReg1.LL = -normlike([0,sqrt(MSE)],res);
    testReg1.AIC = 2*numEstParams-2*testReg1.LL;
    testReg1.BIC = numEstParams*log(numObs)-2*testReg1.LL;
    testReg1.HQC = 2*numEstParams*log(log(numObs))-2*testReg1.LL;
                 
else % Compute only the statistics needed to run the test 
    
    testReg1.res = res;
    
end

end % RUNCREG

%-------------------------------------------------------------------------
% Get table of critical values
function CVTable = getCVTable

switch testCreg

    case 'nc'

        switch testType

            case 't1'

                CVTable = EGCV(:,:,1,1,tableDim);

            case 't2'

                CVTable = EGCV(:,:,1,2,tableDim);

        end

    case 'c'

        switch testType

            case 't1'

                CVTable = EGCV(:,:,2,1,tableDim);

            case 't2'

                CVTable = EGCV(:,:,2,2,tableDim);

        end

    case 'ct'

        switch testType

            case 't1'

                CVTable = EGCV(:,:,3,1,tableDim);

            case 't2'

                CVTable = EGCV(:,:,3,2,tableDim);

        end

    case 'ctt'

        switch testType

            case 't1'

                CVTable = EGCV(:,:,4,1,tableDim);

            case 't2'

                CVTable = EGCV(:,:,4,2,tableDim);

        end

end

end % GETCVTABLE
        
%-------------------------------------------------------------------------
% Perform the residual regression
function [testStat,varargout] = runRReg
        
switch testRreg
    
    case 'ADF'
        
        s = warning('off','all'); % Disable p-value warnings from ADFTEST
        warning('on','econ:adftest:InvalidStatistic') % Enable invalid ADFTEST warning
        
        if needReg2Out
            
            [~,~,testStat,~,testReg2] = adftest(res,'lags',testLags,'test',testType);
                        
            % Make structure similar to PP structure for concatenation:
            testReg2.autoCov = 'PP tests only';
            testReg2.NWEst = 'PP tests only';
            testReg2 = orderfields(testReg2,reg2);
            
            varargout{1} = testReg2;
                        
        else
            
            [~,~,testStat] = adftest(res,'lags',testLags,'test',testType);
        
        end
        
        warning(s)
         
    case 'PP'
        
        s = warning('off','all'); % Disable p-value warnings from PPTEST
        
       	if needReg2Out
            
            [~,~,testStat,~,testReg2] = pptest(res,'lags',testLags,'test',testType);
                        
            varargout{1} = testReg2;
                        
        else
            
            [~,~,testStat] = pptest(res,'lags',testLags,'test',testType);
                    
        end
        
        warning(s)

end

end % RUNRREG

%-------------------------------------------------------------------------
% Get p-values
function testPValue = getPValue
    
% P-values are estimated using two successive 1D interpolations:
%
% 1. Find all critical values associated with the sample size.
% 2. Find the cumulative probability associated with the test statistic.

CVTableRow = interp2(sigLevels,sampSizes,CVTable,sigLevels,numObs,'linear');

% Left-tailed test
    
    if testStat <= CVTableRow(1)
        
        testPValue = sigLevels(1);
        warning('econ:egcitest:LeftTailStatTooSmall',...
                'Test statistic #%d below tabulated critical values: \nminimum p-value = %5.3f reported.',testNum,testPValue)
            
    elseif testStat >= CVTableRow(end)
        
        testPValue = sigLevels(end);
        warning('econ:egcitest:LeftTailStatTooBig',...
                'Test statistic #%d above tabulated critical values: \nmaximum p-value = %5.3f reported.',testNum,testPValue)
            
    else
        
        testPValue = interp1(CVTableRow,sigLevels,testStat,'linear');
            
    end

end % GETPVALUE

%-------------------------------------------------------------------------
% Test the statistic
function [testCValue,testH] = runTest
    
testCValue = interp2(sigLevels,sampSizes,CVTable,testAlpha,numObs,'linear');

if needCValue && (numObs > maxT)
    
	asymptoticRow = CVTable(end,:);
	asymptoticCV = interp1(sigLevels,asymptoticRow,testAlpha,'linear');
    warning('econ:egcitest:ReportAsymptoticCV',...
            'Sample size of the data \nis more than the maximum size %d \nin the table of critical values. \nUsing critical value %.4f at \nmaximum size for the test. Compare \nasymptotic critical value %.4f.',maxT,testCValue,asymptoticCV)
    
end

testH = (testStat < testCValue); % Left-tailed test

end % RUNTEST

%-------------------------------------------------------------------------

end % EGCITEST