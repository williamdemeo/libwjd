% Run the regression:
X = load('datafile');
testT = size(X,1); % effective sample size
N = testT+1;       % original sample size (before differencing once)
testY = X(:,3);
X = X(:,1:2);
[Q,R] = qr(X,'0');
coeff = R\(Q'*testY);
numParams = length(coeff);
yHat = X*coeff;
res = testY-yHat;
SSE = res'*res;
dfe = testT-numParams;
MSE = SSE/dfe;
S = R\eye(numParams);
Cov = S*S'*MSE;

needRegOut=1;
% Write results to the regression record:

if needRegOut % Compute all statistics

    yBar = mean(testY);
    regRes = yHat-yBar;
    SSR = regRes'*regRes;
    dfr = numParams-1;
    dft = testT-1;
    diffRes = diff(res);
    
    testReg.num = N;        
    testReg.size = testT;
%    testReg.names = names;
    testReg.coeff = coeff;
    testReg.se = sqrt(diag(Cov));
    testReg.Cov = Cov;
    testReg.tStats.t = coeff./testReg.se;
    testReg.tStats.pVal = 2*(tcdf(-abs(testReg.tStats.t),dfe));    
    testReg.FStat.F = (SSR/dfr)/(SSE/dfe);
    testReg.FStat.pVal = 1-fcdf(testReg.FStat.F,dfr,dfe);
    testReg.yMu = yBar;
    testReg.ySigma = std(testY);
    testReg.yHat = yHat;
    testReg.res = res;
    testReg.DWStat = (diffRes'*diffRes)/SSE;
    testReg.SSR = SSR;
    testReg.SSE = SSE;
    testReg.SST = SSR+SSE;
    testReg.MSE = MSE;
    testReg.RMSE = sqrt(MSE);
    testReg.RSq = 1-SSE/testReg.SST;
    testReg.aRSq = 1-(SSE/testReg.SST)*(dft/dfe);
%    testReg.LL = -normlike([0,sqrt(MSE)],res);
%    testReg.AIC = 2*numParams-2*testReg.LL;
%    testReg.BIC = numParams*log(testT)-2*testReg.LL;
%    testReg.HQC = 2*numParams*log(log(testT))-2*testReg.LL;
                 
else % Compute only the statistics needed to run the test 
    
    %testReg.names = names;
    testReg.coeff = coeff;
    testReg.Cov = Cov;
    
end

