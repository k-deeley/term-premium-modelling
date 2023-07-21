function decomposition = estimateACM( yields, maturities, numFactors )

arguments
    % Matrix of yields (%).
    yields(:, :) double    
    % Vector of maturities (months).
    maturities(1, :) double
    % Number of factors (principal components) to use in the model.
    numFactors(1, 1) double = 3
end % arguments

%% Prepare the data.
% Convert the yields from percentages to proportions.
percent2proportion = 1 / 100;
yields = yields * percent2proportion;

% Short-term interest rate. This acts as the instantaneous risk-free rate.
% We estimate this using the yield of the first month. In turn, we estimate
% the yield of the first month by taking the first maturity, and 
% normalizing it to convert to a monthly rate.
monthsPerYear = 12;
shortTermInterestRate = yields(:, 1) / (maturities(1) * monthsPerYear);

% Compute the principal components of the yields. The principal component 
% factors are interpreted as follows.
%
% First factor:
%
% * The expected returns of all maturities move together, so a single
% variable x_t (first factor) can describe the time-variation of the
% expected returns.
%
% The remaining factors are, respectively:
%
% * Level
% * Slope
% * Curvature

[yieldsNorm, yieldsMean, yieldsStd] = zscore( yields );
[~, factors] = pca( yieldsNorm, "Centered", false );

% Extract the first "numComponents" principal components.
selectedFactors = factors(:, 1:numFactors);

% Write down the number of observations (times) and the number of
% maturities.
[numTimes, numMaturities] = size( yields );  

%% Step 1: Estimate the P dynamics using a VAR.

% Create and estimate a VAR model for the yield principal components
% (qquation (1) in ACM).
numLags = 1;
VARModel = varm( numFactors, numLags );
VARModel.SeriesNames = "PCA_Factor_" + (1:numFactors);
VARModel.Description = "Vector autoregressive zero-mean " + ...
    "model of order 1 for the yield principal components";

% Set the model constants to be zero.
mu = 0;
VARModel.Constant = mu * ones( numFactors, 1 );

% Estimate the model.
[VARModel, ~, ~, modelResiduals] = estimate( VARModel, selectedFactors );
Phi = VARModel.AR{1};
mu = VARModel.Constant;

% Compute the covariance matrix of the residuals.
Sigma = VARModel.Covariance;

%% Step 2: Estimate excess returns via linear regression.
% Convert yields into log prices. We use the formula:
%
%   y_t^n = -1/n * p_t^n (1)
% 
% where: 
%   * y_t^n is the continuously compounded yield ("log yield")
%   * p_t^n is the log price of n-year discount bond at time t
%   * n is the maturity (Cochrane and Piazzesi, 2005)
% No need to divide by 12; our maturities are already in terms of years
% (unlike the original code).
prices = (-1) * maturities .* yields;

% Compute excess returns (equation (6) in ACM). If we are getting the
% continously compounded yield, the log is not necessary.
% https://quant.stackexchange.com/questions/28426/how-to-derive-the-relationship-between-log-yield-and-log-price
% Here, the rows represent continuous time (monthly data) and the columns 
% represent discrete time (bond maturity periods).
% 
% t - row index (continuous, monthly)
% n - column index (discrete, bond maturity period)
excessReturns = prices(2:end, 1:end-1) - prices(1:end-1, 2:end) - ...
    shortTermInterestRate(1:end-1);
numExcessReturnSeries = width( excessReturns );

% OLS estimator (equation (15) in ACM).
Z = [ones( numTimes-1, 1 ), modelResiduals, selectedFactors(1:end-1, :)];
allCoeffs = Z \ excessReturns;
aCoeffs = allCoeffs(1, :).'; % numMaturities x 1
bCoeffs = allCoeffs(2:numFactors+1, :).'; % numMaturities x numFactors
cCoeffs = allCoeffs(numFactors+2:end, :).'; % numMaturities x numFactors
returnPricingErrors = excessReturns - Z * allCoeffs;
sigma2 = trace( returnPricingErrors*returnPricingErrors' ) / ...
    (numMaturities * numTimes);

%% Step 3: Compute the market prices of risk.
% Construct the B* matrix (below equation (13) in ACM).
Bstar = zeros( numExcessReturnSeries, numFactors^2 );
for k = 1 : numExcessReturnSeries
    currentRow = bCoeffs(k, :);    
    Bstar(k, :) = reshape( currentRow.' * currentRow, 1, [] );        
end % for

% Compute the market prices of risk (the parameters lambda0 and lambda1).
lambda0 = bCoeffs \ (aCoeffs + 0.5 * (Bstar * Sigma(:) + ...
    sigma2 * ones( numExcessReturnSeries, 1 )));
lambda1 = bCoeffs \ cCoeffs;

%% Step 4: Short rate regression.

% Regress the short-term interest rates on a constant and the selected 
% factors. The model is as follows:
% shortTermInterestRate = delta0 + selectedFactors * delta1
deltas = [ones( numTimes, 1 ), selectedFactors] \ shortTermInterestRate;
delta0 = deltas(1);
delta1 = deltas(2:end);

%% Step 5: Bond pricing.
sqrtSigma = sqrtm( Sigma );

% Construct risk-neutral pricing coefficients (equations (25) and (26) in
% ACM).
q0 = mu - lambda0;
qx = Phi - lambda1;

% Preallocate.
aFull = NaN( 1, numMaturities );
bFull = NaN( numFactors, numMaturities );

% Bond prices under P with sigma = 0
aFull_P = NaN( 1, numMaturities );
bFull_P = NaN( numFactors, numMaturities );

% Bond prices under Q with sigma = 0
aFull_Q = NaN( 1, numMaturities );
bFull_Q = NaN( numFactors, numMaturities );

for n = 1 : numMaturities    
    if n == 1        
        aFull(1, n) = -delta0;
        bFull(:, n) = -delta1;
        aFull_P(1,n) = -delta0;
        bFull_P(:,n) = -delta1;
        aFull_Q(1,n) = -delta0;
        bFull_Q(:,n) = -delta1;        
    else        
        aFull(1, n) = aFull(1, n-1) + dot( bFull(:, n-1), q0 ) + ...
            0.5 * (bFull(:, n-1).' * Sigma * bFull(:, n-1) + sigma2) - delta0;
        bFull(:, n) = -delta1 + qx.' * bFull(:, n-1);
        aFull_P(1, n) = aFull_P(1, n-1) + mu.' * bFull_P(:, n-1) - delta0;
        bFull_P(:, n) = -delta1 + Phi.' * bFull_P(:, n-1);
        aFull_Q(1, n) = aFull_Q(1, n-1) - delta0 + q0.' * bFull_Q(:, n-1);
        bFull_Q(:, n) = -delta1 + qx.' * bFull_Q(:, n-1);        
    end % if
end % for

% Convert to yields. These are log prices, so we need to invert the yield
% to log price transformation performed earlier.
proportion2percent = 1 / percent2proportion;
A = (-1) * proportion2percent * aFull ./ maturities;
B = (-1) * proportion2percent * bFull ./ maturities;
AP = (-1) * proportion2percent * aFull_P ./ maturities;
BP = (-1) * proportion2percent * bFull_P ./ maturities;
AQ = (-1) * proportion2percent * aFull_Q ./ maturities;
BQ = (-1) * proportion2percent * bFull_Q ./ maturities;

%% Step 6: Prepare the outputs.
rebasedFactors = proportion2percent * ...
    (yieldsMean(1:numFactors) + yieldsStd(1:numFactors) .* selectedFactors);
lambda_t = sqrtSigma \ (lambda0 + lambda1 * rebasedFactors(1:end-1, :).');
lambda_t = lambda_t.';
decomposition.MarketPriceOfRisk = lambda_t;
onesNumTimes = ones( numTimes, 1 );
decomposition.Fitted = onesNumTimes * A + rebasedFactors * B;
decomposition.Expected = onesNumTimes * AP + rebasedFactors * BP;
decomposition.TermPremium = onesNumTimes * (AQ-AP) + rebasedFactors * (BQ-BP);
decomposition.Convexity = onesNumTimes * (A-AQ);

end