function [decomposition, lambda0, lambda1, delta0, delta1] = ...
    estimateACM( yields, shortTermInterestRate, maturities, numFactors )

arguments
    % Matrix of yields (%).
    yields(:, :) double
    % Vector of short term interest rates (%).
    shortTermInterestRate(:, 1) double
    % Vector of maturities (months).
    maturities(1, :) double
    % Number of factors (principal components) to use in the model.
    numFactors(1, 1) double = 3
end % arguments

%% Prepare the data.

% Convert the short-term interest rates from percentages to proportions.
percent2proportion = 1 / 100;
shortTermInterestRate = shortTermInterestRate * percent2proportion;

% Compute the principal components of the yields.
[~, scores] = pca( yields );

% Extract the first "numComponents" principal components.
selectedScores = scores(:, 1:numFactors);

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
[VARModel, ~, ~, modelResiduals] = estimate( VARModel, selectedScores );
Phi = VARModel.AR{1};
mu = VARModel.Constant;

% Compute the covariance matrix of the residuals.
Sigma = VARModel.Covariance;

%% Step 2: Estimate excess returns via linear regression.
% Convert yields into log prices. In this case we divide by 12 to go from
% the monthly maturities to annual maturities, and by 100 to obtain a 
% proportion.
monthly2annual = 12;
proportion2percent = 100;
fudge = monthly2annual * proportion2percent;
prices = (-1) * maturities .* yields / fudge;

% Compute excess returns (equation (6) in ACM). If we are getting the
% continously compounded yield, the log is not necessary.
% https://quant.stackexchange.com/questions/28426/how-to-derive-the-relationship-between-log-yield-and-log-price
excessReturns = prices(2:end, 1:end-1) - prices(1:end-1, 2:end) - ...
    shortTermInterestRate(1:end-1);
numExcessReturnSeries = width( excessReturns );

% OLS estimator (equation (15) in ACM).
Z = [ones( numTimes-1, 1 ), modelResiduals, selectedScores(1:end-1, :)];
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
% factors.
% shortTermInterestRate = delta0 + selectedScores * delta1
deltas = [ones( numTimes, 1), selectedScores] \ shortTermInterestRate;
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

% Convert to yields. These are log pricesis is the log price, to convert to yield we need to
% divide by -n, in years which is -1/((1:120)/12)*100
A = (-1) * fudge * aFull ./ maturities;
B = (-1) * fudge * bFull ./ maturities;

AP = (-1) * fudge * aFull_P ./ maturities;
BP = (-1) * fudge * bFull_P ./ maturities;

AQ = (-1) * fudge * aFull_Q ./ maturities;
BQ = (-1) * fudge * bFull_Q ./ maturities;

%% Step 6: Prepare the outputs.
lambda_t = sqrtSigma \ (lambda0 + lambda1 * selectedScores(1:end-1, :).');
lambda_t = lambda_t.';
decomposition.lambda_t = lambda_t;
decomposition.yHat = (ones(numTimes,1) * A + selectedScores * B);
decomposition.expected = (ones(numTimes,1) * AP + selectedScores*BP);
decomposition.riskPremium = (ones(numTimes,1) * (AQ-AP) + selectedScores*(BQ-BP));
decomposition.convexity   = (ones(numTimes,1) * (A-AQ));

end