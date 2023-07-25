function [decomposition, lambda0, lambda1, delta0, delta1, tsquared] = estimateACM(yields, stir, maturities, numFactors, nvp )

arguments %(Input)
    yields     %(:,:) timetable            % Timetable with yields in %
    stir       %(:,1) timetable            % Time table with the short term interest rate in %
    maturities (1,:) double               % Maturities (in months)
    numFactors (1,1) double = 3           % Number of factors
    nvp.pcMats (1,:) double = maturities  % Maturities to compute the principal components (in months)
    nvp.rxMats (1,:) double = maturities  % Maturities to compute the excessReturns (in months)
    nvp.model  (1,1) string = "meanzero"  % Type of var model
end

% arguments (Output)
%     decomposition (1,1) struct
%     lambda0 (:,1) double
%     lambda1 (:,:) double
%     delta0 (1,1) double
%     delta1 (:,1) double
%     tsquared double
% end

pcMats = nvp.pcMats;
rxMats = nvp.rxMats;

yieldMatrix = yields{:,:};           % Extract all the yields
shortTermInterestRate = stir{:,:}/100;

[T, ny] = size(yieldMatrix);  

[~, scores, ~, tsquared] = pca( yieldMatrix(:,pcMats) );
%[~, scores, ~, tsquared] = pca( yieldMatrix(:, :) );

% Extract the first "numComponents" principal components.
selectedScores = scores(:, 1:numFactors);

%************************************************************************************************
% Step 1: estimate the P dynamics using a VAR
%************************************************************************************************

%% Create and estimate a VAR model for the yield principal components. (Eq 1)
numLags = 1;
VARModel = varm( numFactors, numLags );
VARModel.SeriesNames = "PCA_Factor_" + (1:numFactors);
VARModel.Description = "Vector autoregressive zero-mean " + ...
    "model of order 1 for the yield principal components";

% Set the model constants to be zero.
switch nvp.model
    case "meanzero"
        mu = 0;
        VARModel.Constant = mu * ones( numFactors, 1 );
    case "benchmark"

end

% Estimate the model.
[VARModel, ~, ~, modelResiduals] = ...
    estimate( VARModel, selectedScores );
Phi = VARModel.AR{1};
mu = VARModel.Constant;

% Compute the covariance matrix of the residuals.
covModelResiduals = VARModel.Covariance;

%************************************************************************************************
% Step 2: Estimate excess return regressions
%************************************************************************************************
% Convert yields into log prices. In this case we divide by 12 to go from
% the monthly maturities to the yearly ones, and by 100 to remove the %. Is
% different than the one below.
prices      = -maturities.*yieldMatrix/1200;

% Compute excess returns - equation (6) in ACM. If we are getting the
% continously compounded yield, the log is not necessary.
% https://quant.stackexchange.com/questions/28426/how-to-derive-the-relationship-between-log-yield-and-log-price
excessReturns = prices(2:end, 1:end-1) - prices(1:end-1, 2:end) - shortTermInterestRate(1:end-1);
excessReturns = excessReturns(:,rxMats-1);
%excessReturns = excessReturns(:, maturities(1:end-1));
numExcessReturnSeries = width( excessReturns );

% OLS estimator - equation (15) in ACM
Z         = [ones(T-1,1), modelResiduals, selectedScores(1:end-1,:)];
allCoeffs = Z\excessReturns;
aCoeffs   = allCoeffs(1, :)'; % NMats x 1
bCoeffs   = allCoeffs(2:numFactors+1, :)'; % NMats x NFactors
cCoeffs   = allCoeffs(numFactors+2:end, :)'; % NMats x NFactors

returnPricingErrors = excessReturns - Z*allCoeffs;
normTrace           = trace(returnPricingErrors*returnPricingErrors')/numel( Z );

%************************************************************************************************
% Step 3: Prices of risk
%************************************************************************************************
% Construct the B* matrix - equation (13) in ACM
Bstar = zeros( numExcessReturnSeries, numFactors^2 );
for k = 1 : numExcessReturnSeries
    Bstar(k, :) = kron( bCoeffs(k, :), bCoeffs(k, :) ).';
end % for

%% Compute the market prices of risk (the parameters lambda0 and lambda1).
lambda0 = bCoeffs \ (aCoeffs + 0.5 * (Bstar * covModelResiduals(:) + ...
    normTrace * ones( numExcessReturnSeries, 1 )));
lambda1 = bCoeffs \ cCoeffs;

%************************************************************************************************
% Short rate regression
%************************************************************************************************
% Regress short rates on a constant and the factors.
deltas = [ones(T,1),selectedScores]\shortTermInterestRate;
delta0 = deltas(1);
delta1 = deltas(2:end);

%************************************************************************************************
% Bond pricing, including decompositions
%************************************************************************************************
% Construct risk-neutral pricing coefficients - equations (25) and (26) in
% ACM
q0 = mu - lambda0;
qx = Phi - lambda1;

% Recursive bond prices - equations (25) and (26) in ACM
aFull = nan(1, ny);
bFull = nan(numFactors, ny);

% Bond prices under P with sigma = 0
aFull_P = nan(1,ny);
bFull_P = nan(numFactors,ny);

% Bond prices under Q with sigma = 0
aFull_Q = nan(1,ny);
bFull_Q = nan(numFactors,ny);

for n = 1:ny
    if n == 1
        aFull(1,n) = -delta0;
        bFull(:,n) = -delta1;

        aFull_P(1,n) = -delta0;
        bFull_P(:,n) = -delta1;

        aFull_Q(1,n) = -delta0;
        bFull_Q(:,n) = -delta1;
    else
        aFull(1,n) = aFull(1,n-1) + dot(bFull(:,n-1),q0) + 0.5*(bFull(:,n-1)'*covModelResiduals*bFull(:,n-1) + normTrace) - delta0;
        bFull(:,n) = -delta1 + qx'*bFull(:,n-1);

        aFull_P(1,n) = aFull_P(1,n-1) - delta0;
        bFull_P(:,n) = -delta1 + Phi'*bFull_P(:,n-1);

        aFull_Q(1,n) = aFull_Q(1,n-1) - delta0 + q0'*bFull_Q(:,n-1);
        bFull_Q(:,n) = -delta1 + qx'*bFull_Q(:,n-1);
    end
end

% Convert to yields. This is the log price, to convert to yield we need to
% divide by -n, in years which is -1/((1:120)/12)*100
%TODO: Divide each column by (1:ny)
A = -1200*aFull./(1:ny);
B = -1200*bFull./repmat((1:ny),numFactors,1);

AP = -1200*aFull_P./(1:ny);
BP = -1200*bFull_P./repmat((1:ny),numFactors,1);

AQ = -1200*aFull_Q./(1:ny);
BQ = -1200*bFull_Q./repmat((1:ny),numFactors,1);

% Fitted yields
decomposition.yHat = (ones(T,1) * A + selectedScores * B);
% yHat = ones(T,1)*A' + selectedScores*B';
decomposition.expected = (ones(T,1) * AP + selectedScores*BP);
decomposition.riskPremium = (ones(T,1) * (AQ-AP) + selectedScores*(BQ-BP));
decomposition.convexity   = (ones(T,1) * (A-AQ));

end