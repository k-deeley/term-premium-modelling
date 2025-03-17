function decomposition = fitACM( yields, shortTermInterestRate, maturities, nameValArgs )
%FITACM Fit the ACM term premium model for the given yields.
%
% Inputs:
% * yields - a timetable of zero-coupon bond yields.
% * shortTermInterestRate - a timetable containing the short-term interest
%   rate.
% * maturities - a numeric vector of bond maturities, expressed in months.
%   This should be on a month-by-month basis for best results.
% * excessReturnMaturities - a numeric vector of bond maturities, expressed
%   in months. The maturities for which excess returns are computed may be
%   different to the maturities available for the provided yield data,
%   e.g., excess returns may be computed every 6 months instead of monthly.
% * factorMaturities - a numeric vector of bond maturities, expressed
%   in months. The maturities used to compute the PCA factors. This may be
%   exclude some maturities from the yield data, eg. the first six months
%   may be excluded.
% * numBootstrapSamples - an integer, specifying the number of bootstrap
%   samples to use if bootstrapping the VAR model. If this is <= 0, then no
%   bootstrapping is performed.
%
% References:
%
% [1] Pricing the Term Structure with Linear Regressions, Tobias Adrian,
% Richard K. Crump, and Emanuel Moench, Federal Reserve Bank of New York
% Staff Reports, no. 340, August 2008; revised April 2013, JEL
% classification: G10, G12.
% [2] Bond Risk Premia, John H. Cochrane and Monika Piazzesi,
% American Economic Review, Vol. 95, No. 1, March 2005 (pp. 138-160).

arguments ( Input )
    yields(:, :) timetable
    shortTermInterestRate(:, 1) timetable
    maturities(1, :) double {mustBePositive, mustBeInteger}
    nameValArgs.ExcessReturnMaturities(1, :) double {mustBePositive, mustBeInteger} = maturities
    nameValArgs.FactorMaturities(1, :) double {mustBePositive, mustBeInteger} = maturities
    nameValArgs.NumBootstrapSamples(1,:) double {mustBeInteger} = -1
    nameValArgs.BootstrapSampleDistribution(1, 1) string = "saved"
end % arguments ( Input )

excessReturnMaturities = nameValArgs.ExcessReturnMaturities;
factorMaturities = nameValArgs.FactorMaturities;

%% Prepare and validate the input data.
% Divide rt by 100.
[yieldDates, y, rt, maturities] = ...
    prepareInputData( yields, shortTermInterestRate, maturities, excessReturnMaturities, factorMaturities );

% Write down the dimensions of the yield matrix.
% We use the notation from [1].
[T, ~] = size( y );
% We distinguish between N for all maturities available (month-by-month)
% and N for the number of excess return maturities (may be every 6 months)
Nall = numel(maturities); % Keep this as all maturities
Nexr = numel(excessReturnMaturities);

% * T is the number of observations/data points in each yield curve, which
% is equal to the number of rows of the yield matrix.
% * N is the number of maturities, which is equal to the number of columns
% of the yield matrix and the length of the maturities vector.
% * y is the yield matrix of size T-by-N.
% * rt is the vector of short-term interest rates of size T-by-1.
% * yieldDates is the vector of dates corresponding to the yields and
% short-term interest rates, and has size T-by-1.

%% Compute the principal components of the yield curve matrix.
% First, de-mean the yields by subtracting the mean value from each yield
% curve (i.e., each column).
yieldMeans = mean( y );
zeroMeanYields = y - yieldMeans;
% Next, compute the principal components, turning off the default
% auto-centering behavior via the "Centered" option.
[~, factors] = pca( zeroMeanYields(:,factorMaturities), "Centered", false );
% Extract the required number of factors.
% * K is the number of principal components used in the model.
% * X is a T-by-K matrix containing the first K principal components of the
% yield matrix.
K = 4;
X = factors(:, 1:K);
% The first three principal component factors are interpreted as follows.
%
% * Level
% * Slope
% * Curvature

%% Estimate equation (1) from [1] with a vector-autoregressive model.
if nameValArgs.NumBootstrapSamples <= 0
    [mu,Phi,Sigma,V] = PDynamics(X);
else
    [mu,Phi,Sigma,V] = PDynamics(X, "bootstrap", "NumBootstrapSamples", nameValArgs.NumBootstrapSamples, "BootstrapSampleDistribution", nameValArgs.BootstrapSampleDistribution);
end

%% Fit the excess returns using equation (14) from [1].
% Equation (14) from [1] regresses the excess returns (rx) on a constant,
% the lagged pricing factors (X), and the contemporaneous pricing factor
% innovations (V) from the VAR model in the previous step.

% First, we need to convert the observed percentage yields into log prices.
% We use the following from [2]:
%
%   y_t^{(n)} = -1/n * p_t^{(n)}
%
% where:
%   * y_t^{(n)} is the continuously compounded yield (also referred to as
%     the "log yield").
%   * p_t^{(n)} is the log price of an n-year discount bond at time t.
%   * n is the maturity.
% If we are computing the continuously compounded yield, taking the log is
% not necessary.
% See also:
% https://quant.stackexchange.com/questions/28426/how-to-derive-the-relationship-between-log-yield-and-log-price

% In this case we divide by 12 to convert from monthly maturities to annual
% maturities, and by 100 to convert from percentages to proportions.
conversionFactor = 12 * 100;
lnP = (-1) * maturities .* y / conversionFactor;
% * lnP is a T-by-Nall matrix of log prices.

% Next, compute the excess returns using equation (6) of [1].
% Here, the rows represent continuous time (e.g., monthly observations)
% and the columns represent discrete time (bond maturity periods).
% Translating equation (6) of [1] to our matrix lnP:
% * t represents the row index and continuous time.
% * n represents the column index and the discrete bond maturity period.

% Compute excess returns from 12 months onwards (exclude first two entries)
% First row and final column all zeros because we're computing the
% difference
rx = zeros( T, Nall );
rx(2:end, 1:end-1) = lnP(2:end, 1:end-1) - lnP(1:end-1, 2:end) - rt(1:end-1);
% * rx is a T-by-Nall matrix of excess returns.

% Subselect so we only use the excessReturnMaturities specified, i.e. every
% six months
rx = rx(:, excessReturnMaturities-1);

% Implement the linear regression expressed in equation (14) of [1].
% Write down the design matrix Z of size T-by-(2*K+1). We need to pad the
% model residuals V (the contemporaneous pricing factor innovations) and
% the lagged pricing factors (X(1:end-1, :)) with zeros in the first row to
% match the required matrix dimensions.
Z = [ones( T, 1 ), [zeros( 1, K ); V], [zeros(1, K); X(1:end-1, :)]];
allCoeffs = Z \ rx; % (2*K+1)-by-Nexr
a = allCoeffs(1, :).'; % Nexr-by-1
beta = allCoeffs(2:K+1, :).'; % Nexr-by-K
c = allCoeffs(K+2:end, :).'; % Nexr-by-K
% Compute the residuals E (a T-by-N matrix).
E = rx - Z * allCoeffs;
% Estimate sigma^2 (a scalar value).
sigma2 = trace( E * E.' ) / (Nexr * T);
assert (Nexr == width(allCoeffs))

%% Estimate the price of risk parameters via cross-sectional regression.
% Construct the B* matrix (defined below equation (13) in [1]).
Bstar = zeros( Nexr, K^2 );
for n = 1 : Nexr % For each maturity...
    betan = beta(n, :); % 1-by-K
    betanOuterProduct = betan.' * betan; % K-by-K
    % Flatten the outer product into a 1-by-K^2 row vector.
    Bstar(n, :) = reshape( betanOuterProduct, 1, [] );
end % for

% Compute the price of risk parameter lambda0 using equation (16) of [1].
vecSigma = Sigma(:); % K^2-by-1
lambda0 = beta \ (a + 0.5 * (Bstar * vecSigma + sigma2));
% * lambda0 is a K-by-1 risk parameter vector.

% Compute the price of risk parameter lambda1 using equation (17) of [1].
lambda1 = beta \ c;
% * lambda1 is a K-by-K matrix of risk parameters.

%% Estimate the short rate regression parameters delta0 and delta1.
% These parameters are used to start the recursive estimation of the bond
% pricing parameters A_n and B_n in the next step.
deltaDesignMatrix = [ones( T, 1 ), X];
% The design matrix has size T-by-(K+1).
deltas = deltaDesignMatrix \ rt; % A (K+1)-by-1 vector of parameters.
delta0 = deltas(1); % A scalar value (1-by-1).
delta1 = deltas(2:end); % A K-by-1 vector of parameters.

%% Estimate the bond prices via recursive estimation.
% Define the following expressions for convenience in the following
% calculations.
q0 = mu - lambda0; % A K-by-1 vector.
qx = Phi - lambda1; % A K-by-K matrix.

% Preallocate space for the outputs. We evaluate three versions of the "A"
% and "B" matrices according to equations (25) and (26) in [1].
%
% * AP and BP correspond directly to equations (25) and (26).
% * AQ and BQ are the risk-neutral versions of equations (25) and (26), in
% which we impose lambda0 = 0 and lambda1 = 0. Recall that lambda0 and
% lambda1 are the risk parameters, so setting these to zero implies that we
% are risk-neutral. The term premium is calculated as the difference
% between the risk-neutral yield and the model-implied fitted yield.
% * AR and BR are the convexity-neutral versions of equations (25) and
% (26), in the sense that we set the convexity term in equation (25)
% to 0 by imposing Sigma = 0 and sigma2 = 0. The convexity is calculated as
% the difference between the convexity-neutral yield and the model-implied
% fitted yield.
AP = NaN( 1, Nall );
BP = NaN( K, Nall );
AQ = NaN( 1, Nall );
BQ = NaN( K, Nall );
AR = NaN( 1, Nall );
BR = NaN( K, Nall );

for n = 1 : Nall % For each bond maturity period...
    if n == 1
        % Initialize.
        AP(1, n) = -delta0;
        BP(:, n) = -delta1;
        AQ(1, n) = -delta0;
        BQ(:, n) = -delta1;
        AR(1, n) = -delta0;
        BR(:, n) = -delta1;
    else
        % Update recursively.
        % A and B.
        AP(1, n) = AP(1, n-1) + q0.' * BP(:, n-1) + ...
            0.5 * (BP(:, n-1).' * Sigma * BP(:, n-1) + sigma2) - delta0;
        BP(:, n) = qx.' * BP(:, n-1) - delta1;
        % AP and BP (setting lambda0 = 0 and lambda1 = 0). If lambda0 = 0
        % then q0 = mu, and if lambda1 = 0 then qx = Phi.
        AQ(1, n) = AQ(1, n-1) + mu.' * BQ(:, n-1) + ...
            0.5 * (BQ(:, n-1).' * Sigma * BQ(:, n-1) + sigma2) - delta0;
        BQ(:, n) = Phi.' * BQ(:, n-1) - delta1;
        % AQ and BQ (setting Sigma = 0 and sigma2 = 0).
        AR(1, n) = AR(1, n-1) + q0.' * BR(:, n-1) - delta0;
        BR(:, n) = qx.' * BR(:, n-1) - delta1;
    end % if
end % for

%% Convert from log prices to yields.
% This implements the inverse of the yield to price conversion performed
% above.
price2yield = @( lnP ) (-1) * conversionFactor * lnP ./ maturities;
AP = price2yield( AP );
BP = price2yield( BP );
AQ = price2yield( AQ );
BQ = price2yield( BQ );
AR = price2yield( AR );
BR = price2yield( BR );

%% Prepare the outputs.
% Define a helper function to convert from a matrix to a timetable.
ts2tt = @( ts ) array2timetable( ts, "RowTimes", yieldDates, ...
    "VariableNames", string( maturities ) );
% Evaluate the market price of risk (equation (5) of [1]).
sqrtSigma = sqrtm( Sigma );
lambda_t = sqrtSigma \ (lambda0 + lambda1 * X.');
lambda_t = lambda_t.';
decomposition.MarketPriceOfRisk = array2timetable( lambda_t, ...
    "RowTimes", yieldDates, ...
    "VariableNames", "PC" + (1:K) );
% The fitted yields (equation (21) of [1]).
onesNumTimes = ones( T, 1 );
fitted = onesNumTimes * AP + X * BP; % A T-by-N matrix.
decomposition.Fitted = ts2tt( fitted );
% The (risk-neutral) expected values.
riskNeutral = onesNumTimes * AQ + X * BQ; % A T-by-N matrix.
decomposition.RiskNeutralExpected = ts2tt( riskNeutral );
% The term premium: this is calculated as the difference between the
% risk-neutral yield and the model-implied fitted yield.
termPremium = onesNumTimes * (AP - AQ) + X * (BP - BQ);
decomposition.TermPremium = ts2tt( termPremium );
% The convexity: this is calculated as the difference between the
% convexity-neutral yield and the model-implied fitted yield.
convexity = onesNumTimes * (AP - AR) + X * (BP - BR);
decomposition.Convexity = ts2tt( convexity );

end % fitACM

function [yieldDates, yieldMatrix, shortTermInterestRate, maturities] = ...
    prepareInputData( yields, shortTermInterestRate, maturities,  ...
    excessReturnMaturities, factorMaturities)
%PREPAREINPUTDATA Prepare and validate the input data.

% Prepare the dates.
yieldDates = yields.Properties.RowTimes;
mustBeA( yieldDates, "datetime" )
mustBeFinite( yieldDates )

% Prepare the yields.
yieldMatrix = yields{:, :};
mustBeA( yieldMatrix, "double" )
mustBeReal( yieldMatrix )
mustBeFinite( yieldMatrix )

% Validate the short term interest rate dates.
rateDates = shortTermInterestRate.Properties.RowTimes;
mustBeA( rateDates, "datetime" )
mustBeFinite( yieldDates )
assert( isequal( yieldDates, rateDates ), "fitACM:DateMismatch", ...
    "The yield dates must be equal to the dates used for the " + ...
    "short-term interest rate." )

% Prepare the short term interest rates.
shortTermInterestRate = shortTermInterestRate{:, 1};
mustBeA( shortTermInterestRate, "double" )
mustBeReal( shortTermInterestRate )
mustBeFinite( shortTermInterestRate )

% Convert the short term interest rates from percentage yields to
% proportions.
pc2prop = 1/100;
% Convert the short term interest rate from annual to monthly
year2month = 1/12;
shortTermInterestRate = shortTermInterestRate * pc2prop * year2month;

% Validate the maturities.
assert( isequal( length( maturities ), width( yieldMatrix ) ), ...
    "fitACM:MaturityMismatch", ...
    "The number of maturities must match the number of yield curves." )

assert( length(excessReturnMaturities) <=  width( yieldMatrix ), ...
    "fitACM:MaturityMismatch", ...
    "The number of excess return maturities cannot exceed the number of yield curves." )

assert( length(factorMaturities) <=  width( yieldMatrix ), ...
    "fitACM:MaturityMismatch", ...
    "The number of factor maturities cannot exceed the number of yield curves." )

assert( excessReturnMaturities(end) <= maturities(end), ...
    "fitACM:MaturityMismatch", ...
    "The final excess return maturity cannot exceed the final maturity."  )

assert( factorMaturities(end) <= maturities(end), ...
    "fitACM:MaturityMismatch", ...
    "The final factor maturity cannot exceed the final maturity."  )

end % prepareInputData