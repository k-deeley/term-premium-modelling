%% ACM Model Estimation.

%% Import the zero coupon yields.
filename = "ZeroCouponYieldsMonthly.xlsx";
maturities = readmatrix( filename, "Range", "B1:K1" );
dates = datetime( readmatrix( filename, "Range", "A2:A330" ), ...
    "ConvertFrom", "excel" );
yields = readmatrix( filename, "Range", "B2:K330" );
varNames = "Maturity_" + maturities + "_years";
T = array2timetable( yields, "RowTimes", dates, ...
    "VariableNames", varNames );

%% Remove rows with missing or negative data.
missingIdx = ismissing( T );
badRows = any( missingIdx, 2 ) | any( T.Variables < 0, 2 );
T(badRows, :) = [];
numObservations = height( T );

%% For each maturity, plot the corresponding yields over time.
figure
ax = axes;
plot( ax, T.Properties.RowTimes, T.Variables )
xlabel( ax, "Date" )
ylabel( ax, "Yield (%)" )
title( ax, "Zero-Coupon UK Yields" )
grid( ax, "on" )
leg = legend( ax, string( maturities ), "NumColumns", 2 );
leg.Title.String = "Maturity (years)";
numMaturities = length( maturities );
ax.ColorOrder = jet( numMaturities );

%% Compute the principal components of the yield matrix.
yieldMatrix = T.Variables;
[coeffs, scores, eigVals] = pca( yieldMatrix );

% Extract the first 4 principal components.
numComponents = 4;
selectedScores = scores(:, 1:numComponents);

% The principal component factors are interpreted as follows.
%
% First factor:
%
% * The expected returns of all maturities move together, so a single
% variable x_t (first factor) cab describe the time-variation of the
% expected returns.
%
% The remaining factors are, respectively:
%
% * Level
% * Slope
% * Curvature

%% Create and estimate a VAR model for the yield principal components.
numLags = 1;
numSeries = numComponents;
VARModel = varm( numSeries, numLags );
VARModel.SeriesNames = "PCA_Factor_" + (1:numSeries);
VARModel.Description = "Vector autoregressive zero-mean " + ...
    "model of order 1 for the yield principal components";

% Set the model constants to be zero.
mu = 0;
VARModel.Constant = mu * ones( numSeries, 1 );

% Only estimate the diagonal terms, and element (1, 4).
VARModel.AR{1} = diag( NaN( 1, 4 ) );
VARModel.AR{1}(1, 4) = NaN;

% Estimate the model.
[VARModel, estimatedStandardErrors, ~, modelResiduals] = ...
    estimate( VARModel, selectedScores );
Phi = VARModel.AR{1};

% Compute the covariance matrix of the residuals.
covModelResiduals = cov( modelResiduals );

%% Estimate the excess returns.

% Convert yields into prices.
percent2proportion = 0.01;
monthly2annual = 1/12;
conversionFactor = percent2proportion * monthly2annual;
prices = maturities .* yieldMatrix * conversionFactor;

% Short-term interest rate. This acts as the instantaneous risk-free rate.
shortTermInterestRate = yieldMatrix(:, 1) * conversionFactor;

% Compute the excess returns. Here, the rows represent continuous time
% (monthly data) and the columns represent discrete time (bond maturity
% periods).
%
% t - row index (continuous, monthly)
% n - column index (discrete, annual)
excessReturns = log( prices(2:end, 1:end-1) ) - ...
    log( prices(1:end-1, 2:end) ) - shortTermInterestRate(1:end-1);
numExcessReturns = height( excessReturns );

%% Regress the log excess returns on the factors.
numExcessReturnSeries = width( excessReturns );
aCoeffs = zeros( numExcessReturnSeries, 1 );
bCoeffs = zeros( numExcessReturnSeries, numSeries );
cCoeffs = zeros( numExcessReturnSeries, numSeries );
returnPricingErrors = zeros( numExcessReturns, numExcessReturnSeries );

for k = 1 : numExcessReturnSeries
    designMatrix = [ones( numExcessReturns, 1), modelResiduals, ...
        selectedScores(2:end, :)];
    currentExcessReturn = excessReturns(:, k);
    currentModel = fitlm( designMatrix, currentExcessReturn, ...
        "Intercept", false );
    currentModelCoeffs = currentModel.Coefficients.Estimate;
    aCoeffs(k) = currentModelCoeffs(1); % 1
    bCoeffs(k, :) = currentModelCoeffs(2:(2+numSeries-1)); % 4
    cCoeffs(k, :) = currentModelCoeffs((2+numSeries):end); % 4
    returnPricingErrors(:, k) = currentModel.Residuals.Raw;
end % for

% Compute the normalized trace of the covariance matrix of the residuals 
% from the linear regression models.
normTrace = trace( returnPricingErrors * returnPricingErrors.' ) / ...
    numel( designMatrix );

%% Compute Bstar from equation (13).
Bstar = zeros( numExcessReturnSeries, numSeries^2 );
for k = 1 : numExcessReturnSeries
    Bstar(k, :) = kron( bCoeffs(k, :), bCoeffs(k, :) ).';
end % for

%% Compute the market prices of risk (the parameters lambda0 and lambda1).
% aCoeffs - 9-by-1
% bCoeffs - 9-by-4
% cCoeffs - 9-by-4
% Bstar - 9-by-16
% covModelResiduals - 4-by-4
% lambda0 - 4-by-1
% lambda1 - 4-by-4
lambda0 = bCoeffs \ (aCoeffs + 0.5 * (Bstar * covModelResiduals(:) + ...
    normTrace * ones( numExcessReturnSeries, 1 )));
lambda1 = bCoeffs \ cCoeffs;

%% Model the short rates as a linear function of the factors.
shortRateModel = fitlm( selectedScores, shortTermInterestRate );
delta0 = shortRateModel.Coefficients.Estimate(1);
delta1 = shortRateModel.Coefficients.Estimate(2:end);

%% Price the bonds recursively.
q0 = mu - lambda0;
qx = Phi - lambda1;
A = NaN( numExcessReturnSeries, 1 );
B = NaN( numExcessReturnSeries, numComponents );

for n = 1 : numExcessReturnSeries
    if n == 1
        A(n, 1) = -delta0;
        B(n, :) = -delta1.';
    else
        A(n, 1) = A(n-1, 1) + B(n-1, :) * q0 + ...
            0.5 * (B(n-1, :) * covModelResiduals * B(n-1, :)' + ...
            normTrace) - delta0;
        B(n, :) = B(n-1, :) * qx - delta1.';
    end % if
end % for