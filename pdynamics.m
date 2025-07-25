function [mu, Phi, Sigma, V] = pdynamics( X, method, namedArgs )
%PDYNAMICS Estimate pricing factor dynamics for the ACM term premium model.
%
%  [mu, Phi, Sigma, V] = pdynamics( X, method ) estimates the pricing
%  factor dynamics for the given matrix X using the method "deterministic"
%  or "bootstrap".
%
%  Input Arguments
%
% * X is a T-by-K matrix containing the first K principal components of the
%   yield matrix, where T is the number of observations in time.
%
% * method represents the method for computing the p-dynamics and is either
%   "deterministic" (in which case no bootstrapping is used to estimate the
%   model parameters) or "bootstrap". The latter method uses bootstrap
%   sampling from the residuals of the fitted VAR model to provide updated
%   estimates of the model parameters. If unspecified, "deterministic" will
%   be used.
%
% Name-Value Arguments
%
% * NumBootstrapSamples: a positive integer specifying the number of
%   bootstrap samples to use when method is "bootstrap". The default value
%   is 1000.
%
% * BootstrapSampleDistribution: the probability distribution to use when
%   sampling the initial model residuals. This can be either be
%   "empirical", which samples with replacement from the inferred residual
%   series, or "fitted", which samples the residuals from kernel-smoothed
%   distributions fitted to the residuals of each factor. The default value
%   is "empirical".
%
% Output Arguments
%
% * mu: real, finite, K-by-1 column vector of VAR model constants.
%
% * Phi: real, finite, K-by-K matrix of VAR model autoregressive
%   coefficients.
%
% * Sigma: real, finite, positive definite K-by-K covariance matrix of VAR
%   model coefficients.
%
% * V: real, finite T-by-K matrix of VAR model residuals (innovations).
%
% See also fitACM

arguments ( Input )
    X(:, :) double {mustBeReal, mustBeFinite}
    method(1, 1) string {mustBeMember( method, ...
        ["deterministic", "bootstrap"] )} = "deterministic"
    namedArgs.NumBootstrapSamples(1, 1) double ...
        {mustBeInteger, mustBePositive} = 1000
    namedArgs.BootstrapSampleDistribution(1, 1) string ...
        {mustBeMember( namedArgs.BootstrapSampleDistribution, ...
        ["empirical", "fitted"] )} = "empirical"
end % arguments ( Input )

arguments ( Output )
    mu(:, 1) double {mustBeReal, mustBeFinite}
    Phi(:, :) double {mustBeReal, mustBeFinite, mustBeSquare}
    Sigma(:, :) double {mustBeReal, mustBeFinite, mustBeSquare}
    V(:, :) double {mustBeReal, mustBeFinite}
end % arguments ( Output )

% Number of time steps and factors (principal components).
[T, K] = size( X, 1:2 );

% Create and estimate a VAR(1) model for the factors.
numLags = 1;
VARModel = varm( K, numLags );

% Calibrate mu = 0 (the constant terms in the VAR model).
mu = zeros( K, 1 );
VARModel.Constant = mu;

% Estimate the model and return the residual (innovation) matrix V, of size
% T-by-K.
[VARModel, ~, ~, V] = estimate( VARModel, X );

% Extract the K-by-K matrix of autoregressive coefficients Phi.
Phi = VARModel.AR{1};

% Extract the K-by-K covariance matrix of model parameters Sigma.
Sigma = VARModel.Covariance;

% If there's no bootstrapping required, we're done.
if method == "default"
    return
end % if

% Extract the named arguments.
bootstrapDist = namedArgs.BootstrapSampleDistribution;
numBootstrapSamples = namedArgs.NumBootstrapSamples;

if bootstrapDist == "empirical"
    residuals = V;
else

    % Preallocate.
    [numResiduals, numSeries] = size( V, 1:2 );
    residualFit = cell( 1, numSeries );
    residuals = NaN( numResiduals, numSeries );

    % Fit a kernel-smoothed distribution to each factor's residuals. Sample
    % the required number of residuals.
    for seriesIdx = 1 : numSeries
        residualFit{seriesIdx} = fitdist( V(:, seriesIdx), "kernel" );
        residuals(:, seriesIdx) = random( residualFit{seriesIdx}, ...
            [numResiduals, 1] );
    end % for

end % if

% PhiBootstrap will contain the estimated Phi matrix for each bootstrap
% sample.
PhiBootstrap = NaN( numBootstrapSamples, numel( Phi ) );

for s = 1 : numBootstrapSamples

    % We will need to simulate for T timesteps, so choose T random
    % indices from X. To be consistent with the definition of V, we have
    % one less observation in V than in X.
    sampleIdx = randi( T-1, T, 1 );

    % Initialize bootstrap sample at time t = 1.
    simX(1, :) = X(sampleIdx(1), :);

    for t = 2 : T
        % Calculate next time step as Phi * X_t + the randomly sampled
        % residual.
        simX(t, :) = mu' + (simX(t-1, :)) * Phi' + ...
            residuals(sampleIdx(t), :);
    end % for

    % Estimate the new Phi matrix using the bootstrap sample.
    simVARModel = varm( K, numLags );
    simVARModel.Constant = zeros( K, 1 );
    simVARModel = estimate( simVARModel, simX );
    simPhi = simVARModel.AR{1};
    PhiBootstrap(s, :) = simPhi(:);

end % for

% Compute the mean value of Phi over all bootstrap samples.
PhiBar = mean( PhiBootstrap );
PhiBar = reshape( PhiBar, size( Phi ) );
bias = PhiBar - Phi;

% Apply the bias correction.
PhiBC = Phi - bias;

% Enforce stationarity by shrinking towards the OLS estimate if necessary.
scaling = 1;
while any( abs( eig( PhiBC ) ) >= 1 )
    % Shrink
    scaling = 0.999 * scaling;
    PhiBC = Phi - scaling * bias;
end % while

% Compute V with the new, bias-corrected Phi.
% We essentially want to compute LHS - Phi X_t to get V,
% i.e. X_t+1 - Phi * X_t.
Xtp1 = X(2:T, :);
Xt = X(1:T-1, :);
simV = Xtp1 - Xt * simPhi';

% Estimate the new covariance matrix.
simSigma = cov( simV );

% Replace the empirical values.
V = simV;
Sigma = simSigma;

end % pdynamics