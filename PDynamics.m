function [mu,Phi,Sigma,V] = PDynamics(X, method, nameValArgs)
%PDynamics: Compute VAR step on factors, optionally using bootstrap
%
% Inputs:
% * X - a T-by-K matrix containing the first K principal components of the
% yield matrix. 
% * method - the method for computing the p dynamics. This should be either
% "default" or "bootstrap". If unspecified, "default" will be used.
% * numBootstrapSamples - a positive integer, specifying the number of bootstrap
% samples to use if bootstrapping the VAR model. Defaults to 100.
% * bootstrapSampleDistribution - distribution to use for sampling residuals; 
% This can be either "saved", to resample with replacement from the given residuals 
% or "fitted" to sample residuals from distributions fitted to the residuals 
% of each factor.
%

arguments
    X (:, :) double
    method (1,1) string {mustBeMember(method,{'default','bootstrap'})} = "default"
    nameValArgs.NumBootstrapSamples(1,:) double {mustBeInteger, mustBePositive} = 1000;
    nameValArgs.BootstrapSampleDistribution (1,1) string {mustBeMember(nameValArgs.BootstrapSampleDistribution,{'saved','fitted'})} = "saved";
end % arguments

% Compute the number of factors
K = width(X);
% Create and estimate a VAR(1) model for the factors.
numLags = 1;
VARModel = varm( K, numLags );

% Calibrate mu = 0 (the constant terms in the VAR model).
mu = zeros( K, 1 );
VARModel.Constant = mu;
% Estimate the model and return the innovations matrix V.
[VARModel, ~, ~, V] = estimate( VARModel, X );
% mu' (1xK) + X (T x K) * Phi' (KxK) + V (T x K)
% * V is a (T-1)-by-K matrix of innovations.
% * Phi is the K-by-K matrix of VAR model coefficients.

% Extract the autoregressive coefficients Phi.
Phi = VARModel.AR{1};
% Extract the covariance matrix Sigma.
% * Sigma is the K-by-K variance-covariance matrix of the model parameters.
Sigma = VARModel.Covariance;

% No bootstrapping; we are done
if method == "default"
    return;
end

bootstrapDist = nameValArgs.BootstrapSampleDistribution;
numBootstrapSamples = nameValArgs.NumBootstrapSamples;

if bootstrapDist == "saved"
    residuals = V;
else
    numResiduals = height(V);
    numSeries = width(V);
    residualFit = cell(1, numSeries);
    residuals = NaN(numResiduals, numSeries);
	% Fit a distribution to each factor's residuals
    for k = 1:numSeries
        residualFit{k} = fitdist(V(:, k), "kernel");
		% Sample required number of residuals
        residuals(:, k) = random(residualFit{k}, [numResiduals, 1]);
    end
end


% phiBootstrap will contain the estimated Phi for each bootstrap
% sample
phiBootstrap = nan(numBootstrapSamples, width(Phi)*height(Phi));
% Total number of time steps in X
T = height(X);

for s = 1:numBootstrapSamples
    % We will need to simulate for T timesteps, so choose T random
    % indices from X.
    % To be consistent with V definition, we have one
    % fewer V than T.
    indexMC = randi(T-1, T, 1);
    % Initialize bootstrap sample at time t=1
    XMC(1,:) = X(indexMC(1), :);
    for t = 2:T
        % Calculate next time step as Phi * X_t + randomly sampled
        % residual
        XMC(t,:) = (XMC(t-1,:)) * Phi' + residuals(indexMC(t), :) + mu';
    end
    
    % Estimate new Phi using bootstrap sample
    VARModelMC = varm( K, numLags );
    VARModelMC.Constant = zeros( K, 1 );
    VARModelMC = estimate( VARModelMC, XMC );
    PhiMC = VARModelMC.AR{1};
    phiBootstrap(s, :) = PhiMC(:);
end

% Compute the average Phi over all bootstrap samples
PhiBar = mean(phiBootstrap); 
PhiBar = reshape(PhiBar, size(Phi));
bias = PhiBar - Phi;

% Apply the bias correction
PhiBC = Phi - bias;

% Enforce stationarity by shrinking towards the OLS estimate
scaling = 1;
while any(abs(eig(PhiBC)) >= 1)
    % Re-scale
    scaling = scaling*0.999;
    PhiBC = Phi - scaling*bias;
end

% Compute V with the new, bias-corrected Phi
% We essentially want to compute LHS - Phi X_t to get V
% i.e. X_t+1 - Phi * X_t 
Xtp1 = X(2:T,:);
Xt = X(1:T-1,:);
VMC = Xtp1 - Xt * PhiMC';
% Defaults to normalizing by num observations -1
% This is closer to the Sigma from varm than normalizing by num
% observations
SigmaMC = cov(VMC);
% Overwrite empirical values to return
V = VMC;
Sigma = SigmaMC;
