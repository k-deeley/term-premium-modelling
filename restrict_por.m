function [parameters,decomposition,returns] = restrict_por(nx,ny,T,nrx,factors,model,parameters,returns)

% Check that the model type is correct
if ~strcmpi(model,'meanzero')
    error('Model must be mean zero model for restrictions on price of risk')
end

% Wrap parameters
thetaHat = [parameters.lambda0(:);parameters.lambda1(:);parameters.beta(:);...
         parameters.hx(:);parameters.sSq;vech(parameters.sigma)];

% Weighting matrix
W = inv(parameters.vTheta);

% Set linear restrictions - only second and third factors priced
H = zeros(nx^2-2*nx,size(thetaHat,1));
H(1:nx,nx+1:2*nx) = eye(nx);
if nx > 3
    H(nx+1:end,4*nx+1:nx^2+nx) = eye(nx*(nx-3));
end

% Get the minimum distance estimator
thetaMD = thetaHat - inv(W)*H'*inv(H*inv(W)*H')*H*thetaHat;

% Re-wrap parameters
parameters.lambda0 = thetaMD(1:nx);
parameters.lambda1 = reshape(thetaMD(nx+1:nx+nx^2),nx,nx);
parameters.beta = reshape(thetaMD(nx+nx^2+1:nx+nx^2+nx*nrx),nx,nrx);
parameters.hx = reshape(thetaMD(nx+nx^2+nx*nrx+1:nx+nx^2+nx*nrx+nx^2),nx,nx);
parameters.sSq = thetaMD(nx+nx^2+nx*nrx+nx^2+1);
sigmaTmp = thetaMD(end-nx*(nx+1)/2+1:end);
parameters.sigma = zeros(nx);
countParams = 1;
for j = 1:nx
    for i = j:nx
        parameters.sigma(i,j) = sigmaTmp(countParams);
        countParams = countParams + 1;
        if j ~= i
            parameters.sigma(j,i) = parameters.sigma(i,j);
        end
    end
end

% Re-do bond pricing, including decompositions
decomposition = bond_pricing(ny,nx,T,factors,parameters.delta0,parameters.delta1,parameters.lambda0,parameters.lambda1,parameters.h0,parameters.hx,parameters.sigma,parameters.sSq);

% Re-do calculation of fitted returns
returns.rxHat = parameters.beta'*(parameters.lambda0*ones(1,T-1)+parameters.lambda1*factors(:,1:T-1))...
    - 0.5.*(parameters.bStar*parameters.sigma(:) + parameters.sSq*ones(nrx,1))*ones(1,T-1)...
    + parameters.beta'*parameters.V;
    
end