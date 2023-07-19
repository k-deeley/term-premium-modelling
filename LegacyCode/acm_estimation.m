function [parameters,decomposition,returns] = acm_estimation(nx,ny,nrx,T,factors,yData,matsAll,rxMats,model,optimiseDeltas)
%************************************************************************************************
% Step 1: estimate the P dynamics using a VAR
%************************************************************************************************
% Estimate equation (1) in ACM
[h0,hx,sigma,V]  = estimate_p_dynamics(nx,T,factors,model);

% Save output
parameters.h0    = h0;
parameters.hx    = hx;
parameters.sigma = sigma;
parameters.V     = V;

%************************************************************************************************
% Step 2: Estimate excess return regressions
%************************************************************************************************
% Convert yields into log prices
pData      = -repmat(matsAll',1,T).*yData/1200;

% Short-term interest rate
stir       = yData(1,:)/1200;

% Compute excess returns - equation (6) in ACM
rx         = pData(1:end-1,2:end) - pData(2:end,1:end-1) - repmat(stir(1,1:end-1),ny-1,1);
rx         = rx(rxMats-1,:);
returns.rx = rx;

% OLS estimator - equation (15) in ACM
Z         = [ones(T-1,1),V',factors(:,1:end-1)']';
allCoeffs = rx/Z;
a         = allCoeffs(:,1);
betaT     = allCoeffs(:,2:nx+1);
c         = allCoeffs(:,nx+2:end);

% Compute the residual variance - text following equation (15) in ACM
e             = rx - allCoeffs*Z;
returns.rxHat = allCoeffs*Z;
sSq           = trace(e*e')/(nrx*T);

% Save output
parameters.sSq  = sSq;
parameters.beta = betaT';

%************************************************************************************************
% Step 3: Prices of risk
%************************************************************************************************
% Construct the B* matrix - equation (13) in ACM
bStar = zeros(nrx,nx^2);
for i = 1:nrx
    bStar(i,:)  = kron(betaT(i,:),betaT(i,:))';
end

% Solve for prices of risk - equations (16) and (17) in ACM
lambda0 = betaT\(a + 0.5.*(bStar*sigma(:)+sSq*ones(nrx,1)));
lambda1 = betaT\c;

% Save output
parameters.lambda0 = lambda0;
parameters.lambda1 = lambda1;
parameters.bStar = bStar;

%************************************************************************************************
% Short rate regression
%************************************************************************************************
% Regress short rates on a constant and the factors.
deltas = [ones(T,1),factors']\stir';
delta0 = deltas(1);
delta1 = deltas(2:end);

% Choose delta0 and delta1 optimally for the whole yield curve.  (NB: This
% extends the ACM method.  They just pick delta0 and delta1 to maximise fit
% to the short rate.  Here we pick them to maximise the fit to all yields.)
if optimiseDeltas == 1
    options = optimset('display','iter');
    deltas = fminsearch(@(deltas)yield_fit(deltas,ny,nx,T,factors,lambda0,lambda1,h0,hx,sigma,sSq,yData),deltas,options);    
    delta0 = deltas(1);
    delta1 = deltas(2:end);
end

% Save output
parameters.delta0 = delta0;
parameters.delta1 = delta1;

%************************************************************************************************
% Bond pricing, including decompositions
%************************************************************************************************
decomposition = bond_pricing(ny,nx,T,factors,delta0,delta1,lambda0,lambda1,h0,hx,sigma,sSq);

end