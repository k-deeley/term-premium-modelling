function vTheta = compute_covariance_matrix(parameters,factors,nx,nrx,T)

%************************************************************************************************
% Extract parameters
%************************************************************************************************
lambda = [parameters.lambda0,parameters.lambda1];
beta = parameters.beta;
sSq = parameters.sSq;
sigma = parameters.sigma;
bStar = parameters.bStar;

%************************************************************************************************
% Define useful matrices
%************************************************************************************************
% (nx+1) x 1 vector with first element one
d1 = zeros(nx+1,1);
d1(1) = 1;

% aBeta matrix defined in Appendix A.1
aBeta = zeros(nrx*nx,nrx);
for i = 1:nrx
    aBeta((i-1)*nx+1:i*nx,i) = beta(:,i);
end

% The duplication matrix such that vec(A) = G*vech(A) and its Moore-Penrose
% pseudo-inverse
G = duplication_matrix(nx);
GPlus = pinv(G);

% Estimate the limiting distributions of the squared factors
% TO DO: CHECK THIS IS THE CORRECT WAY TO ESTIMATE THESE TERMS.
yxx = factors*factors'/T;
yzz = [ones(1,T);factors]*[ones(1,T);factors]'/T;

%************************************************************************************************
% V_lambda
%************************************************************************************************
% V(lambda,tau,1) - equation (102)
vLambdaTau1 = kron(inv(yzz),eye(nx))*kron([0,zeros(1,nx);zeros(nx,1),yxx],sigma)*kron(inv(yzz),eye(nx))';

% V(lambda,tau,2) - equation (104)
vLambdaTau2 = sSq.*kron(inv(yzz),inv(beta*beta'));

% V(lambda,tau,3) - equation (105)
vLambdaTau3 = sSq.*kron(lambda'*inv(sigma)*lambda,inv(beta*beta'));

% V(lambda,tau,4) - equation (106)
vLambdaTau4 = sSq.*kron(d1*d1',inv(beta*beta')*beta*aBeta'*kron(eye(nrx),sigma)*aBeta*beta'*inv(beta*beta'));

% V(lambda,tau,5) - equation (108)
vLambdaTau5 = 0.25.*kron(d1*d1',inv(beta*beta')*beta*bStar*(eye(nx^2)+vecperm(nx,nx))*kron(sigma,sigma)*bStar'*beta'*inv(beta*beta'));

% V(lambda,tau,5) - equation (109)
vLambdaTau6 = 0.5.*sSq^2.*kron(d1*d1',inv(beta*beta')*beta*ones(nrx,1)*ones(1,nrx)*beta'*inv(beta*beta'));

% Sum of the terms
vLambdaTau = vLambdaTau1 + vLambdaTau2 + vLambdaTau3 + vLambdaTau4 + vLambdaTau5 + vLambdaTau6;

% C(lambda,tau) - equation (101)
cLambdaTau = -kron(lambda',inv(beta*beta')*beta)*vecperm(nx,nrx)*(sSq.*kron(eye(nrx),sigma))*(kron(d1,inv(beta*beta')*beta*aBeta'*kron(eye(nrx),sigma)))';

% V_lambda
vLambda = vLambdaTau + cLambdaTau + cLambdaTau';

%************************************************************************************************
% Remaining covariances
%************************************************************************************************
% V(beta) - equation (98) and text on p.46
vBeta = sSq.*kron(eye(nrx),inv(sigma));

% C(lambda,beta) - equation (115)
cLambdaBeta = -sSq.*vecperm(nx+1,nx)*kron(inv(beta*beta')*beta,lambda'*inv(sigma)) + sSq.*kron(d1,inv(beta*beta')*beta*aBeta');

% C(lambda,phi) - equation (117)
tmpC = kron(yzz,sigma);
cLambdaPhi = kron(inv(yzz),eye(nx))*tmpC(:,end-nx^2+1:end)*kron(inv(yxx),eye(nx));

% C(lambda,sSq) - equation (118)
cLambdaSSq = sSq^2.*kron(d1,inv(beta*beta')*beta)*ones(nrx,1);

% C(lambda,sigma) - equation (121)
% cLambdaSigma = 0.5.*kron(d1,inv(beta*beta')*beta*bStar)*G*GPlus*(eye(nx^2)+vecperm(nx,nx))*kron(sigma,sigma)*GPlus';
cLambdaSigma = kron(d1,inv(beta*beta')*beta*bStar)*G*GPlus*kron(sigma,sigma)*GPlus';

% V(phi) - text following equation (121)
vPhi = kron(inv(yxx),sigma);

% V(sSq) - text following equation (121)
vSSq = 2.*sSq^2;

% V(sigma) - text following equation (121)
vSigma = 2.*GPlus*kron(sigma,sigma)*GPlus';

%************************************************************************************************
% Combine into covariance matrix for theta
%************************************************************************************************
% Initialise memory
vTheta = zeros(nx*(nx+1)+nrx*nx+nx^2+1+nx*(nx+1)/2);

% Allocate the variances
vTheta(1:nx*(nx+1),1:nx*(nx+1)) = vLambda;
vTheta(nx*(nx+1)+1:nx*(nx+1)+nx*nrx,nx*(nx+1)+1:nx*(nx+1)+nx*nrx) = vBeta;
vTheta(nx*(nx+1)+nx*nrx+1:nx*(nx+1)+nx*nrx+nx^2,nx*(nx+1)+nx*nrx+1:nx*(nx+1)+nx*nrx+nx^2) = vPhi;
vTheta(nx*(nx+1)+nx*nrx+nx^2+1,nx*(nx+1)+nx*nrx+nx^2+1) = vSSq;
vTheta(nx*(nx+1)+nx*nrx+nx^2+2:end,nx*(nx+1)+nx*nrx+nx^2+2:end) = vSigma;

% Allocate covariances
vTheta(nx*(nx+1)+1:nx*(nx+1)+nx*nrx,1:nx*(nx+1)) = cLambdaBeta';
vTheta(1:nx*(nx+1),nx*(nx+1)+1:nx*(nx+1)+nx*nrx) = cLambdaBeta;
vTheta(nx*(nx+1)+nx*nrx+1:nx*(nx+1)+nx*nrx+nx^2,1:nx*(nx+1)) = cLambdaPhi';
vTheta(1:nx*(nx+1),nx*(nx+1)+nx*nrx+1:nx*(nx+1)+nx*nrx+nx^2) = cLambdaPhi;
vTheta(nx*(nx+1)+nx*nrx+nx^2+1,1:nx*(nx+1)) = cLambdaSSq';
vTheta(1:nx*(nx+1),nx*(nx+1)+nx*nrx+nx^2+1) = cLambdaSSq;
vTheta(nx*(nx+1)+nx*nrx+nx^2+2:end,1:nx*(nx+1)) = cLambdaSigma';
vTheta(1:nx*(nx+1),nx*(nx+1)+nx*nrx+nx^2+2:end) = cLambdaSigma;

% TO DO : CHECK OTHER ELEMENTS SHOULD BE ZERO

end