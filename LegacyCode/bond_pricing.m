function decomposition = bond_pricing(ny,nx,T,factors,delta0,delta1,lambda0,lambda1,h0,hx,sigma,sSq)

% Construct risk-neutral pricing coefficients - equations (25) and (26) in
% ACM
q0 = h0 - lambda0;
qx = hx - lambda1;

% Recursive bond prices - equations (25) and (26) in ACM
aFull = nan(ny,1);
bFull = nan(ny,nx);
for n = 1:ny
    if n == 1
        aFull(n,1) = -delta0;
        bFull(n,:) = -delta1';
    else
        aFull(n,1) = aFull(n-1,1) - delta0 + bFull(n-1,:)*q0 + 0.5*(bFull(n-1,:)*sigma*bFull(n-1,:)' + sSq);
        bFull(n,:) = -delta1' + bFull(n-1,:)*qx;
    end
end

% Convert to yields
A = -1200*aFull./(1:ny)';
B = -1200*bFull./repmat((1:ny)',1,nx);

% Fitted yields
decomposition.yHat = A*ones(1,T) + B*factors;

% Bond prices under P with sigma = 0
aFull = nan(ny,1);
bFull = nan(ny,nx);
for n = 1:ny
    if n == 1
        aFull(n,1) = -delta0;
        bFull(n,:) = -delta1';
    else
        aFull(n,1) = aFull(n-1,1) - delta0 + bFull(n-1,:)*h0;
        bFull(n,:) = -delta1' + bFull(n-1,:)*hx;
    end
end

% Convert to yields
AP = -1200*aFull./(1:ny)';
BP = -1200*bFull./repmat((1:ny)',1,nx);

% Fitted yields
decomposition.expected = AP*ones(1,T) + BP*factors;

% Bond prices under Q with sigma = 0
aFull = nan(ny,1);
bFull = nan(ny,nx);
for n = 1:ny
    if n == 1
        aFull(n,1) = -delta0;
        bFull(n,:) = -delta1';
    else
        aFull(n,1) = aFull(n-1,1) - delta0 + bFull(n-1,:)*q0;
        bFull(n,:) = -delta1' + bFull(n-1,:)*qx;
    end
end

% Convert to yields
AQ = -1200*aFull./(1:ny)';
BQ = -1200*bFull./repmat((1:ny)',1,nx);

% Fitted yields
decomposition.riskPremium = (AQ-AP)*ones(1,T) + (BQ-BP)*factors;
decomposition.convexity   = (A-AQ)*ones(1,T);

end