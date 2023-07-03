function RSS = yield_fit(deltas,ny,nx,T,factors,lambda0,lambda1,h0,hx,sigma,sSq,yData)

% Extract the short rate parameters
delta0 = deltas(1);
delta1 = deltas(2:nx+1);

% Get the fitted prices
decomposition = bond_pricing(ny,nx,T,factors,delta0,delta1,lambda0,lambda1,h0,hx,sigma,sSq);
yHat = decomposition.yHat;

% Get the sum of squared residuals
RSS = 0;
for t = 1:T
    selectMat = ~isnan(yData(:,t));
    selectMat(1) = 0; % We omit the short rate itself from the fitting
    RSS = RSS + (yHat(selectMat,t)-yData(selectMat,t))'*(yHat(selectMat,t)-yData(selectMat,t))/sum(selectMat)/T;
end

end