% Calculate the pricing factors, which are principal components
function factors = compute_factors(yData,nx,T)

% De-mean the yields
meanYields = mean(yData,1);
yDataDemeaned = yData - repmat(meanYields,T,1);

% Compute the factors
[ySVD,~,~] = svd(yDataDemeaned'*yDataDemeaned);
pcLoadings = ySVD(:,1:nx);
factors = yDataDemeaned*pcLoadings;

% Standardise the factors.
[factors,~,sigmaX] = zscore(factors);
pcLoadings = pcLoadings./repmat(sigmaX,size(yData,2),1);

% Enforce average positive loadings.
negLoadings = mean(pcLoadings) < 0;
factors(:,negLoadings) = -factors(:,negLoadings);
pcLoadings(:,negLoadings) = -pcLoadings(:,negLoadings);
factors = factors';

end

