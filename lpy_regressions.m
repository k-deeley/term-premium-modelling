function [beta1,beta2] = lpy_regressions(parameters,decomposition,T,nx,ny,nMax,yData,factors,dataType)

switch dataType
    
    case 'actual'

        % Get the fitted yields
        yHat = decomposition.yHat;
        
        % Compute the term premium using the Dai and Singleton (2002)
        % definition
        termPremium = yHat - decomposition.expected;
        
        % Compute forward rates and the forward premium
        fHat = nan(nMax,T);
        fPremium = nan(nMax,T);
        for t = 1:T
            for n = 1:nMax % Note that fHat(i,t) = f(i-1,t) in Dai and Singleton notation
                if n == 1
                    fHat(n,t) = yHat(n,t);
                    fPremium(n,t) = termPremium(n,t);
                else
                    fHat(n,t) = n*yHat(n,t) - (n-1)*yHat(n-1,t);
                    fPremium(n,t) = n*termPremium(n,t) - (n-1)*termPremium(n-1,t);
                end
            end
        end

        % LPY regressions
        beta1 = nan(2,nMax-1);
        beta2 = nan(2,nMax-1);
        for n = 2:nMax
            LPY1 = (yData(n-1,2:end) - yData(n,1:end-1))';
            LPY2 = (yData(n-1,2:end) - yData(n,1:end-1) - (termPremium(n-1,2:end) - termPremium(n-1,1:end-1)) + fPremium(n,1:end-1)/(n-1))';
            X = [ones(T-1,1),(yData(n,1:end-1)-yData(1,1:end-1))'/(n-1)];
            beta1(:,n-1) = inv(X'*X)*X'*LPY1;
            beta2(:,n-1) = inv(X'*X)*X'*LPY2;
        end        
        
    case 'simulated'
        
        % Set number of simulations
        numSim = 100000;
        
        % Un-wrap parameters
        h0 = parameters.h0;
        hx = parameters.hx;
        sigma = parameters.sigma;
        sSq = parameters.sSq;
        delta0 = parameters.delta0;
        delta1 = parameters.delta1;
        lambda0 = parameters.lambda0;
        lambda1 = parameters.lambda1;
        
        % Simulate some factors
        factorsSim = nan(nx,numSim);
        factorsSim(:,1) = inv(eye(nx)-hx)*h0;
        cholSigma = chol(sigma,'lower');
        for i = 2:numSim
            factorsSim(:,i) = h0 + hx*factorsSim(:,i-1) + cholSigma*randn(nx,1);
        end
        
        % Fitted yields
        decomposition = bond_pricing(ny,nx,numSim,factorsSim,delta0,delta1,lambda0,lambda1,h0,hx,sigma,sSq);

        % Get the fitted yields
        yHat = decomposition.yHat;
        
        % Compute the term premium using the Dai and Singleton (2002)
        % definition
        termPremium = yHat - decomposition.expected;
        
        % Compute forward rates and the forward premium
        fHat = nan(nMax,numSim);
        fPremium = nan(nMax,numSim);
        for t = 1:numSim
            for n = 1:nMax % Note that fHat(i,t) = f(i-1,t) in Dai and Singleton notation
                if n == 1
                    fHat(n,t) = yHat(n,t);
                    fPremium(n,t) = termPremium(n,t);
                else
                    fHat(n,t) = n*yHat(n,t) - (n-1)*yHat(n-1,t);
                    fPremium(n,t) = n*termPremium(n,t) - (n-1)*termPremium(n-1,t);
                end
            end
        end
        
        % Regressions
        beta1 = nan(2,nMax-1);
        beta2 = nan(2,nMax-1);
        for n = 2:nMax
            LPY1 = (yHat(n-1,2:end) - yHat(n,1:end-1))';
            LPY2 = (yHat(n-1,2:end) - yHat(n,1:end-1) - (termPremium(n-1,2:end) - termPremium(n-1,1:end-1)) + fPremium(n,1:end-1)/(n-1))';
            X = [ones(numSim-1,1),(yHat(n,1:end-1)-yHat(1,1:end-1))'/(n-1)];
            beta1(:,n-1) = inv(X'*X)*X'*LPY1;
            beta2(:,n-1) = inv(X'*X)*X'*LPY2;
        end
end
        
end



