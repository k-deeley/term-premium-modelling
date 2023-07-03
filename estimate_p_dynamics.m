function [h0,hx,sigma,V] = estimate_p_dynamics(nx,T,factors,model)

switch model

    % In the fully flexible case, we allow the constant term to be
    % non-zero.
    case 'benchmark'
   
        % Create regressors
        X = [ones(1,T-1);factors(:,1:T-1)];
        Y = factors(:,2:T);

        % Run the regressions
        hTmp = (X'\Y')';
        h0 = hTmp(:,1);
        hx = hTmp(:,2:end);
        V = Y - hTmp*X;
        sigma = cov(V');
        
	% In the benchmark case, we set h0 = 0 to match the sample mean (recall
	% that the factors are normalised to have mean zero).
    case 'meanzero'

        % Create regressors
        X = factors(:,1:T-1);
        Y = factors(:,2:T);

        % Run the regressions
        h0 = zeros(nx,1);
        hx = (X'\Y')';
        V = Y - hx*X;
        sigma = cov(V');
        
    case 'biascorrect' % Follows Bauer, Rudebusch and Wu (2012), Appendix A
        
        % Create regressors
        X = factors(:,1:T-1);
        Y = factors(:,2:T);

        % Run the regressions
        h0 = zeros(nx,1);
        hx = (X'\Y')';
        thetaHat = hx(:);
        V = Y - hx*X;
        
        % Bootstrap
        numMC = 1000;
        thetaMCArray = nan(nx^2,numMC);
        for j = 1:numMC
            % Generating the data for bootstrap
            indexMC = randi(T-1,1,T);
            factorsMC(:,1) = factors(:,indexMC(1));
            for i = 2:T
                 factorsMC(:,i:i) = h0 + hx*factorsMC(:,i-1) + V(:,indexMC(i));
            end

            % Using the estimator on the bootstrap sample
            XMC = factorsMC(:,1:T-1);
            YMC = factorsMC(:,2:T);
            hxMC = (XMC'\YMC')';
            thetaMCArray(:,j) = hxMC(:);
        end
        
        % Compute the bias
        thetaMC = mean(thetaMCArray,2);
        bias = thetaMC - thetaHat;

        % Apply the bias correction
        thetaBC = thetaHat - bias;

        % Enforce stationarity by shrinking towards the OLS estimate
        hx = reshape(thetaBC,nx,nx);
        scaling = 1;
        while any(abs(eig(hx)) >= 1)
            % Re-scale
            scaling = scaling*0.999;
            thetaBC = thetaHat - scaling*bias;

            % Re-compute hx
            hx = reshape(thetaBC,nx,nx);
        end

        % Re-compute the covariance
        V = Y - hx*X;
        sigma = cov(V');

    case 'inverse' % - Follows Bauer, Rudebusch and Wu (2012), Appendix B

        % Create regressors
        X = factors(:,1:T-1);
        Y = factors(:,2:T);

        % Run the regressions
        h0 = zeros(nx,1);
        hx = (X'\Y')';
        thetaHat = hx(:);
        V = Y - hx*X;

        % Set options - as in Bauer, Rudebusch and Wu (2012)
        alpha = 0.2;
        numMC = 50;
        chainLen = 1000;
        burnIn = 500;
        
        % Inverse bootstrap
        thetaChain = nan(nx^2,chainLen);
        theta = thetaHat;
        for k = 1:chainLen
            % Get the bias at this point in the chain
            thetaMCArray = nan(nx^2,numMC);
            for j = 1:numMC
                indexMC = randi(T-1,1,T);
                factorsMC(:,1) = factors(:,indexMC(1));
                for i = 2:T
                    factorsMC(:,i:i) = h0 + hx*factorsMC(:,i-1) + V(:,indexMC(i));
                end
                XMC = factorsMC(:,1:T-1);
                YMC = factorsMC(:,2:T);
                hxMC = (XMC'\YMC')';
                thetaMCArray(:,j) = hxMC(:);
            end
            thetaMC = mean(thetaMCArray,2);
            bias = thetaHat - thetaMC;

            % Construct the next point in the chain
            theta = theta + alpha*bias;
            hx = reshape(theta,nx,nx);
            thetaChain(:,k) = theta;
        end
        thetaBC = mean(thetaChain(:,burnIn+1:end),2);
        hx = reshape(thetaBC,nx,nx);        
        
        % Shrink hx towards zero to ensure stationarity
        while any(abs(eig(hx)) >= 1)
            hx = hx*0.999;
        end
        
        % Estimate covariance
        V = Y - hx*X;
        sigma = cov(V');  
end
        
end