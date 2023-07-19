%************************************************************************************************
% Options
%************************************************************************************************
clear
close all
clc
warning off

% Model options
nx = 5;                             % Number of factors
model = 'meanzero';                 % Type of estimator for P dynamics
                                    %       (1) 'benchmark' -> h0/hx
                                    %       estimated by OLS
                                    %       (2) 'meanzero' -> h0 set to
                                    %       match sample mean (as in ACM)
                                    %       (3) 'biascorrect' -> as (2)
                                    %       with bootstrap bias correction
                                    %       (4) 'inverse' -> as (2) with
                                    %       inverse bias correction of
                                    %       Bauer, Rudebusch and Wu (2012)
optimiseDeltas = 0;                 % Choose short rate parameters to 
                                    % maximise fit to entire yield curve 
                                    % (1) or just short rate (0)
restrictPOR = 0;                    % Impose restrictions on the price of 
                                    % risk (1) or not (0)
startDate = datenum('31-Jan-1985'); % Start date

%************************************************************************************************
% Data
%************************************************************************************************
% Specify maturity vectors
nMax = 120;                         % Maximum maturity (in months)
matsAll = 1:nMax;                   % Vector of all maturities up to nMax
pcMats = 12:120;                    % Maturities used to compute the pricing factors
rxMats = 18:6:120;                  % Maturities used for the return regressions
nrx = size(rxMats,2);               % Number of excess returns in the return regressions

% Read in raw data
load('rawData75.mat');              % Uses VRP data for all except 1-month,
                                    % which is interpolated between policy 
                                    % rate and 3-month T-bill rate
dates = rawData(2:end,1);           % Extract dates
yData = rawData(2:end,2:end)';      % Extract all the yields
yData = yData(:,dates>=startDate);  % Remove any data before the start date
dates = dates(dates>=startDate);    % Remove any dates before the start date
[ny,T] = size(yData);               % Number of time periods and yields

% Compute the pricing factors - uses only the maturities in pcMats
factors = compute_factors(yData(pcMats,:)',nx,T);

%************************************************************************************************
% Estimation
%************************************************************************************************
% Estimate the model
[parameters,decomposition,returns] = acm_estimation(nx,ny,nrx,T,factors,yData,matsAll,rxMats,model,optimiseDeltas);

% Covariance matrix
parameters.vTheta = compute_covariance_matrix(parameters,factors,nx,nrx,T);

% Restrict price of risk
if restrictPOR == 1
    
    % Save unadjusted results
    parametersUnadj = parameters;
    returnsUnadj = returns;
    
    % Impose the restrictions
    [parameters,decomposition,returns] = restrict_por(nx,ny,T,nrx,factors,model,parameters,returns)
    
end

%************************************************************************************************
% LPY regressions
%************************************************************************************************
% Dai and Singleton LPY regressions
% [beta1Data,beta2Data] = lpy_regressions(parameters,decomposition,T,nx,ny,nMax,yData,factors,'actual');
% [beta1Sim,beta2Sim] = lpy_regressions(parameters,decomposition,T,nx,ny,nMax,yData,factors,'simulated');

%************************************************************************************************
% Plot results
%************************************************************************************************
% Plot 3-month decomposition
figure
plot(dates,yData(3,:),'k')
hold on
plot(dates,decomposition.yHat(3,:),'r')
plot(dates,decomposition.expected(3,:),'b')
plot(dates,decomposition.riskPremium(3,:),'g')
plot(dates,decomposition.convexity(3,:),'c')
legend('Actual','Fitted','Expected','Premium','Convexity')
ylabel('Per cent')
datetick('x','yyyy')
title('3-month decomposition')

% Plot 10-year decomposition
figure
plot(dates,yData(end,:),'k')
hold on
plot(dates,decomposition.yHat(end,:),'r')
plot(dates,decomposition.expected(end,:),'b')
plot(dates,decomposition.riskPremium(end,:),'g')
plot(dates,decomposition.convexity(end,:),'c')
legend('Actual','Fitted','Expected','Premium','Convexity')
ylabel('Per cent')
datetick('x','yyyy')
title('10-year decomposition')

% Plot return fit
index = 1;
for i = 1:nrx
    % Start a new figure if old one is full
    if index == 1
        figure
    end
    
    % Plot
    subplot(3,1,index)
    hold on
    plot(dates(1:T-1,1),returns.rx(i,:),'k')
    plot(dates(1:T-1,1),returns.rxHat(i,:),'r')
    datetick('x','yyyy')
    
	% Keep track of how many sub-plots in the current figure
    if index == 3
        index = 1;
    else
        index = index + 1;
    end
end
    
% Plot LPY results
% figure
% plot(beta1Data(2,:),'k')
% hold on
% plot(beta1Sim(2,:),'--k')
% plot(beta2Data(2,:),':k')
% legend('Actual data','Simulated data - LPY(i)','Premia-adjusted - LPY(ii)')


