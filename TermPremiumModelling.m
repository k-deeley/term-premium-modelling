%% ACM Model Estimation.

%% Evaluate the model on the original data.
load("rawData75.mat")
dates = datetime(rawData(2:end, 1), "ConvertFrom", "datenum" );
yieldMatrix = rawData(2:end,2:end);
startDate = datetime('31-Jan-1985');
yieldMatrix = yieldMatrix(dates >= startDate, :);   
dates       = dates(dates >= startDate);            
yields      = array2timetable(yieldMatrix, "RowTimes", dates);
stir        = timetable(yields{:,1}/12, 'RowTimes', yields.Time);
yields = fillmissing(yields, "nearest");
nMax    = 120;              % Maximum maturity (in months)
matsAll = 1:nMax;           % Vector of all maturities up to nMax
decomposition = fitACM( yields, stir, matsAll );

idx = width(yieldMatrix);
figure
plot(dates,yieldMatrix(:, idx),'k')
hold on
plot(dates,decomposition.Fitted{:, idx},'r')
plot(dates,decomposition.RiskNeutralExpected{:, idx},'b')
plot(dates,decomposition.TermPremium{:, idx},'g')
plot(dates,decomposition.Convexity{:, idx},'c')
legend('Actual','Fitted','Expected','Premium','Convexity')
ylabel('(%)')
title('10-year decomposition')

%% Import the zero coupon yields.
filename = "ZeroCouponYieldsMonthlyWithBankRate.xlsx";
%% Select maturities up to 10 years
maturities = readmatrix( filename, "Range", "B1:H1" );
dates = datetime( readmatrix( filename, "Range", "A2:A330" ), ...
    "ConvertFrom", "excel" );
bankRates = readmatrix(filename, "Range", "N2:N330");
yields = readmatrix( filename, "Range", "B2:H330" );
varNames = "Maturity_" + maturities + "_years";
yields = array2timetable( yields, "RowTimes", dates, ...
    "VariableNames", varNames );

%% Select data from 1998 onwards to avoid missing values.
from1998Idx = yields.Time >= datetime( 2010, 1, 1 );
yields = yields(from1998Idx, :);
bankRates = bankRates(from1998Idx);

%% Interpolate the remaining short-term rates.
yields = fillmissing( yields, "linear" );

%% For each maturity, plot the corresponding yields over time.
figure
ax = axes;
plot( ax, yields.Properties.RowTimes, yields.Variables )
xlabel( ax, "Date" )
ylabel( ax, "Yield (%)" )
title( ax, "Zero-Coupon UK Yields (%)" )
grid( ax, "on" )
leg = legend( ax, string( maturities ), "NumColumns", 2 );
leg.Title.String = "Maturity (years)";
numMaturities = length( maturities );
ax.ColorOrder = jet( numMaturities );

%% Fit the ACM model.
% Extract the short-term interest rate.
stir = yields(:, 1);
stir{:, 1} = stir{:, 1} / 6;
% Convert the maturities to months.
yr2mth = 12;
maturitiesMonths = yr2mth * maturities;
% Fit the model.
decomposition = fitACM( yields, stir, maturitiesMonths );

%% Visualize the results.
maturityIdx = 7;

figure
plot( yields.Time, yields{:, maturityIdx}, ...
    "DisplayName", "Observed Yield" )
hold on
plot( decomposition.Fitted.Time, decomposition.Fitted{:, maturityIdx}, ...
    "DisplayName", "Fitted" )
plot( decomposition.RiskNeutralExpected.Time, ...
    decomposition.RiskNeutralExpected{:, maturityIdx}, ...
    "DisplayName", "Risk Neutral Expected" )
plot( decomposition.TermPremium.Time, ...
    decomposition.TermPremium{:, maturityIdx}, ...
    "DisplayName", "Term Premium" )
plot( decomposition.Convexity.Time, ...
    decomposition.Convexity{:, maturityIdx}, ...
    "DisplayName", "Convexity" )
legend


%% Interpolate yields for continuous range n = 1 to 120 months
% Input data is for years 0.5, 1, 2, 3, 5, 7 10
maxMaturity = 10 * 12;
maturitiesRequired = 1:maxMaturity;

% Bank rate is used for n=0
maturitiesAvailable = [0, maturities*12];
yieldsAvailable = array2table([bankRates yields.Variables]); % Convert to table so we can use rowfun
interpolateRow = @(rowyields) interp1(maturitiesAvailable, rowyields, maturitiesRequired);
% Easiest way I could find to stack the outputs of the rowfun is to ask for
% cells and then convert to a matrix; would welcome a better approach!
yieldsAll = cell2mat(rowfun(interpolateRow, yieldsAvailable, 'SeparateInputs', false, "OutputFormat", "cell"));
yieldsAll = array2timetable(yieldsAll, "RowTimes", yields.Time);

% Check that our interpolated values make sense
% Choose 12 random dates to plot
idxs = randi(height(yieldsAvailable), [1,12]);
figure
ax = axes;
plot(ax,maturitiesRequired, yieldsAll{idxs,:}, "DisplayName","Interpolated");
hold on;
plot(ax,maturitiesAvailable, yieldsAvailable{idxs,:}, 'O', "LineStyle", "none", "DisplayName","Observed");
ax.ColorOrder = jet( 12 );


%% Fit the ACM model, use the same stir as before
%stir = yieldsAll(:,1); % This gives quite poor results; approximating with 6
%month yields/6 as before works much better.
decomposition = fitACM( yieldsAll, stir, maturitiesRequired );

%% Visualize the results.
maturityYear = 5;
maturityIdx = yr2mth *maturityYear;

figure
plot( yieldsAll.Time, yieldsAll{:, maturityIdx}, ...
    "DisplayName", "Observed Yield" )
 hold on
 plot( decomposition.Fitted.Time, decomposition.Fitted{:, maturityIdx}, ...
     "DisplayName", "Fitted" )
 plot( decomposition.RiskNeutralExpected.Time, ...
     decomposition.RiskNeutralExpected{:, maturityIdx}, ...
     "DisplayName", "Risk Neutral Expected" )
 plot( decomposition.TermPremium.Time, ...
     decomposition.TermPremium{:, maturityIdx}, ...
     "DisplayName", "Term Premium" )
plot( decomposition.Convexity.Time, ...
    decomposition.Convexity{:, maturityIdx}, ...
    "DisplayName", "Convexity" )
legend




