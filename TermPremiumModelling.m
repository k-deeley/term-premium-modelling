%% ACM Model Estimation.

%% Evaluate the ACM model on the GSW [2007] data set.

% Load the data.
dataStruct = load( "rawData75.mat" );
rawData = dataStruct.rawData;

% Separate the dates and yields.
dates = datetime( rawData(2:end, 1), "ConvertFrom", "datenum" );
yieldMatrix = rawData(2:end, 2:end);

% Extract data from a chosen point onwards.
startDate = datetime( 1985, 1, 31 );
yieldMatrix = yieldMatrix(dates >= startDate, :);
dates = dates(dates >= startDate);

% Convert the yields to a timetable.
yields = array2timetable( yieldMatrix, "RowTimes", dates );
stir = array2timetable( yields{:, 1} / 12, "RowTimes", yields.Time );

% Interpolate missing values.
yields = fillmissing( yields, "nearest" );

% Fit the model to obtain the decomposition.
maxMaturity = 120; % Maximum maturity (in months)
allMaturities = 1 : maxMaturity; % Vector of all maturities
decomposition = fitACM( yields, stir, allMaturities );

%% Visualize the model results.
year2month = 12;
selectedMaturities = [24, 60, 120];
selectedYears = selectedMaturities / year2month;

for k = 1 : length( selectedMaturities )
    figure
    nexttile
    plot( dates, yieldMatrix(:, selectedMaturities(k)), "k" )
    hold on
    plot( dates, decomposition.Fitted{:, selectedMaturities(k)}, "r" )
    plot( dates, decomposition.RiskNeutralExpected{:, selectedMaturities(k)}, "b" )
    plot( dates, decomposition.TermPremium{:, selectedMaturities(k)}, "g" )
    legend( "Actual", "Fitted", "Expected", "Premium" )
    ylabel( "(%)" )
    title( selectedYears(k) + "-year Decomposition (GSW Data)" )
    grid on
end

%% Move to the BoE data - start by importing the zero coupon yields.
filename = "ZeroCouponYieldsMonthlyWithBankRate.xlsx";

% Select maturities up to 10 years.
maturities = readmatrix( filename, "Range", "B1:H1" );
datesNew = datetime( readmatrix( filename, "Range", "A2:A330" ), ...
    "ConvertFrom", "excel" );
bankRates = readmatrix(filename, "Range", "N2:N330");
yieldsNew = readmatrix( filename, "Range", "B2:H330" );
varNames = "Maturity_" + maturities + "_years";
yieldsNew = array2timetable( yieldsNew, "RowTimes", datesNew, ...
    "VariableNames", varNames );

% Select data from January 1996 onwards to avoid missing values.
from1998Idx = yieldsNew.Time >= datetime( 1996, 1, 1 );
yieldsNew = yieldsNew(from1998Idx, :);
bankRates = bankRates(from1998Idx);

% Interpolate the remaining short-term rates.
yieldsNew = fillmissing( yieldsNew, "linear" );

%% For each maturity, plot the corresponding yields over time.
figure
ax = axes;
plot( ax, yieldsNew.Properties.RowTimes, yieldsNew.Variables )
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
stir = yieldsNew(:, 1);
stir{:, 1} = stir{:, 1} / 6;
% Convert the maturities to months.
year2month = 12;
maturitiesMonths = year2month * maturities;
% Fit the model.
decomposition = fitACM( yieldsNew, stir, maturitiesMonths );

%% Visualize the results.
[~, selectedIdx] = ismember( selectedYears, maturities );

for k = 1 : length( selectedMaturities )
    
    figure
    plot( yieldsNew.Time, yieldsNew{:, selectedIdx(k)}, "k", ...
        "DisplayName", "Observed Yield" )
    hold on
    plot( decomposition.Fitted.Time, ...
        decomposition.Fitted{:, selectedIdx(k)}, "r", ...
        "DisplayName", "Fitted" )
    plot( decomposition.RiskNeutralExpected.Time, ...
        decomposition.RiskNeutralExpected{:, selectedIdx(k)}, "b", ...
        "DisplayName", "Risk Neutral Expected" )
    plot( decomposition.TermPremium.Time, ...
        decomposition.TermPremium{:, selectedIdx(k)}, "g", ...
        "DisplayName", "Term Premium" )
    legend
    grid on
    title( selectedYears(k) + "-year Maturity Decomposition (BoE Data)" )
end % for

% %% Interpolate yields for continuous range n = 1 to 120 months
%
% % Input data is for years 0.5, 1, 2, 3, 5, 7, 10.
% maxMaturity = 10 * 12;
% allMaturities = 1 : maxMaturity;
%
% % Convert the maturities to months.
% maturitiesAvailable = [0, maturities * yr2mth];
%
% % The bank rate is used for n = 0. Convert the matrix to a table so we can
% % use rowfun.
% yieldsAvailable = array2table( [bankRates, yieldsNew.Variables] );
% interpolateRow = @(rowYields) interp1( maturitiesAvailable, rowYields, allMaturities, "makima" );
% % Easiest way I could find to stack the outputs of the rowfun is to ask for
% % cells and then convert to a matrix; would welcome a better approach!
% yieldsAll = cell2mat( rowfun( interpolateRow, yieldsAvailable, "SeparateInputs", false, "OutputFormat", "cell" ) );
% yieldsAll = array2timetable( yieldsAll, "RowTimes", yieldsNew.Time );
%
% %% Fit the ACM model, use the same stir as before
% %stir = yieldsAll(:,1); % This gives quite poor results; approximating with 6
% %month yields/6 as before works much better.
% stir = yieldsAll(:, 6);
% stir{:, 1} = stir{:, 1} / 6;
% decomposition = fitACM( yieldsAll, stir, allMaturities );
%
% %% Visualize the results.
% maturityYear = 10;
% maturityIdx = yr2mth * maturityYear;
%
% figure
% plot( yieldsAll.Time, yieldsAll{:, maturityIdx}, ...
%     "DisplayName", "Observed Yield" )
% hold on
% plot( decomposition.Fitted.Time, decomposition.Fitted{:, maturityIdx}, ...
%     "DisplayName", "Fitted" )
% plot( decomposition.RiskNeutralExpected.Time, ...
%     decomposition.RiskNeutralExpected{:, maturityIdx}, ...
%     "DisplayName", "Risk Neutral Expected" )
% plot( decomposition.TermPremium.Time, ...
%     decomposition.TermPremium{:, maturityIdx}, ...
%     "DisplayName", "Term Premium" )
% plot( decomposition.Convexity.Time, ...
%     decomposition.Convexity{:, maturityIdx}, ...
%     "DisplayName", "Convexity" )
% legend