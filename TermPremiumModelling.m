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
stir = array2timetable( yields{:, 1}, "RowTimes", yields.Time );

% Interpolate missing values.
yields = fillmissing( yields, "nearest" );

% Fit the model to obtain the decomposition.
maxMaturity = 120; % Maximum maturity (in months)
maturities = 1 : maxMaturity; % Vector of all maturities
excessReturnMaturities = 6 : 6 : maxMaturity;
decomposition = fitACM( yields, stir, maturities, "excessReturnMaturities", excessReturnMaturities);

% %% Visualize the model results for selected maturities
year2month = 12;
selectedMaturities = [24, 60, 120];
selectedYears = selectedMaturities / year2month;

f = figure;
f.Position(3:4) = [1200,350];
tiledlayout(1, length(selectedMaturities))
for k = 1 : length( selectedMaturities )
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
filename = "ZeroCouponYieldsMonthly_Dummy.xlsx"; %ZeroCouponYieldsMonthly_1992_2023

% Select maturities up to 10 years.
maturities = readmatrix( filename, "Range", "B1:DQ1" );
datesNew = datetime( readmatrix( filename, "Range", "A2:A371" ), ...
    "ConvertFrom", "excel" );
% bankRates = readmatrix(filename, "Range", "R2:R330");
yieldsNew = readmatrix( filename, "Range", "B2:DQ371" );
varNames = "Maturity_" + maturities + "_months";
yieldsNew = array2timetable( yieldsNew, "RowTimes", datesNew, ...
    "VariableNames", varNames );

% Interpolate the remaining short-term rates.
% Interpolate top-down, using adjacent maturities
yieldsNew = fillmissing( yieldsNew, "linear");

%% For each maturity, plot the corresponding yields over time.
figure
ax = axes;
plot( ax, yieldsNew.Properties.RowTimes, yieldsNew{:,6:end} )
xlabel( ax, "Date" )
ylabel( ax, "Yield (%)" )
title( ax, "Zero-Coupon Yields (%)" )
grid( ax, "on" )
numMaturities = length( maturities );
ax.ColorOrder = jet( numMaturities );

%% Fit the ACM model.
% Extract the short-term interest rate.
stir = yieldsNew(:, 1);
stir{:,1} = stir{:,1};

selectedMaturities = [24, 60, 120];
selectedYears = selectedMaturities / year2month;

decomposition = fitACM( yieldsNew, stir, maturities, "excessReturnMaturities", excessReturnMaturities);
colors = {[139 0 0]/255, [0 128 0]/255, 'b'};
%% Visualize the results.
for k = 1 : length( selectedYears )
    figure
    plot( yieldsNew.Time, yieldsNew{:, selectedMaturities(k)}, "k", ...
        "DisplayName", "Observed Yield" )
    hold on
    plot( decomposition.Fitted.Time, ...
        decomposition.Fitted{:, selectedMaturities(k)}, ...
        "Color", [139 0 0]/255, ...
        "DisplayName", "Fitted Yield" )
    plot( decomposition.RiskNeutralExpected.Time, ...
        decomposition.RiskNeutralExpected{:, selectedMaturities(k)}, ...
        "Color", [0 128 0]/255, ...
        "DisplayName", "Risk Neutral Yield (Expected)" )
    plot( decomposition.TermPremium.Time, ...
        decomposition.TermPremium{:, selectedMaturities(k)}, "b", ...
        "DisplayName", "Term Premium" )
    legend
    grid on
    title( selectedYears(k) + "-year Maturity Decomposition " )
end % for

%% Visualize Term Premium for all selectedMaturities
figure
for k = 1 : length( selectedYears )
    p = plot(decomposition.TermPremium.Time, ...
            decomposition.TermPremium{:, selectedMaturities(k)}, "Color", colors{k}, ...
            "DisplayName", selectedYears(k) +"-year" );
    hold on;
end % for
legend
grid on
title( "Term premium for selected maturities " )
hold off;

%% Compute the RMSE between the fitted and the observed yields.
RMSE = 100 * sqrt(mean((decomposition.Fitted.Variables - yieldsNew.Variables).^2));

RMSE = array2table(RMSE, "VariableNames", yieldsNew.Properties.VariableNames, ...
    "RowNames", "RMSE");

selectedRMSE = RMSE(:, [24, 60, 120]);

