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
excessReturnMaturities = 6 : 6 : maxMaturity;
decomposition = fitACM( yields, stir, allMaturities, excessReturnMaturities);

% %% Visualize the model results.
year2month = 12;
selectedMaturities = [24, 60, 120];
selectedYears = selectedMaturities / year2month;

selectedIdx = selectedMaturities / 6;

for k = 1 : length( selectedMaturities )
    figure
    nexttile
    plot( dates, yieldMatrix(:, selectedMaturities(k)), "k" )
    hold on
    plot( dates, decomposition.Fitted{:, selectedIdx(k)}, "r" )
    plot( dates, decomposition.RiskNeutralExpected{:, selectedIdx(k)}, "b" )
    plot( dates, decomposition.TermPremium{:, selectedIdx(k)}, "g" )
    legend( "Actual", "Fitted", "Expected", "Premium" )
    ylabel( "(%)" )
    title( selectedYears(k) + "-year Decomposition (GSW Data)" )
    grid on
end

%% Move to the BoE data - start by importing the zero coupon yields.
filename = "ZeroCouponYieldsMonthly_1992_2023.xlsx";

% Select maturities up to 10 years.
maturities = readmatrix( filename, "Range", "B1:DQ1" );
datesNew = datetime( readmatrix( filename, "Range", "A2:A371" ), ...
    "ConvertFrom", "excel" );
% bankRates = readmatrix(filename, "Range", "R2:R330");
yieldsNew = readmatrix( filename, "Range", "B2:DQ371" );
varNames = "Maturity_" + maturities + "_months";
yieldsNew = array2timetable( yieldsNew, "RowTimes", datesNew, ...
    "VariableNames", varNames );

% Select data from January 1996 onwards to avoid missing values.
% from1998Idx = yieldsNew.Time >= datetime( 1996, 1, 1 );
%yieldsNew = yieldsNew(from1998Idx, :);
% bankRates = yieldsNew(:,1);

% Interpolate the remaining short-term rates.
% Interpolate top-down, using adjacent maturities
yieldsNew = fillmissing( yieldsNew, "linear");

%% For each maturity, plot the corresponding yields over time.
figure
ax = axes;
plot( ax, yieldsNew.Properties.RowTimes, yieldsNew{:,6:end} )
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
stir{:,1} = stir{:,1}/12;
% Maturities are alreay in months; no need to convert
% maturities = maturities(12:end);

% Regimes: 
% start - 30/11/2008; 30/11/2008 - 31/12/2021; 31/12/2021 - end
regimeBounds = [yieldsNew.Time(1), datetime(2008, 11, 30), datetime(2021, 12, 31), yieldsNew.Time(end)];

%% Visualize the results.
selectedMaturities = [24, 60, 120];
selectedYears = selectedMaturities / year2month;
selectedIdx = selectedMaturities / 6;

for r = 1:numel(regimeBounds)-1
	regimeIdx = timerange(regimeBounds(r), regimeBounds(r+1), "intervalType", "open");
    regimeYields = yieldsNew(regimeIdx,:);
    regimeStir = stir(regimeIdx,:);
    % Fit the model to regime
    decomposition = fitACM( regimeYields, regimeStir, maturities, excessReturnMaturities);

    for k = 1 : length( selectedYears )

        figure
		% yieldsNew is indexed for all n. So we should use indices like 24, 60, 120.
        plot( yieldsNew.Time, yieldsNew{:, selectedMaturities(k)}, "k", ...
            "DisplayName", "Observed Yield" )
        hold on
        plot( decomposition.Fitted.Time, ...
            decomposition.Fitted{:, selectedIdx(k)}, "r", ...
            "DisplayName", "Fitted Yield" )
        plot( decomposition.RiskNeutralExpected.Time, ...
            decomposition.RiskNeutralExpected{:, selectedIdx(k)}, "b", ...
            "DisplayName", "Expectation component" )
        plot( decomposition.TermPremium.Time, ...
            decomposition.TermPremium{:, selectedIdx(k)}, "g", ...
            "DisplayName", "Term Premium" )
        legend
        grid on
        title( selectedYears(k) + "-year Maturity Decomposition (BoE Data)" )
    end % for
end % for
