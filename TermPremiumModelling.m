%% ACM Model Estimation.

%% Import the zero coupon yields.
filename = "ZeroCouponYieldsMonthly.xlsx";
maturities = readmatrix( filename, "Range", "B1:K1" );
dates = datetime( readmatrix( filename, "Range", "A2:A330" ), ...
    "ConvertFrom", "excel" );
yields = readmatrix( filename, "Range", "B2:K330" );
varNames = "Maturity_" + maturities + "_years";
T = array2timetable( yields, "RowTimes", dates, ...
    "VariableNames", varNames );

%% Remove rows with missing  data.
missingIdx = ismissing( T );
badRows = any( missingIdx, 2 );
T(badRows, :) = [];
numObservations = height( T );

%% For each maturity, plot the corresponding yields over time.
figure
ax = axes;
plot( ax, T.Properties.RowTimes, T.Variables )
xlabel( ax, "Date" )
ylabel( ax, "Yield (%)" )
title( ax, "Zero-Coupon UK Yields" )
grid( ax, "on" )
leg = legend( ax, string( maturities ), "NumColumns", 2 );
leg.Title.String = "Maturity (years)";
numMaturities = length( maturities );
ax.ColorOrder = jet( numMaturities );

%%
yields = T.Variables;
numFactors = 4;
decomposition = estimateACM( yields, maturities, numFactors );

idx = 7; % 10-year maturity

figure
hold on
plot( T.Time, yields(:, idx), "LineWidth", 2, "DisplayName", "Observed Yield (%)" )
plot( T.Time, decomposition.TermPremium(:, idx), "LineWidth", 2, "DisplayName", "Term Premium" )
plot( T.Time, decomposition.Expected(:, idx), "LineWidth", 2, "DisplayName", "Expected" )
plot( T.Time, decomposition.Fitted(:, idx), "LineWidth", 2, "DisplayName", "Fitted" )
plot( T.Time, decomposition.Convexity(:, idx), "LineWidth", 2, "DisplayName", "Convexity" )
legend