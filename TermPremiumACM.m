%% ACM Model Estimation.

%% Import the data.
zcy = readtimetable( "zeroCurve.xlsx" );

%% Write down the required maturities.
maturities = str2double( extractAfter( ...
    zcy.Properties.VariableNames, "_" ) );
excessReturnMaturities = 6 : 6 : max( maturities );

year2month = 12;
selectedMaturities = [24, 60, 120];
selectedYears = selectedMaturities / year2month;

%% Extract the short-term interest rate (STIR).
stir = zcy(:, 1);

%% Fit the model.
decomposition = fitACM( zcy, stir, maturities, ...
    "excessReturnMaturities", excessReturnMaturities );

%% Plot the corresponding yields over time.
figure
ax = axes;
plot( ax, zcy.Time, zcy.Variables )
xlabel( ax, "Date" )
ylabel( ax, "Yield (%)" )
title( ax, "Zero-Coupon Yields (%)" )
grid( ax, "on" )

numMaturities = width( zcy );
ax.ColorOrder = jet( numMaturities );

%% Visualize the results.
for k = 1 : numel( selectedYears )
    figure
    plot( zcy.Time, zcy{:, selectedMaturities(k)}, "k", ...
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
    title( selectedYears(k) + "-year Maturity Decomposition" )
end % for

%% Visualize Term Premium for all selectedMaturities
colors = {[139 0 0]/255, [0 128 0]/255, 'b'};

figure
for k = 1 : numel( selectedYears )
    p = plot( decomposition.TermPremium.Time, ...
            decomposition.TermPremium{:, selectedMaturities(k)}, ...
            "Color", colors{k}, ...
            "DisplayName", selectedYears(k) + "-year" );
    hold on
end % for
legend
grid on
title( "Term premium for selected maturities " )
hold off

%% Compute the RMSE between the fitted and the observed yields.
RMSE = 100 * rms( decomposition.Fitted.Variables - zcy.Variables );
RMSE = array2table( RMSE, ...
    "VariableNames", zcy.Properties.VariableNames, ...
    "RowNames", "RMSE" );
selectedRMSE = RMSE(:, [24, 60, 120]);