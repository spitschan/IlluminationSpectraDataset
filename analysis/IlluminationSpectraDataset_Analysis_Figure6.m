%% IlluminationSpectraDataset_Analysis_Figure5
%
% Generates Figure 6.
%
% 12/10/2015  spitschan   Wrote it.

fprintf('\n>> Making Figure 6...');

%% Record the daylight/twilight/night regime for each spectrum
for locIndex = [1 2];
    % Save in the data struct
    data(locIndex).regime = NaN*ones(size(data(locIndex).spd, 2), 1);
    data(locIndex).regime(data(locIndex).theIdxNight) = 1;                 % 1 = Night
    data(locIndex).regime(data(locIndex).theIdxTwilightAstronomical) = 2;  % 2 = Astronomical twilight
    data(locIndex).regime(data(locIndex).theIdxTwilightNautical) = 3;      % 3 = Nautical twilight
    data(locIndex).regime(data(locIndex).theIdxTwilightCivic) = 4;         % 4 = Civic twilight
    data(locIndex).regime(data(locIndex).theIdxDay) = 5;                   % 5 = Daylight
    data(locIndex).moonless = zeros(size(data(locIndex).spd, 2), 1);
    data(locIndex).moonless(data(locIndex).theIdx) = 1;                    % Logical value
end

%% Merge the data
data_merged.spd = [data(1).spd data(2).spd];
data_merged.locIndex = [1*ones(size(data(1).spd, 2), 1) ; 2*ones(size(data(2).spd, 2), 1)];
data_merged.regime = [data(1).regime ; data(2).regime];
data_merged.moonless = [data(1).moonless ; data(2).moonless];
data_merged.solarElevationAngleDeg = [data(1).solarElevationDeg data(2).solarElevationDeg];

%% Get the CIE model
B_CIE_360_830 = IlluminationSpectraDataset_Analysis_GetModel('CIE', (360:1:830)', true);

%% 1) Fit the CIE model
B = B_CIE_360_830;
[spds_norm, spds_fit, w, gof, residuals] = IlluminationSpectraDataset_Analysis_FitModel(wls, data_merged.spd, (360:1:830)', B, 360, 830);
data_merged.gof_CIE = gof; 
data_merged.resid_CIE = residuals;

% Extract mean residuals for daylight
MU_day = mean(data_merged.resid_CIE(:, (data_merged.regime == 5) & data_merged.moonless), 2);

%% 2) Fit the CIE model + mean of daylight residuals
B = [B_CIE_360_830 MU_day];
[spds_norm, spds_fit, w, gof, residuals] = IlluminationSpectraDataset_Analysis_FitModel(wls, data_merged.spd, (360:1:830)', B, 360, 830);
data_merged.gof_CIE_ResidDayMean = gof; 
data_merged.resid_CIE_ResidDayMean = residuals;

% Extract mean residuals for civic twilight
MU_civicTwilight = mean(data_merged.resid_CIE_ResidDayMean(:, (data_merged.regime == 4) & data_merged.moonless), 2);

%% 3) Fit the CIE model + mean of daylight residuals + mean of civic twilight residuals
B = [B_CIE_360_830 MU_day MU_civicTwilight];
[spds_norm, spds_fit, w, gof, residuals] = IlluminationSpectraDataset_Analysis_FitModel(wls, data_merged.spd, (360:1:830)', B, 360, 830);
data_merged.gof_CIE_ResidDayMean_ResidCivicTwilight = gof; 
data_merged.resid_CIE_ResidDayMean_ResidCivicTwilight = residuals;

% Extract mean residuals for atronomical twilight
MU_astronomicalTwilight1 = mean(data_merged.resid_CIE_ResidDayMean_ResidCivicTwilight(:, (data_merged.regime == 2) & (data_merged.moonless) & (data_merged.locIndex == 1)), 2);
MU_astronomicalTwilight2 = mean(data_merged.resid_CIE_ResidDayMean_ResidCivicTwilight(:, (data_merged.regime == 2) & (data_merged.moonless) & (data_merged.locIndex == 2)), 2);

%% Fit the CIE model + mean of daylight residuals
% Location 1
theIdx1 = find(data_merged.locIndex == 1);
B = [B_CIE_360_830 MU_day MU_civicTwilight MU_astronomicalTwilight1];
[spds_norm, spds_fit, w, gof, residuals] = IlluminationSpectraDataset_Analysis_FitModel(wls, data_merged.spd(:, theIdx1), (360:1:830)', B, 360, 830);
data_merged.gof_CIE_ResidDayMean_ResidCivicTwilight_ResidAstronomicalTwilight(theIdx1) = gof;
data_merged.resid_CIE_ResidDayMean_ResidCivicTwilight_ResidAstronomicalTwilight(:, theIdx1) = residuals;

% Location 2
theIdx2 = find(data_merged.locIndex == 1);
B = [B_CIE_360_830 MU_day MU_civicTwilight MU_astronomicalTwilight1];
[spds_norm, spds_fit, w, gof, residuals] = IlluminationSpectraDataset_Analysis_FitModel(wls, data_merged.spd(:, theIdx2), (360:1:830)', B, 360, 830);
data_merged.gof_CIE_ResidDayMean_ResidCivicTwilight_ResidAstronomicalTwilight(theIdx2) = gof;
data_merged.resid_CIE_ResidDayMean_ResidCivicTwilight_ResidAstronomicalTwilight(:, theIdx2) = residuals;

%% Plot basis functions
gcFig6 = figure;
wlRangeStartIdxModel = find(wls == 360);
wlRangeEndIdxModel = find(wls == 830);
validWlRngIdxModel = wlRangeStartIdxModel:wlRangeEndIdxModel;

% Panel A
B_CIE3 = [B_CIE_360_830 MU_day MU_civicTwilight MU_astronomicalTwilight1 MU_astronomicalTwilight2];
for i = 1:size(B_CIE3, 2)
    B_CIE3(:, i) = B_CIE3(:, i)/norm(B_CIE3(~isnan(B_CIE3(:, i)), i));
end
subplot(2, 4, 1);
plot(wls(validWlRngIdxModel), B_CIE3(:, 4), 'Color', [14 127 46]/255); hold on
plot(wls(validWlRngIdxModel), B_CIE3(:, 5), 'Color', [115 81 191]/255); hold on
xlabel('Wavelength [nm]'); ylabel('Relative value');
pbaspect([1 1 1]);
xlim([230 890]); pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off;
ylim([-0.25 0.25]);

subplot(2, 4, 2);
plot(wls(validWlRngIdxModel), B_CIE3(:, 6), 'Color', [0.85 0.22 0.16]); hold on;
xlabel('Wavelength [nm]'); ylabel('Relative value');
pbaspect([1 1 1]);
xlim([230 890]); pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off;
ylim([-0.25 0.25]);
plot([558 558], [0.22 0.22], '.r', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');

% Plot the night data from Warrant et al.
dataPath = fullfile(fileparts(mfilename('fullpath')), 'data');
tmp = csvread(fullfile(dataPath, 'spd_Zabriskie', 'data.csv') ,1);
[~, sortIndex] = sort(tmp(:, 1));
tmp = [tmp(sortIndex, 1) tmp(sortIndex, 2)];
wls_Warrant = (300:700)';
spd_Warrant = QuantaToEnergy(wls_Warrant, interp1(tmp(:, 1), tmp(:, 2), wls_Warrant));
plot(wls_Warrant, -0.2+0.08*(spd_Warrant/max(spd_Warrant)), '-k');
plot([280 280], [-0.2 -0.08], '-k');

subplot(2, 4, 6);
plot(wls(validWlRngIdxModel), B_CIE3(:, 7), 'Color', [0.26 0.49 0.76]); hold on;
xlabel('Wavelength [nm]'); ylabel('Relative value');
pbaspect([1 1 1]);
xlim([230 890]); pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off;
ylim([-0.25 0.25]);
plot([819 819], [0.22 0.22], '.r', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
plot([570 615], [0.22 0.22], '-r');

% Plot the light pollution data from Warrant et al.
dataPath = fullfile(fileparts(mfilename('fullpath')), 'data');
tmp = csvread(fullfile(dataPath, 'spd_LightPollution', 'data.csv') ,1);
[~, sortIndex] = sort(tmp(:, 1));
tmp = [tmp(sortIndex, 1) tmp(sortIndex, 2)];
wls_Warrant = (300:700)';
spd_Warrant = QuantaToEnergy(wls_Warrant, interp1(tmp(:, 1), tmp(:, 2), wls_Warrant));
plot(wls_Warrant, -0.2+0.08*(spd_Warrant/max(spd_Warrant)), '-k');
plot([280 280], [-0.2 -0.08], '-k');

% Fit all model variants
for s = 1:size(data_merged.spd, 2)
    %% CIE only
    switch data_merged.locIndex(s)
        case 1
            B0 = B_CIE3(:, 1:3);
        case 2
            B0 = B_CIE3(:, 1:3);
    end
    [~, ~, ~, data_merged.gof_CIE(s)] = IlluminationSpectraDataset_Analysis_FitModel(wls, data_merged.spd(:, s), (360:1:830)', B0, 360, 830);
    
    %% CIE with day
    switch data_merged.locIndex(s)
        case 1
            B0 = B_CIE3(:, 1:4);
        case 2
            B0 = B_CIE3(:, 1:4);
    end
    [~, ~, ~, data_merged.gof_CIEMUDay(s)] = IlluminationSpectraDataset_Analysis_FitModel(wls, data_merged.spd(:, s), (360:1:830)', B0, 360, 830);
    
    %% CIE with day + civic
    switch data_merged.locIndex(s)
        case 1
            B0 = B_CIE3(:, 1:5);
        case 2
            B0 = B_CIE3(:, 1:5);
    end
    [~, ~, ~, data_merged.gof_CIEMUDayMUCivic(s)] = IlluminationSpectraDataset_Analysis_FitModel(wls, data_merged.spd(:, s), (360:1:830)', B0, 360, 830);
    
    %% CIE with day + civic + astronomical
    switch data_merged.locIndex(s)
        case 1
            B0 = B_CIE3(:, 1:6);
        case 2
            B0 = B_CIE3(:, [1:5 7]);
    end
    [~, ~, w, data_merged.gof_CIEMUDayMUCivicMUAstronomical(s), ~, data_merged.scalars_CIEMUDayMUCivicMUAstronomical(s)] = IlluminationSpectraDataset_Analysis_FitModel(wls, data_merged.spd(:, s), (360:1:830)', B0, 360, 830);
    data_merged.w_CIEMUDayMUCivicMUAstronomical(:, s) = w;
end

% Save the basis functions
B_CIE3R = B_CIE3(:, 1:6);
B_CIE3C = B_CIE3(:, [1:5 7]);
wls_CIE3R = wls(validWlRngIdxModel);
wls_CIE3C = wls(validWlRngIdxModel);
dataResultsPath = fullfile(fileparts(mfilename('fullpath')), 'data', 'B_CIE3x');
save(fullfile(dataResultsPath, 'B_CIE3x.mat'), 'B_CIE3R', 'B_CIE3C', 'wls_CIE3R', 'wls_CIE3C');
dlmwrite(fullfile(dataResultsPath, 'B_CIE3R.csv'), [wls_CIE3R B_CIE3R], 'delimiter', ',', 'precision',10);
dlmwrite(fullfile(dataResultsPath, 'B_CIE3C.csv'), [wls_CIE3C B_CIE3C], 'delimiter', ',', 'precision',10);

%% Bin the GOF values
theIdx1 = find(data_merged.moonless & (data_merged.locIndex == 1)) ;
[~, ~, ~, ~, gofMean_CIE1, ~, gofSD_CIE1] = bindataflex(data_merged.solarElevationAngleDeg(theIdx1), data_merged.gof_CIE(theIdx1), solarElevationLevels);
[~, ~, ~, ~, gofMean_CIEMUDay1, ~, gofSD_CIEMUDay1] = bindataflex(data_merged.solarElevationAngleDeg(theIdx1), data_merged.gof_CIEMUDay(theIdx1), solarElevationLevels);
[~, ~, ~, ~, gofMean_CIEMUDayMUCivic1, ~, gofSD_CIEMUDayCivic1] = bindataflex(data_merged.solarElevationAngleDeg(theIdx1), data_merged.gof_CIEMUDayMUCivic(theIdx1), solarElevationLevels);
[~, ~, ~, ~, gofMean_CIEMUDayMUCivicMUAstronomical1, ~, gofSD_CIEMUDayCivicMUAstronomical1] = bindataflex(data_merged.solarElevationAngleDeg(theIdx1), data_merged.gof_CIEMUDayMUCivicMUAstronomical(theIdx1), solarElevationLevels);

theIdx2 = find(data_merged.moonless & (data_merged.locIndex == 2)) ;
[~, ~, ~, ~, gofMean_CIE2, ~, gofSD_CIE2] = bindataflex(data_merged.solarElevationAngleDeg(theIdx2), data_merged.gof_CIE(theIdx2), solarElevationLevels);
[~, ~, ~, ~, gofMean_CIEMUDay2, ~, gofSD_CIEMUDay2] = bindataflex(data_merged.solarElevationAngleDeg(theIdx2), data_merged.gof_CIEMUDay(theIdx2), solarElevationLevels);
[~, ~, ~, ~, gofMean_CIEMUDayMUCivic2, ~, gofSD_CIEMUDayCivic2] = bindataflex(data_merged.solarElevationAngleDeg(theIdx2), data_merged.gof_CIEMUDayMUCivic(theIdx2), solarElevationLevels);
[~, ~, ~, ~, gofMean_CIEMUDayMUCivicMUAstronomical2, ~, gofSD_CIEMUDayCivicMUAstronomical2] = bindataflex(data_merged.solarElevationAngleDeg(theIdx2), data_merged.gof_CIEMUDayMUCivicMUAstronomical(theIdx2), solarElevationLevels);

% Run the analysis for Figure S2 to get the GOF when fitting just the night
% spectrum.
IlluminationSpectraDataset_Analysis_FigureS2;

% Plot the R^2 values
subplot(2, 4, 3);
lineProps.col{1} = locRGB{1};
lineProps.width = 2;
mseb(solarElevationLevels, gofMean_CIEMUDayMUCivicMUAstronomical1, gofSD_CIEMUDayCivicMUAstronomical1, lineProps); hold on;
plot(solarElevationLevels, gofMean_CIE1, 'LineWidth', 0.75, 'Color', [142 142 142]/255);
plot(solarElevationLevels, gofMean_CIEMUDay1, 'LineWidth', 0.75, 'Color', [188 151 150]/255);
plot(solarElevationLevels, gofMean_CIEMUDayMUCivic1, 'LineWidth', 0.75, 'Color', [226 132 125]/255);
plot([min(solarElevationLevels) max(solarElevationLevels)], [data(1).gof_nightMean data(1).gof_nightMean]); 
ylims = [0 1];
plot([0 0], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-6 -6], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-12 -12], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-18 -18], yLims, '--', 'Color', [0.7 0.7 0.7]); % Orientation lines
pbaspect([1 1 1]);
xlim([-40 60]); set(gca, 'XTick', [-40:20:60]);
set(gca, 'TickDir', 'out'); box off;
ylim([0 1]);
xlabel('Solar elevation [deg]');
ylabel('R^2');

% Also fit the DiCarlo & Wandell (2000) and the Granada datasets
spd_DiCarloWandell = IlluminationSpectraDataset_Analysis_GetDataset('DiCarlo&Wandell', wls);
spd_Granada = IlluminationSpectraDataset_Analysis_GetDataset('Granada', wls);
[~, ~, ~, gof_DiCarloWandell] = IlluminationSpectraDataset_Analysis_FitModel(wls, spd_DiCarloWandell, (360:1:830)', B_CIE3R, 380, 780);
[~, ~, ~, gof_Granada] = IlluminationSpectraDataset_Analysis_FitModel(wls, spd_Granada, (360:1:830)', B_CIE3R, 380, 780);
errorbar(59, mean(gof_DiCarloWandell), std(gof_DiCarloWandell), '-k');
plot(59, mean(gof_DiCarloWandell), '-sk', 'MarkerFaceColor', 'k');
errorbar(52, mean(gof_Granada), std(gof_Granada), '-k');
plot(52, mean(gof_Granada), '-sk', 'MarkerFaceColor', 'k');

% Plot the R^2 values
subplot(2, 4, 7);
lineProps.col{1} = locRGB{2};
lineProps.width = 2;
mseb(solarElevationLevels, gofMean_CIEMUDayMUCivicMUAstronomical2, gofSD_CIEMUDayCivicMUAstronomical2, lineProps); hold on;
plot(solarElevationLevels, gofMean_CIE2, 'LineWidth', 0.75, 'Color', [142 142 142]/255);
plot(solarElevationLevels, gofMean_CIEMUDay2, 'LineWidth', 0.75, 'Color', [125 159 183]/255);
plot(solarElevationLevels, gofMean_CIEMUDayMUCivic2, 'LineWidth', 0.75, 'Color', [86 147 193]/255);
plot([min(solarElevationLevels) max(solarElevationLevels)], [data(2).gof_nightMean data(2).gof_nightMean]);
plot([0 0], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-6 -6], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-12 -12], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-18 -18], yLims, '--', 'Color', [0.7 0.7 0.7]); % Orientation lines
pbaspect([1 1 1]);
xlim([-40 60]); set(gca, 'XTick', [-40:20:60]);
set(gca, 'TickDir', 'out'); box off;
ylim([0 1]);
xlabel('Solar elevation [deg]');
ylabel('R^2')

% Also fit the DiCarlo & Wandell (2000) and the Granada datasets
spd_DiCarloWandell = IlluminationSpectraDataset_Analysis_GetDataset('DiCarlo&Wandell', wls);
spd_Granada = IlluminationSpectraDataset_Analysis_GetDataset('Granada', wls);
[~, ~, ~, gof_DiCarloWandell] = IlluminationSpectraDataset_Analysis_FitModel(wls, spd_DiCarloWandell, (360:1:830)', B_CIE3C, 380, 780);
[~, ~, ~, gof_Granada] = IlluminationSpectraDataset_Analysis_FitModel(wls, spd_Granada, (360:1:830)', B_CIE3C, 380, 780);
errorbar(59, mean(gof_DiCarloWandell), std(gof_DiCarloWandell), '-k');
plot(59, mean(gof_DiCarloWandell), '-sk', 'MarkerFaceColor', 'k');
errorbar(52, mean(gof_Granada), std(gof_Granada), '-k');
plot(52, mean(gof_Granada), '-sk', 'MarkerFaceColor', 'k');

%% Bin the w values
theIdx1 = find(data_merged.moonless & (data_merged.locIndex == 1)) ;
for k = 1:6
    [~, ~, ~, ~, w_mean1(:, k), ~, w_SD1(:, k)] = bindataflex(data_merged.solarElevationAngleDeg(theIdx1), ...
        data_merged.w_CIEMUDayMUCivicMUAstronomical(k, theIdx1), solarElevationLevels);
end

theIdx2 = find(data_merged.moonless & (data_merged.locIndex == 2)) ;
for k = 1:6
    [~, ~, ~, ~, w_mean2(:, k), ~, w_SD2(:, k)] = bindataflex(data_merged.solarElevationAngleDeg(theIdx2), ...
        data_merged.w_CIEMUDayMUCivicMUAstronomical(k, theIdx2), solarElevationLevels);
end

% Plot the weights
subplot(2, 4, 4);
plot(solarElevationLevels, w_mean1(:, 4), 'LineWidth', 1, 'Color', [14 127 46]/255); hold on;
plot(solarElevationLevels, w_mean1(:, 5), 'LineWidth', 1, 'Color', [115 81 191]/255); hold on;
plot(solarElevationLevels, w_mean1(:, 6), 'LineWidth', 1, 'Color', locRGB{1}); hold on;
yLims = [-0.5 1.5];
plot([0 0], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-6 -6], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-12 -12], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-18 -18], yLims, '--', 'Color', [0.7 0.7 0.7]); % Orientation lines
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off;
xlim([-40 60]); set(gca, 'XTick', [-40:20:60]);
ylim([-0.5 1.5]);
xlabel('Solar elevation [deg]'); ylabel('Component loading');

% Plot the weights
subplot(2, 4, 8);
plot(solarElevationLevels, w_mean2(:, 4), 'LineWidth', 1, 'Color', [14 127 46]/255); hold on;
plot(solarElevationLevels, w_mean2(:, 5), 'LineWidth', 1, 'Color', [115 81 191]/255); hold on;
plot(solarElevationLevels, w_mean2(:, 6), 'LineWidth', 1, 'Color', locRGB{2}); hold on;
yLims = [-0.5 1.5];
plot([0 0], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-6 -6], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-12 -12], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-18 -18], yLims, '--', 'Color', [0.7 0.7 0.7]); % Orientation lines
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off;
xlim([-40 60]); set(gca, 'XTick', [-40:20:60]);
ylim([-0.5 1.5]);
xlabel('Solar elevation [deg]'); ylabel('Component loading');

% Save the graph
set(gcFig6, 'PaperPosition', [0 0 10 5]);
set(gcFig6, 'PaperSize', [10 5]);
saveas(gcFig6, fullfile(resultsPath, ['IlluminationSpectraDataset_Analysis_Figure6.pdf']), 'pdf');
close(gcFig6);

fprintf('Done.');

%% Bin the w values with finer spacing
theIdx1 = find(data_merged.moonless & (data_merged.locIndex == 1)) ;
for k = 1:6
    [~, ~, ~, ~, w_mean1_fine(:, k), ~, w_SD1_fine(:, k)] = bindataflex(data_merged.solarElevationAngleDeg(theIdx1), ...
        data_merged.w_CIEMUDayMUCivicMUAstronomical(k, theIdx1), solarElevationLevels);
end
[~, ~, ~, ~, scaleFactor_mean1_fine, ~, scaleFactor_SD1_fine] = bindataflex(data_merged.solarElevationAngleDeg(theIdx1), ...
        data_merged.scalars_CIEMUDayMUCivicMUAstronomical(theIdx1), solarElevationLevels);

theIdx2 = find(data_merged.moonless & (data_merged.locIndex == 2)) ;
for k = 1:6
    [~, ~, ~, ~, w_mean2_fine(:, k), ~, w_SD2_fine(:, k)] = bindataflex(data_merged.solarElevationAngleDeg(theIdx2), ...
        data_merged.w_CIEMUDayMUCivicMUAstronomical(k, theIdx2), solarElevationLevels);
end
[~, ~, ~, ~, scaleFactor_mean2_fine, ~, scaleFactor_SD2_fine] = bindataflex(data_merged.solarElevationAngleDeg(theIdx2), ...
        data_merged.scalars_CIEMUDayMUCivicMUAstronomical(theIdx2), solarElevationLevels);
    
dlmwrite(fullfile(dataResultsPath, 'w_CIE3R.csv'), [solarElevationLevels' scaleFactor_mean1_fine' scaleFactor_SD1_fine' w_mean1_fine(:, 1) w_SD1_fine(:, 1) w_mean1_fine(:, 2) w_SD1_fine(:, 2) w_mean1_fine(:, 3) w_SD1_fine(:, 3) w_mean1_fine(:, 4) w_SD1_fine(:, 4) w_mean1_fine(:, 5) w_SD1_fine(:, 5) w_mean1_fine(:, 6) w_SD1_fine(:, 6)], 'delimiter', ',', 'precision',10);
dlmwrite(fullfile(dataResultsPath, 'w_CIE3C.csv'), [solarElevationLevels' scaleFactor_mean2_fine' scaleFactor_SD2_fine' w_mean2_fine(:, 1) w_SD2_fine(:, 1) w_mean2_fine(:, 2) w_SD2_fine(:, 2) w_mean2_fine(:, 3) w_SD2_fine(:, 3) w_mean2_fine(:, 4) w_SD2_fine(:, 4) w_mean2_fine(:, 5) w_SD2_fine(:, 5) w_mean2_fine(:, 6) w_SD2_fine(:, 6)], 'delimiter', ',', 'precision',10);