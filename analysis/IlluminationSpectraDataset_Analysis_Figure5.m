%% IlluminationSpectraDataset_Analysis_Figure4
%
% Generates Figure 5.
%
% 12/10/2015  spitschan   Wrote it.

fprintf('\n>> Making Figure 5...');
gcFig5 = figure;

%% Load and show the CIE daylight basis functions
B_CIE_360_830 = IlluminationSpectraDataset_Analysis_GetModel('CIE', (360:1:830)', true);
B_CIE_380_780 = IlluminationSpectraDataset_Analysis_GetModel('CIE', (380:1:780)', true);


% Plot the CIE basis functions
subplot(2, 3, 1);
plot((360:1:830)', B_CIE_360_830(:, 1), '-', 'Color', [0 0 0]); hold on;
plot((360:1:830)', B_CIE_360_830(:, 2), '-', 'Color', [0.54 0.54 0.54]); 
plot((360:1:830)', B_CIE_360_830(:, 3), '-', 'Color', [0.35 0.35 0.35]);
xlabel('Wavelength [nm]'); ylabel('Relative value');
pbaspect([1 1 1]);
xlim([230 890]);
ylim([-0.05 0.1]);
pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off;

%% Define some parameters
subplotIndexGOF = [2 5];
subplotIndexWeights = [3 6];
locRGB{1} = [0.85 0.22 0.16];
locRGB{2} = [0.26 0.49 0.76];

for locIndex = [1 2]
    % Fit the CIE model
    [spds_norm, spds_fit, w, gof, residuals] = IlluminationSpectraDataset_Analysis_FitModel(wls, data(locIndex).spd, (360:1:830)', B_CIE_360_830, 360, 830);
    
    data(locIndex).allSpectra_resid = residuals;
    data(locIndex).gofCIE3Components = gof;
    
    % Bin the data
    % Define our relevant spectra
    [~, ~, gofMin, gofMax, gofMean, gofMedian, gofSD] = bindataflex(data(locIndex).solarElevationDeg(data(locIndex).theIdx), gof(data(locIndex).theIdx), solarElevationLevels);
    
    % Plot the R^2 values
    subplot(2, 3, subplotIndexGOF(locIndex));
  
    lineProps.col{1} = locRGB{locIndex};
    lineProps.width = 1;
    mseb(solarElevationLevels, gofMean, gofSD, lineProps); hold on;
    
    yLims = [0 1];
    plot([0 0], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-6 -6], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-12 -12], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-18 -18], yLims, '--', 'Color', [0.7 0.7 0.7]); % Orientation lines
    pbaspect([1 1 1]);
    xlim([-40 60]); set(gca, 'XTick', [-40:20:60]);
    set(gca, 'TickDir', 'out'); box off;
    ylim([0 1]);
    xlabel('Solar elevation [deg]');
    ylabel('R^2')
    
    % Plot the weights
    subplot(2, 3, subplotIndexWeights(locIndex));
    [~, ~, ~, ~, w1Mean, ~, w1SD] = bindataflex(data(locIndex).solarElevationDeg(data(locIndex).theIdx), w(1, data(locIndex).theIdx), solarElevationLevels);
    [~, ~, ~, ~, w2Mean, ~, w2SD] = bindataflex(data(locIndex).solarElevationDeg(data(locIndex).theIdx), w(2, data(locIndex).theIdx), solarElevationLevels);
    [~, ~, ~, ~, w3Mean, ~, w3SD] = bindataflex(data(locIndex).solarElevationDeg(data(locIndex).theIdx), w(3, data(locIndex).theIdx), solarElevationLevels);
    
    lineProps.col{1} = [0 0 0]; 
    lineProps.col{2} = [0.54 0.54 0.54];
    lineProps.col{3} = [0.34 0.34 0.34];
    lineProps.edgestyle = '--';
    lineProps.width = 1;
    mseb(solarElevationLevels, [w1Mean ; w2Mean ; w3Mean], [w1SD ; w2SD ; w3SD], lineProps); hold on;    
    yLims = [-0.5 11.5];
    plot([0 0], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-6 -6], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-12 -12], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-18 -18], yLims, '--', 'Color', [0.7 0.7 0.7]); % Orientation lines
    
    pbaspect([1 1 1]);
    xlim([-40 60]); set(gca, 'XTick', [-40:20:60]);
    set(gca, 'TickDir', 'out'); box off;
    ylim([-0.5 1.5]);
    xlabel('Solar elevation [deg]');
    ylabel('Component loading');
end

% Also fit the DiCarlo & Wandell (2000) and the Granada datasets
spd_DiCarloWandell = IlluminationSpectraDataset_Analysis_GetDataset('DiCarlo&Wandell', wls);
spd_Granada = IlluminationSpectraDataset_Analysis_GetDataset('Granada', wls);

[~, ~, ~, gof_DiCarloWandell] = IlluminationSpectraDataset_Analysis_FitModel(wls, spd_DiCarloWandell, (380:1:780)', B_CIE_380_780, 380, 780);
[~, ~, ~, gof_Granada] = IlluminationSpectraDataset_Analysis_FitModel(wls, spd_Granada, (380:1:780)', B_CIE_380_780, 380, 780);

for f = subplotIndexGOF
    subplot(2, 3, f);
    errorbar(59, mean(gof_DiCarloWandell), std(gof_DiCarloWandell), '-k');
    plot(59, mean(gof_DiCarloWandell), '-sk', 'MarkerFaceColor', 'k');
    errorbar(52, mean(gof_Granada), std(gof_Granada), '-k');
    plot(52, mean(gof_Granada), '-sk', 'MarkerFaceColor', 'k');
end

% Save the figure
set(gcFig5, 'PaperPosition', [0 0 8 5]);
set(gcFig5, 'PaperSize', [8 5]);
saveas(gcFig5, fullfile(resultsPath, ['IlluminationSpectraDataset_Analysis_Figure5.pdf']), 'pdf');
close(gcFig5);

fprintf('Done.');