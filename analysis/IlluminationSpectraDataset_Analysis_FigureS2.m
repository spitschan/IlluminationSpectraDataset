%% IlluminationSpectraDataset_Analysis_FigureS2
%
% Generates Figure S2.
%
% 12/10/2015  spitschan   Wrote it.

fprintf('\n>> Making Figure S2...');
gcFigS2 = figure;

wlRangeStartIdxModel = find(wls == 360);
wlRangeEndIdxModel = find(wls == 830);
validWlRngIdxModel = wlRangeStartIdxModel:wlRangeEndIdxModel;

%% Fit the night data with their mean
for locIndex = [1 2]
    % Get the mean night spectrum
    B_night = mean(data(locIndex).spd(:, data(locIndex).theIdxNight), 2);
    B_night = B_night/norm(B_night(validWlRngIdxModel, :));
    
    % Iterate over the relevant spectra
    [~, ~, ~, gof] = IlluminationSpectraDataset_Analysis_FitModel(wls, data(locIndex).spd, wls, B_night, 360, 830);
    
    % Bin the GOF
    [~, ~, gofMin, gofMax, gofMean, gofMedian, gofSD] = bindataflex(data(locIndex).solarElevationDeg(data(locIndex).theIdx), gof(data(locIndex).theIdx), solarElevationLevels);
    
    % Plot the R^2 values
    subplot(2, 2, 2*(locIndex-1)+1);
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
    
    data(locIndex).gof_nightMean = mean(gof(data(locIndex).theIdxNight));
    
    % Plot the R^2 values
    subplot(2, 2, 2*(locIndex-1)+2);
    plot(wls(validWlRngIdxModel), B_night(validWlRngIdxModel, :), 'Color', locRGB{locIndex});
    pbaspect([1 1 1]);
    xlabel('Wavelength [nm[');
    ylabel('Relative value');
    xlim([230 890]); pbaspect([1 1 1]); set(gca, 'TickDir', 'out'); box off;
    ylim([0 0.25]);
end

set(gcFigS2, 'PaperPosition', [0 0 5 5]);
set(gcFigS2, 'PaperSize', [5 5]);
saveas(gcFigS2, fullfile(resultsPath, ['IlluminationSpectraDataset_Analysis_FigureS2.pdf']), 'pdf');
close(gcFigS2);

fprintf('Done.');