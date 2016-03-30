%% IlluminationSpectraDataset_Analysis_FigureS3
%
% Generates Figure S3.
%
% 12/10/2015  spitschan   Wrote it.

fprintf('\n>> Making Figure S3...');
gcFigS3 = figure;

%% Get the basis functions
B_Granada = IlluminationSpectraDataset_Analysis_GetModel('Granada', wls, true);

subplotIndexGOF = [1 3];
for locIndex = [1 2]
    %% Fit the Granada model
    [spds_norm, spds_fit, w, gof, residuals] = IlluminationSpectraDataset_Analysis_FitModel(wls, data(locIndex).spd, wls, B_Granada, 380, 780);
    
    % Bin the data
    [~, ~, gofMin, gofMax, gofMean, gofMedian, gofSD] = bindataflex(data(locIndex).solarElevationDeg(data(locIndex).theIdx), gof(data(locIndex).theIdx), solarElevationLevels);
    
    % Plot the R^2 values
    subplot(2, 2, subplotIndexGOF(locIndex));
    lineProps.col{1} = [0 0 0];
    lineProps.width = 2;
    mseb(solarElevationLevels, gofMean, gofSD, lineProps); hold on;
    
    %% Fit our model
    load B_CIE3x.mat
    switch locIndex
        case 1
            B_CIE3x = B_CIE3R;
        case 2
            B_CIE3x = B_CIE3C;
    end
    [spds_norm, spds_fit, w, gof, residuals] = IlluminationSpectraDataset_Analysis_FitModel(wls, data(locIndex).spd, (360:1:830)', B_CIE3x, 380, 780);
    
    % Bin the data
    [~, ~, gofMin, gofMax, gofMean, gofMedian, gofSD] = bindataflex(data(locIndex).solarElevationDeg(data(locIndex).theIdx), gof(data(locIndex).theIdx), solarElevationLevels);
    
    % Plot the R^2 values
    subplot(2, 2, subplotIndexGOF(locIndex));
    plot(solarElevationLevels, gofMean, '-', 'Color', locRGB{locIndex}); hold on;
    plot([0 0], [0 1], '-k'); plot([-6 -6], [0 1], '-k'); plot([-12 -12], [0 1], '-k'); plot([-18 -18], [0 1], '-k'); % Orientation lines
    pbaspect([1 1 1]);
    xlim([-40 60]); set(gca, 'XTick', [-40:20:60]);
    set(gca, 'TickDir', 'out'); box off;
    ylim([0 1]);
    xlabel('Solar elevation [deg]');
    ylabel('R^2')
end

% Save the graph
set(gcFigS3, 'PaperPosition', [0 0 5 5]);
set(gcFigS3, 'PaperSize', [5 5]);
saveas(gcFigS3, fullfile(resultsPath, ['IlluminationSpectraDataset_Analysis_FigureS3.pdf']), 'pdf');
close(gcFigS3);

fprintf('Done.');