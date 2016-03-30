%% IlluminationSpectraDataset_Analysis_Figure2
%
% Generates Figure 2.
%
% 12/10/2015  spitschan   Wrote it.

fprintf('\n>> Making Figure 2...');
gcFig2 = figure;

for locIndex = [1 2]
    %% UV-A
    totalIrradianceRangeUVA = [280 315];
    startIdx = find(wls == totalIrradianceRangeUVA(1));
    endIdx = find(wls == totalIrradianceRangeUVA(2));
    totalIrradianceUVA = nansum(data(locIndex).spd(startIdx:endIdx, data(locIndex).theIdx));
    solarElevationDeg = data(locIndex).solarElevationDeg(data(locIndex).theIdx);
    [~, ~, ~, ~, logTotalIrradianceUVAMean, ~, logTotalIrradianceUVASD] = bindataflex(solarElevationDeg(~(totalIrradianceUVA == 0)), log10(totalIrradianceUVA(~(totalIrradianceUVA == 0))), solarElevationLevels);
    data(locIndex).totalIrradianceUVA = totalIrradianceUVA;
    
    %% UV-B
    totalIrradianceRangeUVB = [315 400];
    startIdx = find(wls == totalIrradianceRangeUVB(1));
    endIdx = find(wls == totalIrradianceRangeUVB(2));
    totalIrradianceUVB = nansum(data(locIndex).spd(startIdx:endIdx, data(locIndex).theIdx));
    solarElevationDeg = data(locIndex).solarElevationDeg(data(locIndex).theIdx);
    [~, ~, ~, ~, logTotalIrradianceUVBMean, ~, logTotalIrradianceUVBSD] = bindataflex(solarElevationDeg(~(totalIrradianceUVB == 0)), log10(totalIrradianceUVB(~(totalIrradianceUVB == 0))), solarElevationLevels);
    data(locIndex).totalIrradianceUVB = totalIrradianceUVB;
    
    %% VIS
    totalIrradianceRangeVIS = [400 800];
    startIdx = find(wls == totalIrradianceRangeVIS(1));
    endIdx = find(wls == totalIrradianceRangeVIS(2));
    totalIrradianceVIS = nansum(data(locIndex).spd(startIdx:endIdx, data(locIndex).theIdx));
    solarElevationDeg = data(locIndex).solarElevationDeg(data(locIndex).theIdx);
    [~, ~, ~, ~, logTotalIrradianceVISMean, ~, logTotalIrradianceVISSD] = bindataflex(solarElevationDeg(~(totalIrradianceVIS == 0)), log10(totalIrradianceVIS(~(totalIrradianceVIS == 0))), solarElevationLevels);
    data(locIndex).totalIrradianceVIS = totalIrradianceVIS;
    data(locIndex).logTotalIrradianceVISMean = logTotalIrradianceVISMean;
    
    maxVISIrr(locIndex) = max(logTotalIrradianceVISMean);
    minVISIrr(locIndex) = min(logTotalIrradianceVISMean);
    
    %% VIS
    yLims = [-8 4];
    dotMarkerSize1 = 2;
    subplot(2, 3, (locIndex-1)*3+1);
    lineProps.col{1} = [0.5647 0.7333 0.5647];
    lineProps.width = 1;
    mseb(solarElevationLevels, logTotalIrradianceVISMean, logTotalIrradianceVISSD, lineProps); hold on;
    xlim([-40 60]); ylim(yLims);
    set(gca, 'YTick', [-8:4]);
    set(gca, 'XTick', [-40:20:60]);
    pbaspect([1 1 1]);  box off;
    set(gca, 'TickDir', 'out');
    xlabel('Solar elevation [deg]');
    ylabel('Total irradiance [W m2]');
    plot([0 0], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-6 -6], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-12 -12], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-18 -18], yLims, '--', 'Color', [0.7 0.7 0.7]); % Orientation lines
    
    %% Add reference lines
    if locIndex == 2
        plot([-35 70], [minVISIrr(1) minVISIrr(1)], '--k');
        plot([-35 70], [maxVISIrr(1) maxVISIrr(1)], '--k');
    end
    
    subplot(2, 3, (locIndex-1)*3+2);
    % Evening
    if strcmp(theLocations(locIndex).name, 'CSSP')
        dates = [1 3 4 6 7 8];
    elseif strcmp(theLocations(locIndex).name, 'DRL')
        dates = [1 3 4 5 6 8];
    end

    theIdx = intersect(find(ismember(data(locIndex).date, dates)), data(locIndex).theIdx);
    totalIrradianceVISEvening = nansum(data(locIndex).spd(startIdx:endIdx, theIdx));
    solarElevationDeg = data(locIndex).solarElevationDeg(theIdx);
    [~, ~, ~, ~, data(locIndex).logTotalIrradianceBinnedEveningMean, ~, data(locIndex).logTotalIrradianceBinnedEveningSD] = bindataflex(solarElevationDeg, log10(totalIrradianceVISEvening), solarElevationLevels);

    % Morning
    if strcmp(theLocations(locIndex).name, 'CSSP')
        dates = [2 5 9];
    elseif strcmp(theLocations(locIndex).name, 'DRL')
        dates = [2 7 9];
    end
    
    theIdx = intersect(find(ismember(data(locIndex).date, dates)), data(locIndex).theIdx);
    totalIrradianceVISMorning = nansum(data(locIndex).spd(startIdx:endIdx, theIdx));
    solarElevationDeg = data(locIndex).solarElevationDeg(theIdx);
    [~, ~, ~, ~, data(locIndex).logTotalIrradianceBinnedMorningMean, ~, data(locIndex).logTotalIrradianceBinnedMorningSD] = bindataflex(solarElevationDeg, log10(totalIrradianceVISMorning), solarElevationLevels);
    
    lineProps.col{1} = [246 74 71]/255; lineProps.width = 1;
    mseb(solarElevationLevels, data(locIndex).logTotalIrradianceBinnedMorningMean, data(locIndex).logTotalIrradianceBinnedMorningSD, lineProps); hold on;
    lineProps.col{1} = [127 139 109]/255; lineProps.width = 1;
    mseb(solarElevationLevels, data(locIndex).logTotalIrradianceBinnedEveningMean, data(locIndex).logTotalIrradianceBinnedEveningSD, lineProps); hold on;
    plot([0 0], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-6 -6], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-12 -12], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-18 -18], yLims, '--', 'Color', [0.7 0.7 0.7]); % Orientation lines
    
    ylim([-8 4]);
    set(gca, 'YTick', [-8:4]);
    xlim([-40 60]);
    set(gca, 'XTick', [-40:20:60]);
    pbaspect([1 1 1]);  box off;
    set(gca, 'TickDir', 'out');
    xlabel('Solar elevation [deg]');
    ylabel('Total irradiance [W m2]');
    
    %% UV-A
    dotMarkerSize1 = 2;
    subplot(2, 3, (locIndex-1)*3+3);
    lineProps.col{1} = [147 112 219]/255;
    lineProps.width = 1;
    h1 = mseb(solarElevationLevels, logTotalIrradianceUVAMean, logTotalIrradianceUVASD, lineProps); hold on;
    
    %% UV-B
    dotMarkerSize1 = 2;
    lineProps.width = 1;
    lineProps.col{1} = [65 105 225]/255;
    h2 = mseb(solarElevationLevels, logTotalIrradianceUVBMean, logTotalIrradianceUVBSD, lineProps);
    plot([0 0], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-6 -6], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-12 -12], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-18 -18], yLims, '--', 'Color', [0.7 0.7 0.7]); % Orientation lines
    xlim([-40 60]);  ylim([-8 4]);
    set(gca, 'YTick', [-8:4]);
    set(gca, 'XTick', [-40:20:60]);
    pbaspect([1 1 1]);  box off;
    set(gca, 'TickDir', 'out');
    xlabel('Solar elevation [deg]');
    ylabel('Total irradiance [W m2]');
end
set(gcFig2, 'PaperPosition', [0 0 10 7]);
set(gcFig2, 'PaperSize', [10 7]);
saveas(gcFig2, fullfile(resultsPath, ['IlluminationSpectraDataset_Analysis_Figure2.pdf']), 'pdf');
close(gcFig2);
fprintf('Done.');