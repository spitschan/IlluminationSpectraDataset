%% IlluminationSpectraDataset_Analysis_Figure1
%
% Generates Figure 1.
%
% 12/10/2015  spitschan   Wrote it.

fprintf('\n>> Making Figure 1...');

for locIndex = [1 2] % Iterate over locations
    figure;
    theRGB = parula(length(solarElevationLevels));
    normIdx = find(wls == 555);
    
    % First, we pull out the night spectra
    if strcmp(theLocations(locIndex).name, 'CSSP')
        data(locIndex).spdAirglow = data(locIndex).spd(:, data(locIndex).theIdxNight);
        data(locIndex).spdAirglowCSSP = data(locIndex).spdAirglow;
        datesAirglowCSSP = dateID(data(locIndex).theIdxNight);
        data(locIndex).spdAirglow_mean = nanmean(data(locIndex).spdAirglowCSSP, 2);
        data(locIndex).spdAirglow_meanNorm = data(locIndex).spdAirglow_mean./data(locIndex).spdAirglow_mean(normIdx);
    elseif strcmp(theLocations(locIndex).name, 'DRL')
        data(locIndex).spdAirglow = data(locIndex).spd(:, data(locIndex).theIdxNight);
        data(locIndex).spdAirglowDRL = data(locIndex).spdAirglow;
        datesAirglowDRL = dateID(data(locIndex).theIdxNight);
        lightPollution = mean(data(locIndex).spdAirglow, 2);
        data(locIndex).spdAirglow_mean = nanmean(data(locIndex).spdAirglowDRL, 2);
        data(locIndex).spdAirglow_meanNorm = data(locIndex).spdAirglow_mean./data(locIndex).spdAirglow_mean(normIdx);
    end
    
    % Night
    subplot(2, 5, 1);
    plot(wls, data(locIndex).spdAirglow_meanNorm, 'Color', [0 0 0], 'LineWidth', 1); hold on;
    plot([555 555], [1 1], 'ok', 'MarkerFaceColor', 'k'); hold on;
    xlim([230 890]); pbaspect([1 1 1]);  ylim([0 4.5]); % used to be 4.2
    box off; set(gca, 'TickDir', 'out');
    title({'Night' 'theta < -18°'});
    switch locIndex
        case 1
            plot([558 558], [2 2], '.r', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
        case 2
            plot([819 819], [3.5 3.5], '.r', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
            plot([570 615], [3.5 3.5], '-r');
    end
    
    % Night (log)
    subplot(2, 5, 6);
    
    plot(wls, log10(data(locIndex).spdAirglow_mean), 'Color', [0 0 0], 'LineWidth', 1); hold on;
    xlim([230 890]); pbaspect([1 1 1]); plot([260 260], [-10 0], '-k');
    plot([235 260], [-10 -10], '-k'); plot([235 260], [-5 -5], '-k'); plot([235 260], [0 0], '-k');
    axis off; set(gca, 'TickDir', 'out'); ylim([-10 0]);
    
    % Then, we pull out some moonless spectra
    % First, we pull out the night spectra
    if strcmp(theLocations(locIndex).name, 'CSSP')
        idx = data(1).theIdx;
    elseif strcmp(theLocations(locIndex).name, 'DRL')
        idx = data(2).theIdx;
    end
    
    data(locIndex).spdMoonless = [data(locIndex).spd(:, idx)];
    solarElevationTwilight = data(locIndex).solarElevationDeg(:, idx);
    nBands = size(data(locIndex).spdMoonless, 1);
    
    % Bin the spectra
    for i = 1:nBands
        [~, ~, ~, ~, data(locIndex).spdMoonlessAvg(i, :)] = bindataflex(solarElevationTwilight, data(locIndex).spdMoonless(i, :), solarElevationLevels);
    end
    
    % Normalize the spectra
    normIdx = find(wls == 555);
    for i = 1:size(data(locIndex).spdMoonlessAvg, 2)
        data(locIndex).spdMoonlessAvgNorm(:, i) = data(locIndex).spdMoonlessAvg(:, i)/data(locIndex).spdMoonlessAvg(normIdx, i);
    end
    
    % Astronomical twilight
    theIdx = find(solarElevationLevels >= -18 & solarElevationLevels < -12);
    theAstronomicalTwilightRGB = hex2rgb(['00004D' ; '000099' ; '0000FF' ;  '1A1AFF' ; '3333FF' ; '6666FF'])/255;
    c = 1;
    for i = 1:length(theIdx)
        % Astronomical twilight (normalized)
        subplot(2, 5, 2);
        plot(wls, data(locIndex).spdMoonlessAvgNorm(:, theIdx(i)), 'Color', theAstronomicalTwilightRGB(c, :)); hold on;
        
        % Astronomical twilight (log)
        subplot(2, 5, 7);
        plot(wls, log10(data(locIndex).spdMoonlessAvg(:, theIdx(i))), 'Color', theAstronomicalTwilightRGB(c, :)); hold on;
        c = c+1;
    end
    subplot(2, 5, 2);
    plot([555 555], [1 1], 'ok', 'MarkerFaceColor', 'k');
    pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');  xlim([230 890]); ylim([0 4.5]);
    title({'Astronomical twilight' '-18° < theta < -12°'});
    subplot(2, 5, 7);
    xlim([230 890]); pbaspect([1 1 1]); plot([260 260], [-10 0], '-k');
    plot([235 260], [-10 -10], '-k'); plot([235 260], [-5 -5], '-k'); plot([235 260], [0 0], '-k');
    axis off; set(gca, 'TickDir', 'out'); ylim([-10 0]);
    
    % Nautical twilight
    theNauticalTwilightRGB = hex2rgb(['004C80' ; '005C99' ; '006BB3' ;  '007ACC' ; '008AE6' ; '0099FF'])/255;
    c = 1;
    theIdx = find(solarElevationLevels >= -12 & solarElevationLevels < -6);
    for i = 1:length(theIdx)
        % Nautical twilight (normalized)
        subplot(2, 5, 3);
        plot(wls, data(locIndex).spdMoonlessAvgNorm(:, theIdx(i)), 'Color', theNauticalTwilightRGB(c, :)); hold on;
        
        % Nautical twilight (log)
        subplot(2, 5, 8);
        plot(wls, log10(data(locIndex).spdMoonlessAvg(:, theIdx(i))), 'Color', theNauticalTwilightRGB(c, :)); hold on;
        c = c+1;
    end
    subplot(2, 5, 3);
    plot([555 555], [1 1], 'ok', 'MarkerFaceColor', 'k');
    pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');  xlim([230 890]); ylim([0 4.5]);
    title({'Nautical twilight' '-12° < theta < -6°'});
    subplot(2, 5, 8);
    xlim([230 890]); pbaspect([1 1 1]); plot([260 260], [-10 0], '-k');
    plot([235 260], [-10 -10], '-k'); plot([235 260], [-5 -5], '-k'); plot([235 260], [0 0], '-k');
    axis off; set(gca, 'TickDir', 'out'); ylim([-10 0]);
    
    
    % Civic twilight
    theCivicTwilightRGB = hex2rgb(['006060' ; '006E6E' ; '007D7D' ;  '009999' ; '00BEBE' ; '00E3E3'])/255;
    c = 1;
    theIdx = find(solarElevationLevels >= -6 & solarElevationLevels < 0);
    for i = 1:length(theIdx)
        % Civic twilight (normalized)
        subplot(2, 5, 4);
        plot(wls, data(locIndex).spdMoonlessAvgNorm(:, theIdx(i)), 'Color', theCivicTwilightRGB(c, :)); hold on;
        
        % Civic twilight (log)
        subplot(2, 5, 9);
        plot(wls, log10(data(locIndex).spdMoonlessAvg(:, theIdx(i))), 'Color', theCivicTwilightRGB(c, :)); hold on;
        c = c+1;
    end
    subplot(2, 5, 4);
    plot([555 555], [1 1], 'ok', 'MarkerFaceColor', 'k');
    pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');  xlim([230 890]); ylim([0 4.5]);
    title({'Civic twilight' '-6° < theta < -0°'});
    subplot(2, 5, 9);
    xlim([230 890]); pbaspect([1 1 1]); plot([260 260], [-10 0], '-k');
    plot([235 260], [-10 -10], '-k'); plot([235 260], [-5 -5], '-k'); plot([235 260], [0 0], '-k');
    axis off; set(gca, 'TickDir', 'out'); ylim([-10 0]);
    
    % Day
    theIdx = find(solarElevationLevels >= 0);
    theRGB = summer(length(theIdx));
    c = 1;
    for i = 1:length(theIdx)
        % Day (normalized(
        subplot(2, 5, 5);
        plot(wls, data(locIndex).spdMoonlessAvgNorm(:, theIdx(i)), 'Color', theRGB(c, :)); hold on;
        
        % Day (log)
        subplot(2, 5, 10);
        plot(wls, log10(data(locIndex).spdMoonlessAvg(:, theIdx(i))), 'Color', theRGB(c, :)); hold on;
        c = c+1;
    end
    subplot(2, 5, 5);
    plot([555 555], [1 1], 'ok', 'MarkerFaceColor', 'k');
    pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');  xlim([230 890]); ylim([0 4.5]);
    title({'Day' 'theta > 0°'});
    subplot(2, 5, 10);
    xlim([230 890]); pbaspect([1 1 1]); plot([260 260], [-10 0], '-k');
    plot([235 260], [-10 -10], '-k'); plot([235 260], [-5 -5], '-k'); plot([235 260], [0 0], '-k');
    axis off; set(gca, 'TickDir', 'out'); ylim([-10 0]);
    
    % Save figure
    set(gcf, 'PaperPosition', [0 0 12 6]);
    set(gcf, 'PaperSize', [12 6]);
    switch locIndex
        case 1
            saveas(gcf, fullfile(resultsPath, 'IlluminationSpectraDataset_Analysis_Figure1A.pdf'), 'pdf');
        case 2
            saveas(gcf, fullfile(resultsPath, 'IlluminationSpectraDataset_Analysis_Figure1B.pdf'), 'pdf');
    end
    close(gcf);
end

% Save out the color maps
if ~isdir(fullfile(resultsPath, 'ColorMaps'));
   mkdir(fullfile(resultsPath, 'ColorMaps'));
end
colormap(theAstronomicalTwilightRGB); colorbar; axis off;
set(gcf, 'PaperPosition', [0 0 2 6]);
set(gcf, 'PaperSize', [2 6]);
saveas(gcf, fullfile(resultsPath, 'ColorMaps', 'IlluminationSpectraDataset_Analysis_Figure1_ColorMapAstronomicalTwilight.pdf'), 'pdf');
close(gcf);

colormap(theNauticalTwilightRGB); colorbar; axis off;
set(gcf, 'PaperPosition', [0 0 2 6]);
set(gcf, 'PaperSize', [2 6]);
saveas(gcf, fullfile(resultsPath, 'ColorMaps', 'IlluminationSpectraDataset_Analysis_Figure1_ColorMapNauticalTwilight.pdf'), 'pdf');
close(gcf);

colormap(theCivicTwilightRGB); colorbar; axis off;
set(gcf, 'PaperPosition', [0 0 2 6]);
set(gcf, 'PaperSize', [2 6]);
saveas(gcf, fullfile(resultsPath, 'ColorMaps', 'IlluminationSpectraDataset_Analysis_Figure1_ColorMapCivicTwilight.pdf'), 'pdf');
close(gcf);

colormap(theRGB); colorbar; axis off;
set(gcf, 'PaperPosition', [0 0 2 6]);
set(gcf, 'PaperSize', [2 6]);
saveas(gcf, fullfile(resultsPath, 'ColorMaps', 'IlluminationSpectraDataset_Analysis_Figure1_ColorMapDay.pdf'), 'pdf');
close(gcf);

fprintf('Done.');