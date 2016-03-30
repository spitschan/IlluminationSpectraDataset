%% IlluminationSpectraDataset_Analysis_FigureS1
%
% Generates Figure S1.
%
% 12/10/2015  spitschan   Wrote it.


fprintf('\n>> Making Figure S1...');

for locIndex = [1 2]
    gcFigS1 = figure;
    % Find out how many spectral bands there are
    nBands = size(data(locIndex).spd, 1);
    fractionIlluminatedStart = 0:0.1:0.9;
    fractionIlluminatedEnd = (0:0.1:1.0)+0.1;
    
    c = 1;
    theRGB = parula(5);
    normIdx = find(wls == 555);
    for f = 1:length(fractionIlluminatedStart)
        theIdxMoon = find(data(locIndex).solarElevationDeg < -18 & data(locIndex).lunarFractionIlluminated ...
            >= fractionIlluminatedStart(f) & data(locIndex).lunarFractionIlluminated < fractionIlluminatedEnd(f));
        spdNightFracIll(:, f) = nanmean(data(locIndex).spd(:, theIdxMoon), 2);
        
        subplot(1, 2, 1);
        h(c) = plot(wls, spdNightFracIll(:, f)/spdNightFracIll(normIdx, f), '-', 'Color', theRGB(c, :)); hold on;
        subplot(1, 2, 2);
        plot(wls, log10(spdNightFracIll(:, f)), '-', 'Color', theRGB(c, :)); hold on;
        if ~all(isnan(spdNightFracIll(:, f)))
            c = c+1;
        end
    end
    
    % Add the night sky spectrum for city
    % locIndex = 2;
    % subplot(1, 2, 2);
    % plot(wls, log10(mean(data(locIndex).spd(:, data(locIndex).theIdxNight), 2)), '-', 'Color', [0.5 0.5 0.5]);
    %
    % (data(locIndex).lunarFractionIlluminated > 0.9) & (data(locIndex).solarElevationDeg < -18)
    % subplot(1, 2, 1);
    % plot(data(locIndex).spd(:, (data(locIndex).lunarFractionIlluminated > 0.9) & (data(locIndex).solarElevationDeg < -18)))
    % title('City full moon')
    % ylim([0 0.00003]);
    % subplot(1, 2, 2);
    % plot(data(locIndex).spd(:, (data(locIndex).lunarFractionIlluminated < 0.6) & (data(locIndex).solarElevationDeg < -18)))
    % title('City frac < 0.6');
    % ylim([0 0.00003]);
    
    fracStarts = fractionIlluminatedStart(~all(isnan(spdNightFracIll)));
    fracEnds = fractionIlluminatedEnd(~all(isnan(spdNightFracIll)));
    for m = 1:c-1
        theLegendEntries{m} = [num2str(fracStarts(m), '%.2f') '-' num2str(fracEnds(m), '%.2f')];
    end
    
    subplot(1, 2, 1);
    pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');  xlim([230 890]); ylim([0 4.5]);
    title('Relative night spectra');
    xlabel('Wavelength [nm]');
    ylabel('Relative irradiance');
    legend(h, theLegendEntries{1}, theLegendEntries{2}, theLegendEntries{3}, theLegendEntries{4}, 'Location', 'NorthWest') ; legend boxoff;
    
    subplot(1, 2, 2);
    pbaspect([1 1 1]); box off; set(gca, 'TickDir', 'out');  xlim([230 890]); ylim([-10 0]);
    title('log night spectra');
    xlabel('Wavelength [nm]');
    ylabel('log irradiance');
    xlim([230 890]); pbaspect([1 1 1]); plot([260 260], [-10 0], '-k');
    plot([235 260], [-10 -10], '-k'); plot([235 260], [-5 -5], '-k'); plot([235 260], [0 0], '-k');
    axis off; set(gca, 'TickDir', 'out'); ylim([-10 0]);
    
    % Save figure
    set(gcFigS1, 'PaperPosition', [0 0 6 3]);
    set(gcFigS1, 'PaperSize', [6 3]);
    switch locIndex
        case 1
            saveas(gcFigS1, fullfile(resultsPath, 'IlluminationSpectraDataset_Analysis_FigureS1A.pdf'), 'pdf');
        case 2
            saveas(gcFigS1, fullfile(resultsPath, 'IlluminationSpectraDataset_Analysis_FigureS1B.pdf'), 'pdf');
    end
    close(gcFigS1);
end
fprintf('Done.');