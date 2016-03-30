%% IlluminationSpectraDataset_Analysis_Figure3
%
% Generates Figure 3.
%
% 12/10/2015  spitschan   Wrote it.

fprintf('\n>> Making Figure 3...');
gcFig3 = figure;

%% Load CIE functions.
S = WlsToS(wls);
load T_xyz1931
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,WlsToS((380:1:780)'));

% Find relevant wl range
startIdx = find(wls == 380);
endIdx = find(wls == 780);

for locIndex = [1 2]
    % Calculate XYZ, luminance and chromaticity x y.
    xyzVals = T_xyz * data(locIndex).spd(startIdx:endIdx, :);
    chromaticityX = xyzVals(1, :)./sum(xyzVals);
    chromaticityY = xyzVals(2, :)./sum(xyzVals);
    data(locIndex).chromaticityX = chromaticityX;
    data(locIndex).chromaticityY = chromaticityY;
    
    [~, ~, ~, ~, data(locIndex).chromXValsAvg, ~, data(locIndex).chromXValsSD] = bindataflex(data(locIndex).solarElevationDeg(data(locIndex).theIdx), data(locIndex).chromaticityX(data(locIndex).theIdx), solarElevationLevels);
    [~, ~, ~, ~, data(locIndex).chromYValsAvg, ~, data(locIndex).chromYValsSD] = bindataflex(data(locIndex).solarElevationDeg(data(locIndex).theIdx), data(locIndex).chromaticityY(data(locIndex).theIdx), solarElevationLevels);
end

% Horse shoe scatter plot
locRGB{1} = [0.85 0.22 0.16];
locRGB{2} = [0.26 0.49 0.76];
subplot(3, 3, 5);
plot(data(1).chromaticityX, data(1).chromaticityY, '.', 'Color', locRGB{1}); hold on;
plot(data(2).chromaticityX, data(2).chromaticityY, '.', 'Color', locRGB{2});

load B_cieday
B_cieday = SplineRaw(S_cieday, B_cieday, (380:1:780)');
B_cieday_raw = B_cieday;
for i = 1:3
    B_cieday(:, i) = B_cieday(:, i)/norm(B_cieday(:, i));
end
[~,xd,yd] = GenerateCIEDay([4000:25000],B_cieday);
[~,xd65,yd65] = GenerateCIEDay(6500,B_cieday);
plot(xd, yd, '-k');

box off; set(gca, 'TickDir', 'out');
set(gca, 'XTick', [0.2 0.3 0.4 0.5]);
set(gca, 'YTick', [0.2 0.3 0.4 0.5]);
xlim([0.15 0.5]); ylim([0.15 0.5]);
pbaspect([1 1 1]);

% Histogram for x chromaticity
subplot(3, 3, 2);
xy = linspace(0, 1, 75);
% Day
histctmp = histc(data(1).chromaticityX(data(1).theIdxDay), xy);
stairs(xy, 0.8*(histctmp/sum(histctmp)), 'Color', locRGB{1}); hold on;
histctmp = histc(data(2).chromaticityX(data(2).theIdxDay), xy);
stairs(xy, 0.8*(histctmp/sum(histctmp)), 'Color', locRGB{2});
plot([mean(data(1).chromaticityX(data(1).theIdxDay))], [0.6], 's', 'Color', locRGB{1});
plot([mean(data(2).chromaticityX(data(2).theIdxDay))], [0.6], 's', 'Color', locRGB{2});

% Twilight
histctmp = histc(data(1).chromaticityX([data(1).theIdxTwilightAstronomical data(1).theIdxTwilightNautical data(1).theIdxTwilightCivic]), xy);
stairs(xy, 1+0.8*(histctmp/sum(histctmp)), 'Color', locRGB{1}); hold on;
histctmp = histc(data(2).chromaticityX([data(2).theIdxTwilightAstronomical data(2).theIdxTwilightNautical data(2).theIdxTwilightCivic]), xy);
stairs(xy, 1+0.8*(histctmp/sum(histctmp)), 'Color', locRGB{2});
plot([mean(data(1).chromaticityX([data(1).theIdxTwilightAstronomical data(1).theIdxTwilightNautical data(1).theIdxTwilightCivic]))], [1.6], 's', 'Color', locRGB{1});
plot([mean(data(2).chromaticityX([data(2).theIdxTwilightAstronomical data(2).theIdxTwilightNautical data(2).theIdxTwilightCivic]))], [1.6], 's', 'Color', locRGB{2});

% Night
histctmp = histc(data(1).chromaticityX(data(1).theIdxNight), xy);
stairs(xy, 2+0.8*(histctmp/sum(histctmp)), 'Color', locRGB{1}); hold on;
histctmp = histc(data(2).chromaticityX(data(2).theIdxNight), xy);
stairs(xy, 2+0.8*(histctmp/sum(histctmp)), 'Color', locRGB{2});
plot([mean(data(1).chromaticityX(data(1).theIdxNight))], [2.6], 's', 'Color', locRGB{1});
plot([mean(data(2).chromaticityX(data(2).theIdxNight))], [2.6], 's', 'Color', locRGB{2});

plot([xd65 xd65], [0 3], '-k');
xlim([0.15 0.5]); ylim([0 3]); pbaspect([1 1 1]); box off; axis off;
set(gca, 'XTick', [0.2 0.3 0.4 0.5]);

% Histogram for x chromaticity
subplot(3, 3, 6);
xy = linspace(0, 1, 75);
% Day
histctmp = histc(data(1).chromaticityY(data(1).theIdxDay), xy);
stairs(xy, 0.8*(histctmp/sum(histctmp)), 'Color', locRGB{1}); hold on;
histctmp = histc(data(2).chromaticityY(data(2).theIdxDay), xy);
stairs(xy, 0.8*(histctmp/sum(histctmp)), 'Color', locRGB{2});
plot([mean(data(1).chromaticityY(data(1).theIdxDay))], [0.6], 's', 'Color', locRGB{1});
plot([mean(data(2).chromaticityY(data(2).theIdxDay))], [0.6], 's', 'Color', locRGB{2});

% Twilight
histctmp = histc(data(1).chromaticityY([data(1).theIdxTwilightAstronomical data(1).theIdxTwilightNautical data(1).theIdxTwilightCivic]), xy);
stairs(xy, 1+0.8*(histctmp/sum(histctmp)), 'Color', locRGB{1}); hold on;
histctmp = histc(data(2).chromaticityY([data(2).theIdxTwilightAstronomical data(2).theIdxTwilightNautical data(2).theIdxTwilightCivic]), xy);
stairs(xy, 1+0.8*(histctmp/sum(histctmp)), 'Color', locRGB{2});
plot([mean(data(1).chromaticityY([data(1).theIdxTwilightAstronomical data(1).theIdxTwilightNautical data(1).theIdxTwilightCivic]))], [1.6], 's', 'Color', locRGB{1});
plot([mean(data(2).chromaticityY([data(2).theIdxTwilightAstronomical data(2).theIdxTwilightNautical data(2).theIdxTwilightCivic]))], [1.6], 's', 'Color', locRGB{2});

% Night
histctmp = histc(data(1).chromaticityY(data(1).theIdxNight), xy);
stairs(xy, 2+0.8*(histctmp/sum(histctmp)), 'Color', locRGB{1}); hold on;
histctmp = histc(data(2).chromaticityY(data(2).theIdxNight), xy);
stairs(xy, 2+0.8*(histctmp/sum(histctmp)), 'Color', locRGB{2});
plot([mean(data(1).chromaticityY(data(1).theIdxNight))], [2.6], 's', 'Color', locRGB{1});
plot([mean(data(2).chromaticityY(data(2).theIdxNight))], [2.6], 's', 'Color', locRGB{2});

plot([yd65 yd65], [0 3], '-k');
xlim([0.15 0.5]); ylim([0 3]); pbaspect([1 1 1]); box off; axis off;
set(gca, 'XTick', [0.2 0.3 0.4 0.5]);
set(gca, 'YDir', 'Reverse');
view(-90,90);

% X chromaticity progressions
yLims = [0.15 0.5];
subplot(3, 3, 8);
lineProps.col{1} = locRGB{1};
lineProps.width = 1;
mseb(solarElevationLevels, data(1).chromXValsAvg, data(1).chromXValsSD, lineProps); hold on;
lineProps.col{1} = locRGB{2};
lineProps.width = 1;
mseb(solarElevationLevels, data(2).chromXValsAvg, data(2).chromXValsSD, lineProps); hold on;
plot([0 0], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-6 -6], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-12 -12], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-18 -18], yLims, '--', 'Color', [0.7 0.7 0.7]); % Orientation lines

ylim(yLims); pbaspect([1 1 1]); box off;
xlim([-40 60]);
set(gca, 'XTick', [-40:20:60]);
set(gca, 'YTick', [0.2 0.3 0.4 0.5]);
set(gca, 'TickDir', 'out');

% Y chromaticity progressions
subplot(3, 3, 4);
lineProps.col{1} = locRGB{1};
lineProps.width = 1;
mseb(solarElevationLevels, data(1).chromYValsAvg, data(1).chromYValsSD, lineProps); hold on;
lineProps.col{1} = locRGB{2};
lineProps.width = 1;
mseb(solarElevationLevels, data(2).chromYValsAvg, data(2).chromYValsSD, lineProps); hold on;
plot([0 0], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-6 -6], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-12 -12], yLims, '--', 'Color', [0.7 0.7 0.7]); plot([-18 -18], yLims, '--', 'Color', [0.7 0.7 0.7]); % Orientation lines

ylim(yLims); pbaspect([1 1 1]); box off;
xlim([-40 60]);
set(gca, 'XTick', [-40:20:60]);
set(gca, 'YTick', [0.2 0.3 0.4 0.5]);
set(gca, 'TickDir', 'out');

% Save figure
set(gcFig3, 'PaperPosition', [0 0 8 8]);
set(gcFig3, 'PaperSize', [8 8]);
saveas(gcFig3, fullfile(resultsPath, ['IlluminationSpectraDataset_Analysis_Figure3.pdf']), 'pdf');
close(gcFig3);

fprintf('Done.');
