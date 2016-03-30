%% IlluminationSpectraDataset_Analysis_FigureS4
%
% Generates Figure S4.
%
% 12/10/2015  spitschan   Wrote it.

fprintf('\n>> Making Figure S4...');
gcFigS4 = figure;

%% Get the spds
spd_DiCarloWandell = IlluminationSpectraDataset_Analysis_GetDataset('DiCarlo&Wandell', wls);
spd_Granada = IlluminationSpectraDataset_Analysis_GetDataset('Granada', wls);
spd_new = data_merged.spd(:, (data_merged.regime == 5));

%% Get the basis functions
B_CIE = IlluminationSpectraDataset_Analysis_GetModel('CIE', wls, true);
B_Granada = IlluminationSpectraDataset_Analysis_GetModel('Granada', wls, true);

%% Fit to DiCarlo & Wandell data set
[~, ~, ~, gof_DiCarloWandell_CIE] = IlluminationSpectraDataset_Analysis_FitModel(wls, spd_DiCarloWandell, wls, B_CIE, 380, 780);
[~, ~, ~, gof_DiCarloWandell_Granada] = IlluminationSpectraDataset_Analysis_FitModel(wls, spd_DiCarloWandell, wls, B_Granada, 380, 780);
[~, ~, ~, gof_DiCarloWandell_CIE3R] = IlluminationSpectraDataset_Analysis_FitModel(wls, spd_DiCarloWandell, (360:1:830)', B_CIE3R, 380, 780);
[~, ~, ~, gof_DiCarloWandell_CIE3C] = IlluminationSpectraDataset_Analysis_FitModel(wls, spd_DiCarloWandell, (360:1:830)', B_CIE3C, 380, 780);

% Summarize
gof_DiCarloWandell = [gof_DiCarloWandell_CIE ; gof_DiCarloWandell_Granada ; gof_DiCarloWandell_CIE3R ; gof_DiCarloWandell_CIE3C];
gof_DiCarloWandell_mean = mean(gof_DiCarloWandell, 2);
gof_DiCarloWandell_std = std(gof_DiCarloWandell, [], 2);

%% Fit to Granada data set
[~, ~, ~, gof_Granada_CIE] = IlluminationSpectraDataset_Analysis_FitModel(wls, spd_Granada, wls, B_CIE, 380, 780);
[~, ~, ~, gof_Granada_Granada] = IlluminationSpectraDataset_Analysis_FitModel(wls, spd_Granada, wls, B_Granada, 380, 780);
[~, ~, ~, gof_Granada_CIE3R] = IlluminationSpectraDataset_Analysis_FitModel(wls, spd_Granada, (360:1:830)', B_CIE3R, 380, 780);
[~, ~, ~, gof_Granada_CIE3C] = IlluminationSpectraDataset_Analysis_FitModel(wls, spd_Granada, (360:1:830)', B_CIE3C, 380, 780);

% Summarize
gof_Granada = [gof_Granada_CIE ; gof_Granada_Granada ; gof_Granada_CIE3R ; gof_Granada_CIE3C];
gof_Granada_mean = mean(gof_Granada, 2);
gof_Granada_std = std(gof_Granada, [], 2);

%% Fit to our dataset
[~, ~, ~, gof_new_CIE] = IlluminationSpectraDataset_Analysis_FitModel(wls, spd_new, wls, B_CIE, 380, 780);
[~, ~, ~, gof_new_Granada] = IlluminationSpectraDataset_Analysis_FitModel(wls, spd_new, wls, B_Granada, 380, 780);
[~, ~, ~, gof_new_CIE3R] = IlluminationSpectraDataset_Analysis_FitModel(wls, spd_new, (360:1:830)', B_CIE3R, 380, 780);
[~, ~, ~, gof_new_CIE3C] = IlluminationSpectraDataset_Analysis_FitModel(wls, spd_new, (360:1:830)', B_CIE3C, 380, 780);

% Summarize
gof_new = [gof_new_CIE ; gof_new_Granada ; gof_new_CIE3R ; gof_new_CIE3C];
gof_new_mean = mean(gof_new, 2);
gof_new_std = std(gof_new, [],  2);

%% Plot as bar graph
M = [gof_Granada_mean gof_DiCarloWandell_mean gof_new_mean];
Ms = [gof_Granada_std gof_DiCarloWandell_std gof_new_std];
offsetX = linspace(-0.25, 0.25, 4);
theColors = [0.4 0.4 0.4 ; 0.3 0.8 0.3 ; 205/255 35/255 31/255 ; 51/255 105/255 180/255];

for i = 1:3
    for j = 1:4
       h(j) = bar(i+offsetX(j),  M(j, i), 0.15, 'FaceColor', theColors(j, :), 'EdgeColor', 'w'); hold on;
       errorbar(i+offsetX(j),  M(j, i), Ms(j, i), '-k');
    end
end
plot([0 4], [1 1], '-', 'Color', [0 0 0]);  
xlim([0 4]);
ylim([0 1.1]);
xlabel('Mean R2');

set(gca, 'TickDir', 'out'); box off;
set(gca, 'XTick', [1 2 3]); set(gca, 'XTickLabel', {'Granada', 'D&W', 'This dataset'});
pbaspect([1 1 1]);
legend(h, 'CIE', 'Granada', 'CIE3R', 'CIE3C', 'Location', 'EastOutside'); legend boxoff;
title('Daylight fits');

% Save the graph
set(gcFigS4, 'PaperPosition', [0 0 5 5]);
set(gcFigS4, 'PaperSize', [5 5]);
saveas(gcFigS4, fullfile(resultsPath, ['IlluminationSpectraDataset_Analysis_FigureS4.pdf']), 'pdf');
close(gcFigS4);

fprintf('Done.');