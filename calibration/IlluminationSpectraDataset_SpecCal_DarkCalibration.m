function darkSpectrum = IlluminationSpectraDataset_SpecCal_DarkCalibration(cal)
% darkSpectrum = IlluminationSpectraDataset_SpecCal_GetCalStruct(cal)
%
% Performs calibration of the dark spectrum
%
% 3/31/2015   spitschan   Commented.

% Load in all the files
basePath = fullfile(cal.calDataPath, 'dark');
theFiles = dir([basePath '/*USB2000+' cal.whichMeter '*']);

% Define the integration times
theIntegrationTimes = [1000 2000 4000 8000 16000 32000 64000 128000 256000 ...
    512000 1000000 2000000 4000000 8000000 16000000 32000000 60000000];
theIntegrationTimeLabels = {'1 ms', '2 ms', '4 ms', '8 ms', '16 ms', '32 ms', '64 ms', ...
    '128 ms', '256 ms', '512 ms', '1 s', '2 s', '4 s', '8 s', '16 s', '32 s', '60 s'};

% Set up some empty arrays
for i = 1:length(theIntegrationTimes)
    temp{i} = [];
    spd{i} = [];
end

% Load one file for reference
load(fullfile(basePath, theFiles(1).name));
wls = meas.wls-cal.wlShift;

% Load in the files
for f = 1:length(theFiles)
    load(fullfile(basePath, theFiles(f).name));
    
    % Determine the integration time
    ind = find(theIntegrationTimes == meas.integrationTime);
    
    % Throw out any of these '0' spds
    if ~all(meas.spd == 0)
        temp{ind} = [temp{ind} mean(meas.boardTempInCelsius)];
        spd{ind} = [spd{ind} mean(meas.spd, 2)];
    end
end

figure;
for i = 1:length(theIntegrationTimes)
    % Get the temp. range
    minTemp = min(temp{i});
    maxTemp = max(temp{i});
    
    % Sort
    [temp{i}, idx] = sort(temp{i});
    spd{i} = spd{i}(:, idx);
    
    % We only want unique temperature values in the temperature array, so
    % we deal with this here.
    uniqueTemps = unique(temp{i});
    
    spdTemp = [];
    for j = 1:length(uniqueTemps)
        tempIndices = temp{i} == uniqueTemps(j);
        spdTemp(j, :) = mean(spd{i}(:, tempIndices), 2);
    end
    
    % Propagate the temporary spd through
    temp{i} = uniqueTemps;
    spd{i} = spdTemp;
    
    meanSpd{i} = mean(spd{i}, 2)';
    
    % Fit the polynomial
    p = polyfit(temp{i}, meanSpd{i}, 4);
    x1{i} = linspace(minTemp, maxTemp, 1000);
    f1{i} = polyval(p, x1{i});
    
    % Plot
    subtightplot(6, 3, i);
    plot(temp{i}, meanSpd{i}/65536, '.k', 'MarkerSize', 1); hold on
    plot(x1{i}, f1{i}/65536, 'Color', 'r', 'LineWidth', 1);
    
    if i > 14
        xlim([10 60]); ylim([0 1]);
        plot([10 10], [0 1], '-k');
        plot([10 60], [0 0], '-k');
    else
        xlim([10 60]); ylim([0 0.3]);
        plot([10 10], [0 0.25], '-k');
        plot([10 60], [0 0], '-k');
    end
    title(theIntegrationTimeLabels{i});
    pbaspect([1 1 1]);
    set(gca, 'TickDir', 'out'); box off;
    axis off;
end
set(gcf, 'PaperPosition', [0 0 5 15])
set(gcf, 'PaperSize', [5 15]); %Set the paper to have width 5 and height 5.
saveas(gcf, fullfile(cal.plotPath, ['IlluminationSpectraDataset_Calibration_FigureS5_DarkNoise_' cal.whichMeter '.pdf']), 'pdf');
close(gcf);

%% Interpolate the temperature measurements
temps = 10:0.5:50;
for i = 1:length(theIntegrationTimes)
    for j = 1:length(temps);
        darkSpectrum.darkTemperature{i}(j) = temps(j);
        darkSpectrum.darkSpd{i}(:, j) =  interp1(temp{i}, spd{i}, temps(j));
        darkSpectrum.darkMean{i}(j) = interp1(x1{i}, f1{i}, temps(j));
    end
end

% Plot the matrix
figure;
theTempsToPlot = [20:5:40];
theColsRGB = cool(length(theTempsToPlot));
theIntegrationTimesToPlot = [32000 64000 128000 256000 ...
    512000 1000000 2000000 4000000 8000000 16000000 32000000 60000000];
c = 1;
for i = 1:length(theIntegrationTimesToPlot)
    integrationTimeIdx = theIntegrationTimes == theIntegrationTimesToPlot(i);
    for j = 1:length(theTempsToPlot)
        subtightplot(length(theIntegrationTimesToPlot), length(theTempsToPlot), c)
        tempIdx = temps == theTempsToPlot(j);
        plot(wls, darkSpectrum.darkSpd{integrationTimeIdx}(:, tempIdx), 'r', 'LineWidth', 0.5); hold on;
        pbaspect([1 1 1]);
        xlim([125 1100]);
        ylim([-2000 70000]);
        axis off;
        plot([150 1100], [0 0], '-k', 'LineWidth', 0.75);
        if j == 1
            plot([150 150], [0 70000], '-k', 'LineWidth', 0.75);
            plot([125 150], [65536 65536], '-k', 'LineWidth', 0.75)
        end
        c = c+1;
    end
end
set(gcf, 'PaperPosition', [0 0 5 15])
set(gcf, 'PaperSize', [5 15]); %Set the paper to have width 5 and height 5.
saveas(gcf, fullfile(cal.plotPath, ['IlluminationSpectraDataset_Calibration_FigureS5_DarkNoiseSpectra_' cal.whichMeter '.pdf']), 'pdf');
close(gcf);

% Assemble a dark spectrum struct
darkSpectrum.integrationTimes = theIntegrationTimes;
darkSpectrum.wls = wls;