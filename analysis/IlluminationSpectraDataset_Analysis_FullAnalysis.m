function IlluminationSpectraDataset_Analysis_FullAnalysis
%% IlluminationSpectraDataset_Analysis_FullAnalysis
%
% Performs analysis on data set
%
% 7/3/2014    spitschan   Written.
% 3/31/2015   spitschan   Commented.

% Housekeeping
close all; clearvars;

% Determine some paths for set up
basePath = fileparts(mfilename('fullpath'));
dataPath = fullfile(basePath, '../data');
resultsPath = fullfile(basePath, 'results');

% Now check what we have in the measurements folder. The hierarchy is typically:
% data/<location>/<date>
theLocations(1).name = 'CSSP';
theLocations(2).name = 'DRL';

% Iterate over the locations
%   1 = Rural (Cherry Springs State Park)
%   2 = City (Philadelphia, PA)
for locIndex = [1 2];
    % Load the data from the CSV files
    rawData = csvread(fullfile(dataPath, [theLocations(locIndex).name '_spectra.csv']), 1, 0); % Skip the first line, which contains a string
    
    % Extract the data
    fileID = fopen(fullfile(dataPath, [theLocations(locIndex).name '_spectra.csv']));
    C = textscan(fileID,'%s', 'Delimiter', ',');
    
    % Extract the data from the matrix
    wls = (280:840)';             % Define the wavelength sampling between 150 and 1100 nm with 1 nm spacing
    spd = rawData(5:end, :);    % Extract the spds
    dateID = rawData(1, :);                     % Extract the date no.
    solarElevationDeg = rawData(2, :);          % Extract the solar elevation
    lunarElevationDeg = rawData(3, :);           % Extract the lunar elevation
    lunarFractionIlluminated = rawData(4, :);   % Extract the fraction illuminated
    
    % Separately obtain the date stamps from the first line
    fid = fopen(fullfile(dataPath, [theLocations(locIndex).name '_spectra.csv']));
    tmp = fread(fid, '*char');
    fclose(fid);
    entries = regexp(tmp', '\n', 'split');
    dates = regexp(entries{1}, ',', 'split');
    
    % Display information about the data set
    uniqueDates = unique(dateID);
    fprintf('\n>> Location %g', locIndex);
    for i = 1:length(uniqueDates)
        tmp = find(dateID == uniqueDates(i));
        fprintf('\n>> Date %g: %g spectra', uniqueDates(i), length(tmp));
        fprintf('\n>> Date %g: Range solar elevation [°]: %.2f - %.2f', uniqueDates(i), min(solarElevationDeg(tmp)), max(solarElevationDeg(tmp)));
        fprintf('\n>> Date %g: Range lunar elevation [°]: %.2f - %.2f', uniqueDates(i), min(lunarElevationDeg(tmp)), max(lunarElevationDeg(tmp)));
        fprintf('\n>> Date %g: Range lunar fraction illuminated [prop.]: %.2f - %.2f', uniqueDates(i), min(lunarFractionIlluminated(tmp)), max(lunarFractionIlluminated(tmp)));
        fprintf('\n>> Date %g: Start/end time: %s / %s', uniqueDates(i), dates{tmp(1)}, dates{tmp(end)});
    end
    fprintf('\n');
        
    % Save the information out
    switch locIndex
        case 1
            locationName = 'Rural';
            fid = fopen(fullfile(resultsPath, ['IlluminationSpectraDataset_Analysis_Table1.csv']), 'w');
        case 2
            locationName = 'City';
            fid = fopen(fullfile(resultsPath, ['IlluminationSpectraDataset_Analysis_Table1.csv']), 'a');
    end
    for i = 1:length(uniqueDates)
        tmp = find(dateID == uniqueDates(i));
        fprintf(fid, '%s,%s,%s,%.2f - %.2f,%.2f - %.2f,%.2f - %.2f,%g\n', locationName, ...
           dates{tmp(1)}, dates{tmp(end)}, min(solarElevationDeg(tmp)), max(solarElevationDeg(tmp)), ...
           min(lunarElevationDeg(tmp)), max(lunarElevationDeg(tmp)), ...
           min(lunarFractionIlluminated(tmp)), max(lunarFractionIlluminated(tmp)), ...
           length(tmp));
    end
    fclose(fid);
    
    %% Filter the data
    spdraw = spd;
    
    % Set up the solar elevation levels for binning.
    solarElevationLevels = [-30:1:72];
    
    data(locIndex).spd = spd;
    data(locIndex).spdraw = spdraw;
    data(locIndex).solarElevationDeg = solarElevationDeg;
    data(locIndex).lunarElevationDeg = lunarElevationDeg;
    data(locIndex).lunarFractionIlluminated = lunarFractionIlluminated;
    switch locIndex
        case 1
            data(locIndex).lunarFractionIlluminatedCriterion = 0.3;
        case 2
            data(locIndex).lunarFractionIlluminatedCriterion = 1;
    end
    
    data(locIndex).date = dateID;    
    data(locIndex).theIdxNight = find(data(locIndex).solarElevationDeg < -18 & data(locIndex).lunarFractionIlluminated < data(locIndex).lunarFractionIlluminatedCriterion); % Astronomical;
    data(locIndex).theIdxTwilightAstronomical = find(data(locIndex).solarElevationDeg >= -18 & data(locIndex).solarElevationDeg < -12 & data(locIndex).lunarFractionIlluminated < data(locIndex).lunarFractionIlluminatedCriterion); % Astronomical
    data(locIndex).theIdxTwilightNautical = find(data(locIndex).solarElevationDeg >= -12 & data(locIndex).solarElevationDeg < -6 & data(locIndex).lunarFractionIlluminated < data(locIndex).lunarFractionIlluminatedCriterion); % Nautical
    data(locIndex).theIdxTwilightCivic = find(data(locIndex).solarElevationDeg >= -6 & data(locIndex).solarElevationDeg < 0 & data(locIndex).lunarFractionIlluminated < data(locIndex).lunarFractionIlluminatedCriterion); % Civic
    data(locIndex).theIdxDay = find(data(locIndex).solarElevationDeg >= 0 & data(locIndex).lunarFractionIlluminated < data(locIndex).lunarFractionIlluminatedCriterion); % Day
    data(locIndex).theIdx = [data(locIndex).theIdxNight data(locIndex).theIdxTwilightAstronomical data(locIndex).theIdxTwilightNautical data(locIndex).theIdxTwilightCivic data(locIndex).theIdxDay];
end % Iterate over locations

solarElevationLevels = -30:1:75;
IlluminationSpectraDataset_Analysis_Figure1; % Make Fig. 1

solarElevationLevels = -30:2:75;
IlluminationSpectraDataset_Analysis_Figure2; % Make Fig. 2
IlluminationSpectraDataset_Analysis_Figure3; % Make Fig. 3
IlluminationSpectraDataset_Analysis_Figure4; % Make Fig. 4
IlluminationSpectraDataset_Analysis_Figure5; % Make Fig. 5
IlluminationSpectraDataset_Analysis_Figure6; % Make Fig. 6

IlluminationSpectraDataset_Analysis_FigureS1; % Make Fig. S1
IlluminationSpectraDataset_Analysis_FigureS2; % Make Fig. S2
IlluminationSpectraDataset_Analysis_FigureS3; % Make Fig. S3
IlluminationSpectraDataset_Analysis_FigureS4; % Make Fig. S4