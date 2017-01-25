%% IlluminationSpectraDataset_FullMoon
%
% Calculates the illuminance of the full moon data in the Spitschan et al.
% (2016) data set.
%
% Requires: Psychtoolbox-3
%
% 1/25/2017     ms      Written.

% Housekeeping
close all; clearvars;

% Determine some paths for set up
basePath = fileparts(mfilename('fullpath'));
dataPath = fullfile(basePath, '../data');
resultsPath = fullfile(basePath, 'results', 'fullmoon');
if ~isdir(resultsPath);
    mkdir(resultsPath);
end
dataFile = 'CSSP_spectra.csv';

% Load the data from the CSV files
rawData = csvread(fullfile(dataPath, dataFile), 1, 0); % Skip the first line, which contains a string

% Extract the data
fileID = fopen(fullfile(dataPath, dataFile));
C = textscan(fileID,'%s', 'Delimiter', ',');

% Extract the data from the matrix
wls = (280:840)';             % Define the wavelength sampling between 150 and 1100 nm with 1 nm spacing
spd = rawData(5:end, :);    % Extract the spds
dateID = rawData(1, :);                     % Extract the date no.
solarElevationDeg = rawData(2, :);          % Extract the solar elevation
lunarElevationDeg = rawData(3, :);           % Extract the lunar elevation
lunarFractionIlluminated = rawData(4, :);   % Extract the fraction illuminated

% Separately obtain the date stamps from the first line
fid = fopen(fullfile(dataPath, dataFile));
tmp = fread(fid, '*char');
fclose(fid);
entries = regexp(tmp', '\n', 'split');
dates = regexp(entries{1}, ',', 'split');
data.spd = spd;
data.solarElevationDeg = solarElevationDeg;
data.lunarElevationDeg = lunarElevationDeg;
data.lunarFractionIlluminated = lunarFractionIlluminated;

%% Pull out the full moon spectra
theIdx = find(data.lunarFractionIlluminated  > 0.90 & data.solarElevationDeg < -18);
moonSpds = data.spd(:, theIdx);
moonFraction = lunarFractionIlluminated(theIdx);
lunarElevation = data.lunarElevationDeg(theIdx);
solarElevation = data.solarElevationDeg(theIdx);

%% Load v-lambda (from Psychtoolbox-3)
load T_xyz1931;
magicFactor = 683;
T_vLambda = SplineCmf(S_xyz1931,magicFactor*T_xyz1931(2,:),wls);
clear T_xyz1931 S_xyz1931

%% Plot
figure;
plot(lunarElevation, T_vLambda*moonSpds./sind(lunarElevation), '-k'); hold on;
plot(lunarElevation, T_vLambda*moonSpds./sind(lunarElevation), '-sk', 'MarkerFaceColor', 'white'); hold on;

set(gca, 'TickDir', 'out');
pbaspect([1 1 1]); box off;
xlabel('Lunar elevation [deg]');
ylabel('Illuminance [lux]'); title({'Lunar illuminance' 'Full moon, fraction illuminated = 0.99' 'Solar elevation <-18 deg'});

% Save the figure
set(gcf, 'PaperPosition', [0 0 5 5]);
set(gcf, 'PaperSize', [5 5]);
saveas(gcf, fullfile(resultsPath, 'Spitschan2016_FullMoon_Illuminance'),'png');

% Save CSV files
M1 = [solarElevation' lunarElevation' moonFraction'];
csvwrite(fullfile(resultsPath, 'Spitschan2016_FullMoon_Aux.csv'), M1);
csvwrite(fullfile(resultsPath, 'Spitschan2016_FullMoon_Spectra.csv'), moonSpds);
csvwrite(fullfile(resultsPath, 'Spitschan2016_FullMoon_Wls.csv'), wls);