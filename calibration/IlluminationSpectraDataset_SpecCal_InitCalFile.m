function IlluminationSpectraDataset_SpecCal_InitCalFile(whichMeter)
% IlluminationSpectraDataset_SpecCal_InitCalFile
%
% Makes a calibration file for the USB2000+ spectrometers.
%
% 3/31/2015   spitschan   Commented

% Housekeeping
close all;

% Figure out the paths
basePath = fileparts(mfilename('fullpath'));
addpath(basePath);
calDataPath = fullfile(basePath, 'calibrationdata');
calPath = fullfile(basePath, 'cals');
plotPath = fullfile(calPath, 'plots');

% Make directories if they don't exist
if ~isdir(calPath);
    mkdir(calPath);
end

if ~isdir(plotPath)
    mkdir(plotPath);
end

%% Initialize the struct
cal = struct;
cal.calPath = calPath;
cal.calDataPath = calDataPath;
cal.plotPath = plotPath;
cal.whichMeter = whichMeter;
cal.wls = (150:1100)';
cal.S = SToWls(cal.wls);

%% Get the correction factor for the wavelengths
cal.wlShift = IlluminationSpectraDataset_SpecCal_WavelengthCalibration(cal);

%% Run dark calibration
cal.dark = IlluminationSpectraDataset_SpecCal_DarkCalibration(cal);

%% Run the absolute calibration
cal = IlluminationSpectraDataset_SpecCal_SensitivityCalibration(cal);

%% Save the calibration
SaveCalFile(cal, ['OO_USB2000+' whichMeter], calPath);
close all;