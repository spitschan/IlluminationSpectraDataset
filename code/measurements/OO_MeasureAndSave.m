function OO_MeasureAndSave(spectrometerObj, wrapper, spectrometerIndex, channelIndex, integrationTime, nAverage, darkFlag, fileName, metaData)
% OO_MeasureAndSave(spectrometerObj, wrapper, spectrometerIndex, channelIndex, integrationTime, nAverage, darkFlag, fileName)
%
% Measured and saves the spectrum. darkFlag determines whether the measurement was a dark measurement (i.e. spectrometer capped or not).
%
% 6/26/2014     ms      Written.

%% Get invoke time
invokeTime = now;
fprintf('* Measurement invoke time: %s\n', datestr(invokeTime, 31));

%% Get path and name of file
[scriptPath, scriptFileName] = fileparts(mfilename('fullpath'));

% Set output path for saving
outPath = fullfile(scriptPath, 'measurements');
if ~isdir(outPath)
    mkdir(outPath);
end

%% Set up low-level parameters
PROMPT_FOR_COMMENTS = false;
SAVE_ASCII          = true;
ASCII_FILE_EXT      = '.csv';
SAVE_CAL_STRUCT     = true;

%% Set up meta-data for the measurement.
% Copy over the low-level parameters
meas.describe.low_level.PROMPT_FOR_COMMENTS     = PROMPT_FOR_COMMENTS;
meas.describe.low_level.SAVE_ASCII              = SAVE_ASCII;
meas.describe.low_level.ASCII_FILE_EXT          = ASCII_FILE_EXT;
meas.describe.low_level.SAVE_CAL_STRUCT         = SAVE_CAL_STRUCT;
meas.describe.spectrometerIndex                 = spectrometerIndex;
meas.describe.channelIndex                      = channelIndex;
meas.describe.darkFlag                          = darkFlag;
meas.describe.sampleID                          = '1L';
meas.describe.filter                            = metaData.filter;

% Set up basic information
meas.describe.userID         = 'Manuel Spitschan / University of Pennsylvania';
meas.describe.scriptPath     = scriptPath;
meas.describe.scriptFileName = scriptFileName;
meas.describe.outPath        = outPath;
meas.describe.date           = datestr(invokeTime, 'yyyymmdd');
meas.describe.invokeTime     = datestr(invokeTime, 30);
meas.describe.comment        = '';

% Pull out the meta data
meas.describe.position         = metaData.position;
meas.describe.spectrometerInfo = metaData.spectrometerInfo;

% Store computer information
if ispc
    [~, tmp] = system('ver');
    meas.describe.machineInfo.ver         = tmp;
end
if isunix
    [~, tmp] = system('uname -a');
    meas.describe.machineInfo.ver         = tmp;
end

meas.describe.machineInfo.computerStr = computer;
meas.describe.machineInfo.archStr     = computer('arch');

%% Prompt for a comment before the measurement
if PROMPT_FOR_COMMENTS
    meas.describe.commentPreMeas = GetWithDefault('Comment before measurement', []);
else
    meas.describe.commentPreMeas = [];
end

%% Set integration time.
meas.integrationTime = integrationTime;
invoke(spectrometerObj, 'setIntegrationTime', spectrometerIndex, channelIndex, meas.integrationTime);
meas.integrationTimeSet = invoke(spectrometerObj, 'getIntegrationTime', spectrometerIndex, channelIndex);

%% Enable correct for detector non-linearity.
meas.setCorrectForDetectorNonlinearityFlag = false;
invoke(spectrometerObj, 'setCorrectForDetectorNonlinearity', spectrometerIndex, channelIndex, meas.setCorrectForDetectorNonlinearityFlag);

%% Enable correct for electrical dark.
meas.setCorrectForElectricalDarkFlag = false;
invoke(spectrometerObj, 'setCorrectForElectricalDark', spectrometerIndex, channelIndex, meas.setCorrectForElectricalDarkFlag);

%% Take measurements
% Get a board temperature measurement.
meas.boardTempInCelsius(1) = OO_GetBoardTemperature(wrapper);

% Take the measurement.
for i = 1:nAverage
    [wls, spd(:, i)] = OO_GetSpectrum(spectrometerObj, spectrometerIndex, channelIndex);
end

% Tuck away the data
meas.wls = wls;
meas.spd = spd;

% Get a second board temperature measurement.
meas.boardTempInCelsius(2) = OO_GetBoardTemperature(wrapper);

% Also take a dark measurement with the cap on.
if darkFlag
    input('*** PUT CAP ON THE SPECTROMETER TO AND PRESS ENTER TO PROCEED ***');
    
    % Get a board temperature measurement.
    meas.boardTempInCelsiusDark(1) = OO_GetBoardTemperature(wrapper);
    
    % Take the measurement.
    for i = 1:nAverage
        [wlsDark, spdDark(:, i)] = OO_GetSpectrum(spectrometerObj, spectrometerIndex, channelIndex);
    end
    
    % Tuck away the data
    meas.wlsDark = wlsDark;
    meas.spdDark = spdDark;
    
    % Get a second board temperature measurement.
    meas.boardTempInCelsiusDark(2) = OO_GetBoardTemperature(wrapper);
end

%% Make a plot to look at.
figure(gcf);
maxCount = 65535;
plot(wls, spd); hold on;
plot([150 900], [maxCount*0.50 maxCount*0.50], '-k');
plot([150 900], [maxCount*0.85 maxCount*0.85], '--k');
plot([150 900], [maxCount*0.85*0.50 maxCount*0.85*0.50], '--k');
plot([150 900], [maxCount maxCount], '-k'); hold off;
xlabel('Wavelength [nm]'); ylabel('Intensity count'); pbaspect([1 1 1]); xlim([150 900]); ylim([-2000 maxCount+3000]);
title({[num2str(integrationTime) ' us'] ; [num2str(integrationTime/1000) ' ms'] ; [num2str(integrationTime/(1000*1000)) ' s'] ; ['Board temp: ' num2str(meas.boardTempInCelsius(1)) '/' num2str(meas.boardTempInCelsius(2))]});
whitebg('red');
drawnow;

%% Prompt for a comment after the measurement
if PROMPT_FOR_COMMENTS
    meas.describe.commentPostMeas = GetWithDefault('Comment before measurement', []);
else
    meas.describe.commentPostMeas = [];
end

%% Save out
%meas.describe.saveFileNameASCII = OO_SaveMeasAscii(meas);
meas.describe.saveFileNameASCII = '';
OO_SaveMeasMATLAB(meas, fileName);
end
