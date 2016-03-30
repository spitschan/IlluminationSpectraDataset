%% Initialize the spectrometer
clear all; close all;
[spectrometerObj, wrapper] = OO_Init();

% Channel and spectrometer indices are always 0.
spectrometerIndex = 0;
channelIndex = 0;

% Get some information.
verboseInfo = true;
OO_GetInfo(spectrometerObj, spectrometerIndex, verboseInfo);

% Open interval measurement mode
intervalInSeconds = 60;
OO_IntervalMeasurements(spectrometerObj, wrapper, spectrometerIndex, channelIndex, intervalInSeconds);
% 
% 
% %% Open the scope
% OO_Scope(spectrometerObj, spectrometerIndex, channelIndex);
% 
% %% Take measurement (dark)
% integrationTime = 200;
% darkFlag = true;
% OO_MeasureAndSave(spectrometerObj, spectrometerIndex, channelIndex, integrationTime, darkFlag)
% 
% %% Take measurement (light)
% darkFlag = false;
% OO_MeasureAndSave(spectrometerObj, spectrometerIndex, channelIndex, integrationTime, darkFlag)
% 
% %% Loop
% while true
% % Take measurement
%     integrationTime = 60*1000*1000;
%     OO_MeasureAndSave(spectrometerObj, spectrometerIndex, channelIndex, integrationTime, darkFlag);
% end
% 
% %% 60 second dark measurement
% integrationTime = 60*1000*1000;
% darkFlag = true;
% OO_MeasureAndSave(spectrometerObj, spectrometerIndex, channelIndex, integrationTime, darkFlag)
% 
% %% Loop
% while true
% % Take measurement
%     integrationTime = 60*1000*1000;
%     OO_MeasureAndSave(spectrometerObj, spectrometerIndex, channelIndex, integrationTime, darkFlag);
% end