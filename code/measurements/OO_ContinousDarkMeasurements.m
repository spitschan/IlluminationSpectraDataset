function OO_IntervalMeasurements(spectrometerObj, wrapper, spectrometerIndex, channelIndex)
% OO_Scope(spectrometerObj, wrapper, spectrometerIndex, channelIndex)
%
% Scope function.
%
% 6/26/14     ms      Written.

% Continously take and display spectra
% Set integration time.
minIntegrationTime = 1000; % microseconds;
maxIntegrationTime = 60*1000*1000; % microsecond = 60 seconds.
integrationTime = minIntegrationTime; % Start with minimum integration time
invoke(spectrometerObj, 'setIntegrationTime', spectrometerIndex, channelIndex, integrationTime);

% Enable correct for detector non-linearity.
setCorrectForDetectorNonlinearityFlag = false;
invoke(spectrometerObj, 'setCorrectForDetectorNonlinearity', spectrometerIndex, channelIndex, setCorrectForDetectorNonlinearityFlag);

% Enable correct for electrical dark.
setCorrectForElectricalDarkFlag = false;
invoke(spectrometerObj, 'setCorrectForElectricalDark', spectrometerIndex, channelIndex, setCorrectForElectricalDarkFlag);

% Discrete settings for integration time are preferrable. They are:
%         1,000 microseconds = 1 milisecond
%         2,000 microseconds = 2 miliseconds
%         4,000 microseconds = 4 miliseconds
%         8,000 microseconds = 8 miliseconds
%        16,000 microseconds = 16 miliseconds
%        32,000 microseconds = 32 miliseconds
%        64,000 microseconds = 64 miliseconds
%       128,000 microseconds = 128 miliseconds
%       256,000 microseconds = 256 miliseconds
%       512,000 microseconds = 512 miliseconds
%     1,000,000 microseconds = 1 second
%     2,000,000 microseconds = 2 seconds
%     4,000,000 microseconds = 4 seconds
%     8,000,000 microseconds = 8 seconds
%    16,000,000 microseconds = 16 seconds
%    32,000,000 microseconds = 32 seconds
%    60,000,000 microseconds = 60 seconds
theIntegrationTimes = [1000 2000 4000 8000 16000 32000 64000 128000 256000 ...
    512000 1000000 2000000 4000000 8000000 16000000 32000000 60000000];

% Generate a sound.
t = linspace(0, 1, 10000);
y = sin(330*2*pi*t);

% Generate meta data
% Set up position information. This has to be done manually typically.
metaData.position.name.full         = 'Cherry Springs State Park/Cabin';
metaData.position.name.abbrev       = 'CSSP/Cabin';
metaData.position.latDecimalDegrees = [];
metaData.position.lonDecimalDegrees = [];
metaData.position.elevationMeters   = [];

% Get and store spectrometer information
[numOfSpectrometers, spectrometerName, spectrometerSerialNumber] = OO_GetInfo(spectrometerObj, spectrometerIndex, false);
metaData.spectrometerInfo.name.full          = 'OceanOptics USB2000+';
metaData.spectrometerInfo.name.abbrev        = 'OO_USB2000+';
appendType = GetWithDefault('Which spectrometer [a = high sensitivity, b = low sensitivity?', 'a');
metaData.spectrometerInfo.name.abbrev = [metaData.spectrometerInfo.name.abbrev appendType];
metaData.spectrometerInfo.serialNumber       = spectrometerSerialNumber;
metaData.spectrometerInfo.modelType          = spectrometerName;
metaData.spectrometerInfo.numOfSpectrometers = numOfSpectrometers;
metaData.spectrometerInfo.calibrationDate    = '';
metaData.filter = '';


% Some house keeping.
theIndex = 1;
keepLooping = true;

t0 = GetSecs;
abs0 = t0;
while (keepLooping)
    for i = 1:length(theIntegrationTimes);
        theIntegrationTime = theIntegrationTimes(i);

        % Set integration time   
        fprintf('*** Integration time: %g ... ', theIntegrationTime);
        invoke(spectrometerObj, 'setIntegrationTime', spectrometerIndex, channelIndex, theIntegrationTime);
        fprintf('set.\n');
        
        % Figure out how many spectra we take. We take 10 measurements below and including 4
        % seconds, 6 for 8 seconds, 3 for 16 seconds and 1 for 32 and 60
        % seconds;
        nAverage = 1;
        
            fileName = [datestr(now, 30) '_' metaData.spectrometerInfo.name.abbrev '_dark_' num2str(theIntegrationTime/1000) 'ms'];
            OO_MeasureAndSave(spectrometerObj, wrapper, spectrometerIndex, channelIndex, theIntegrationTime, nAverage, false, fileName, metaData);
            
            % Play a sound
            sound(y, 40000);

    end
end