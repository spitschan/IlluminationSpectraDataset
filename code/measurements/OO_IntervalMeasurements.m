function OO_IntervalMeasurements(spectrometerObj, wrapper, spectrometerIndex, channelIndex, intervalInSeconds)
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

% Generate meta data
% Set up position information. This has to be done manually typically.
metaData.position.name.full         = 'Cherry Springs State Park';
metaData.position.name.abbrev       = 'CSSP';
metaData.position.latDecimalDegrees = [];
metaData.position.lonDecimalDegrees = [];
metaData.position.elevationMeters   = [];

% No filter used.
metaData.filter = '';

% Get and store spectrometer information
[numOfSpectrometers, spectrometerName, spectrometerSerialNumber] = OO_GetInfo(spectrometerObj, spectrometerIndex, false);
metaData.spectrometerInfo.name.full          = 'OceanOptics USB2000+';
metaData.spectrometerInfo.name.abbrev        = 'OO_USB2000+';
metaData.spectrometerInfo.serialNumber       = spectrometerSerialNumber;
metaData.spectrometerInfo.modelType          = spectrometerName;
metaData.spectrometerInfo.numOfSpectrometers = numOfSpectrometers;
metaData.spectrometerInfo.calibrationDate    = '';

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

% Set up keyboard handling
% Wait for a character keypress.
FlushEvents;

% Display some information
disp('Press q to exit');

% Some house keeping.
theIndex = 1;
keepLooping = true;

t0 = GetSecs;
abs0 = t0;
while (keepLooping)
    if integrationTime <  8000000
        [wls, spd] = OO_GetSpectrum(spectrometerObj, spectrometerIndex, channelIndex);
        figure(gcf);
        maxCount = 65024;
        plot(wls, spd); hold on;
        plot([150 900], [maxCount*0.50 maxCount*0.50], '-k');
        plot([150 900], [maxCount*0.85 maxCount*0.85], '--k');
        plot([150 900], [maxCount*0.85*0.50 maxCount*0.85*0.50], '--k');
        plot([150 900], [maxCount maxCount], '-k'); hold off;
        xlabel('Wavelength [nm]'); ylabel('Intensity count'); pbaspect([1 1 1]); xlim([150 900]); ylim([-2000 maxCount+3000]);
        title({[num2str(integrationTime) ' us'] ; [num2str(integrationTime/1000) ' ms'] ; [num2str(integrationTime/(1000*1000)) ' s']});
        whitebg('white');
        drawnow;
    end
    
    % Check if a key press is available
    if CharAvail
        % Get the key
        theKey = GetChar(false, true);
        
        if (theKey == 'q')
            keepLooping = false;
            fprintf('> Final integration time: %f', integrationTime);
        end
        
        % Increase/decrease integration time, making 10 the smallest)
        if (theKey == '1')
            theIndex = theIndex - 1;
        end
        
        if (theKey == '2')
            theIndex = theIndex + 1;
        end
        
        % Check that we are not out of bounds
        if theIndex < 1
            theIndex = 1;
        end
        
        if theIndex > length(theIntegrationTimes)
            theIndex = length(theIntegrationTimes);
        end
        integrationTime = theIntegrationTimes(theIndex);
    end
    % Set integration time
    invoke(spectrometerObj, 'setIntegrationTime', spectrometerIndex, channelIndex, integrationTime);
    
    % Figure out how many spectra we take. We take 10 measurements below and including 4
    % seconds, 6 for 8 seconds, 3 for 16 seconds and 1 for 32 and 60
    % seconds;
    if integrationTime <=  4000000
        nAverage = 10;
    elseif integrationTime == 8000000
        nAverage = 5;
    elseif integrationTime == 16000000
        nAverage = 3;
    else
        nAverage = 1;
    end
    
    t1 = GetSecs;
    if t1 - t0 > intervalInSeconds
        t0 = GetSecs;
        OO_MeasureAndSave(spectrometerObj, wrapper, spectrometerIndex, channelIndex, integrationTime, nAverage, false, [], metaData);
        
        % Play a sound
        sound(y, 40000);
        
    end
    
end
end