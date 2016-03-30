function OO_Scope(spectrometerObj, spectrometerIndex, channelIndex)
% OO_Scope(spectrometerObj, spectrometerIndex, channelIndex)
%
% Scope function.
%
% 6/26/14     ms      Written.

% Continously take and display spectra
% Set integration time.
integrationTime = 100;
invoke(spectrometerObj, 'setIntegrationTime', spectrometerIndex, channelIndex, integrationTime);

% Enable correct for detector non-linearity.
setCorrectForDetectorNonlinearityFlag = falee;
invoke(spectrometerObj, 'setCorrectForDetectorNonlinearity', spectrometerIndex, channelIndex, setCorrectForDetectorNonlinearityFlag);

% Enable correct for electrical dark.
setCorrectForElectricalDarkFlag = false;
invoke(spectrometerObj, 'setCorrectForElectricalDark', spectrometerIndex, channelIndex, setCorrectForElectricalDarkFlag);

% Set up the step sizes with which the integration time can be changed.
fineStepSize = 1000;
mediumStepSize = 10000;
coarseStepSize = 100000;

% Set up keyboard handling
% Wait for a character keypress.
FlushEvents;

% Display some information
disp('Press q to exit');

keepLooping = true;
while (keepLooping)
    
    if CharAvail
        % Get the key
        theKey = GetChar(false, true);
        
        if (theKey == 'q')
            keepLooping = false;
            fprintf('> Final integration time: %f', integrationTime);
        end
        
        % Increase/decrease integration time, making 10 the smallest)
        if (theKey == '1')
            integrationTime = max(integrationTime-coarseStepSize, 1);
        end
        
        if (theKey == '2')
            integrationTime = max(integrationTime+coarseStepSize, 1);
        end
        
        if (theKey == '3')
            integrationTime = max(integrationTime-mediumStepSize, 1);
        end
        
        if (theKey == '4')
            integrationTime = max(integrationTime+mediumStepSize, 1);
        end
        
        if (theKey == '5')
            integrationTime = max(integrationTime-fineStepSize, 1);
        end
        
        if (theKey == '6')
            integrationTime = max(integrationTime+fineStepSize, 1);
        end
        
        % Cap at 60s.
        if integrationTime > 60000*1000
            integrationTime = 60000*1000;
        end
    end
    maxCount = 65024;
    % Set integration time
    invoke(spectrometerObj, 'setIntegrationTime', spectrometerIndex, channelIndex, integrationTime);
    [wls, spd] = OO_GetSpectrum(spectrometerObj, spectrometerIndex, channelIndex);
    figure(gcf);
    plot(wls, spd); xlabel('Wavelength [nm]'); ylabel('Intensity count'); pbaspect([1 1 1]); ylim([-2000 maxCount+3000]);
    title(['Integration time: ' num2str(integrationTime/1000) ' ms']);
    drawnow;
end
end