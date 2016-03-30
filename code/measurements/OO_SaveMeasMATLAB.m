function OO_SaveMeasMATLAB(meas, fileName)
% OO_SaveMeasMATLAB(meas)
%
% Saves out the 'meas' struct.
%
% 6/26/14     ms      Written.

try
    % Construct the file name
    if isempty(fileName)
        fileName = [meas.describe.invokeTime '_s' meas.describe.spectrometerInfo.name.abbrev '_p' meas.describe.position.name.abbrev];
    end
    
    % Check if the file exists. If it does, we do not want to overwrite.
    if exist(fullfile(meas.describe.outPath, [fileName '.mat'], 'file')) == 2
        fileName = [fileName '_SAFE' '.mat'];
    else
        fileName = [fileName '.mat'];
    end
    saveFileNameMATLAB = fullfile(meas.describe.outPath, fileName);
    
    % Save it, now.
    save(saveFileNameMATLAB, 'meas');
    
    % Report on it
    fprintf('> Saved MATLAB struct to %s\n', saveFileNameMATLAB);
catch e
    rethrow(e);
end
end