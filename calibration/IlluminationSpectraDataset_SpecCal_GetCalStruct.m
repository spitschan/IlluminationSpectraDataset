function cal = IlluminationSpectraDataset_SpecCal_GetCalStruct(whichMeter)
% cal = IlluminationSpectraDataset_SpecCal_GetCalStruct(whichMeter);
%
% Returns a calibration structure
%
% 3/31/2015   spitschan       Commented.

% Figure out the paths
basePath = fileparts(mfilename('fullpath'));
calPath = fullfile(basePath, 'cals');
plotPath = fullfile(calPath, 'plots');

% Make directories if they don't exist
if ~isdir(calPath);
    mkdir(calPath);
end

if ~isdir(plotPath)
    mkdir(plotPath);
end

% Prompt for the calibration type to be used
calTypes = {'a', 'b'};
if isempty(whichMeter)
    while true
        fprintf('- Available calibration types:\n');
        
        for i = 1:length(calTypes)
            fprintf('%d: %s\n', i, calTypes{i});
        end
        
        x = GetInput('Selection', 'number', 1);
        if x >= 1 && x <= length(calTypes)
            break;
        end
    end
    whichMeter = calTypes{x};
end

% Initialize the struct
cal = LoadCalFile(['OO_USB2000+' whichMeter], [], calPath);