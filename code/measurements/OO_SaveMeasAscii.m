function saveFileNameASCII = OO_SaveMeasAscii(meas)
% saveFileNameASCII = OO_SaveMeasAscii(meas)
%
% Saves out the measurement as ASCII.
%
% 6/26/14     ms      Written.

try
    % Construct the file name
    fileName = [meas.describe.invokeTime '_s' meas.describe.spectrometerInfo.name.abbrev '_p' meas.describe.position.name.abbrev];
    
    % Check if the file exists. If it does, we do not want to overwrite.
    if exist(fullfile(meas.describe.outPath, [fileName meas.describe.low_level.ASCII_FILE_EXT], 'file')) == 2
        fileName = [fileName '_SAFE' meas.describe.low_level.ASCII_FILE_EXT];
    else
        fileName = [fileName meas.describe.low_level.ASCII_FILE_EXT];
    end
    fid = fopen(fullfile(meas.describe.outPath, fileName), 'w');
    fprintf(fid,'#BEGINHEADER\n');
    fprintf(fid,'# Program: %s\n', meas.describe.scriptFileName);
    fprintf(fid,'# Path: %s\n', meas.describe.scriptPath);
    fprintf(fid,'# Save path: %s\n', meas.describe.outPath);
    fprintf(fid,'# User ID: %s\n', meas.describe.userID);
    fprintf(fid,'# Date: %s\n', meas.describe.date);
    fprintf(fid,'# Invoke time: %s\n', meas.describe.invokeTime);
    fprintf(fid,'#\n');
    fprintf(fid,'# Position name (full): %s\n', meas.describe.position.name.full);
    fprintf(fid,'# Position name (abbrev): %s\n', meas.describe.position.name.abbrev);
    fprintf(fid,'# Lat [decimal degrees]: %f\n', meas.describe.position.latDecimalDegrees);
    fprintf(fid,'# Lon [decimal degrees]: %f\n', meas.describe.position.lonDecimalDegrees);
    fprintf(fid,'# Elevation [meters]: %f\n', meas.describe.position.elevationMeters);
    fprintf(fid,'#\n');
    fprintf(fid,'# Integration time: %f\n', meas.integrationTime);
    fprintf(fid,'# Integration time [read]: %f\n', meas.integrationTimeSet);
    fprintf(fid,'# setCorrectForDetectorNonlinearityFlag: %f\n', num2str(meas.setCorrectForDetectorNonlinearityFlag));
    fprintf(fid,'# setCorrectForElectricalDarkFlag: %f\n', num2str(meas.setCorrectForElectricalDarkFlag));
    fprintf(fid,'# Dark measurement :%f\n', num2str(meas.describe.darkFlag));
    fprintf(fid,'#\n');
    fprintf(fid,'# Board temperature (1) :%f\n', meas.boardTempInCelsius(1));
    fprintf(fid,'# Board temperature (2) :%f\n', meas.boardTempInCelsius(2));
    fprintf(fid,'#\n');
    fprintf(fid,'# Spec name (full): %s\n', meas.describe.spectrometerInfo.name.full);
    fprintf(fid,'# Spec name (abbrev): %s\n', meas.describe.spectrometerInfo.name.abbrev);
    fprintf(fid,'# Spec serial number: %s\n', meas.describe.spectrometerInfo.serialNumber);
    fprintf(fid,'# Spec model type: %s\n', meas.describe.spectrometerInfo.modelType);
    fprintf(fid,'# Spec calibration date: %s\n', meas.describe.spectrometerInfo.calibrationDate);
    fprintf(fid,'#\n');
    fprintf(fid,'# Machine version: %s', meas.describe.machineInfo.ver);
    fprintf(fid,'# Machine computer string: %s\n', meas.describe.machineInfo.computerStr);
    fprintf(fid,'# Machine arch string: %s\n', meas.describe.machineInfo.archStr);
    fprintf(fid,'#\n');
    fprintf(fid,'# PROMPT_FOR_COMMENTS: %s\n', num2str(meas.describe.low_level.PROMPT_FOR_COMMENTS));
    fprintf(fid,'# SAVE_ASCII: %s\n', num2str(meas.describe.low_level.SAVE_ASCII));
    fprintf(fid,'# ASCII_FILE_EXT: %s\n', num2str(meas.describe.low_level.ASCII_FILE_EXT));
    fprintf(fid,'# SAVE_CAL_STRUCT: %s\n', num2str(meas.describe.low_level.SAVE_CAL_STRUCT));
    fprintf(fid,'#\n');
    fprintf(fid,'# Comment (general): %s\n', meas.describe.comment);
    fprintf(fid,'# Comment (pre-measurement): %s\n', meas.describe.commentPreMeas);
    fprintf(fid,'# Comment (post-measurement): %s\n', meas.describe.commentPostMeas);
    fprintf(fid,'#ENDHEADER\n');
    fprintf(fid, '%f\t%f\n', [meas.wls meas.spd]);
    fclose(fid);
    
    % Tuck away where we save the file.
    saveFileNameASCII = fullfile(meas.describe.outPath, fileName);
    fprintf('> Saved ASCII struct to %s\n', saveFileNameASCII);
catch e
    rethrow(e);
end
end