function [numOfSpectrometers, spectrometerName, spectrometerSerialNumber] = OO_GetInfo(spectrometerObj, spectrometerIndex, verbose)
% [numOfSpectrometers, spectrometerName, spectrometerSerialNumber] = OO_GetInfo(spectrometerObj, spectrometerIndex, verbose)
% Get some info on the spectrometers.
%
% 6/26/14     ms      Written.

% Get number of spectrometers connected.
numOfSpectrometers = invoke(spectrometerObj, 'getNumberOfSpectrometersFound');

% Get spectrometer name.
spectrometerName = invoke(spectrometerObj, 'getName', spectrometerIndex);

% Get spectrometer serial number.
spectrometerSerialNumber = invoke(spectrometerObj, 'getSerialNumber', spectrometerIndex);
if verbose
    display(['> Found ' num2str(numOfSpectrometers) ' Ocean Optics spectrometer(s).']);
    display(['> Model Name : ' spectrometerName])
    display(['> Model S/N  : ' spectrometerSerialNumber]);
end
end