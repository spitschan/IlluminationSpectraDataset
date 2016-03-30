function [wls, spd] = OO_GetSpectrum(spectrometerObj, spectrometerIndex, channelIndex)
% [wls, spd] = OO_GetSpectrum(spectrometerObj, spectrometerIndex, channelIndex)
%
% Obtain the spectrum. Wavelengths are returned as first output arg.
%
% 6/26/14     ms      Written.

%% Take measurements
wls = invoke(spectrometerObj, 'getWavelengths', spectrometerIndex, channelIndex);
spd = invoke(spectrometerObj, 'getSpectrum', spectrometerIndex);
end