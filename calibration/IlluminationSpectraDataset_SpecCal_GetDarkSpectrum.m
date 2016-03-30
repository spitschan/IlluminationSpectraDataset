function darkSpd = IlluminationSpectraDataset_SpecCal_GetDarkSpectrum(cal, reqIntegrationTimeSecs, reqTemperatureDegreesCelsius)
% darkSpd = IlluminationSpectraDataset_SpecCal_GetDarkSpectrum(cal, reqIntegrationTimeSecs, reqTemperatureDegreesCelsius)
%
% Returns a mean dark spectrum appropriate for the integration time and temperature.
%
% 3/31/2015   spitschan   Commented.

% Find the index in the struct that corresponds to the input integration
% time
reqIntegrationTimeMuSecs = 1000*1000*reqIntegrationTimeSecs;
[~, idx] = min(abs(cal.dark.integrationTimes-reqIntegrationTimeMuSecs));

% Throw an error if we do not find that integration time.
if isempty(idx)
   error('Integration time not found'); 
end

% Linearly interpolate. % This is for the mean
darkSpd = interp1(cal.dark.darkTemperature{idx}, cal.dark.darkSpd{idx}', reqTemperatureDegreesCelsius, 'linear', NaN)';