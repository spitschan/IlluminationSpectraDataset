function [omniAbsSpectrum wls] = IlluminationSpectraDataset_SpecCal_RawToAbsolute(cal, omniRawSpectrum, omniRawWls, omniIntegrationTime, omniBoardTempInCelsius, filterSpan)
% [omniRelSpectrum,omniWls] = IlluminationSpectraDataset_SpecCal_RawToAbsolute(cal, omniRawSpectrum, omniRawWls, omniIntegrationTime, omniBoardTempInCelsius, filterSpan)
%
% Convert a raw spectrum measured with the omni to a calibrated relative 
% spectrum.  Uses the calibration structure.  Retuned spectrum is splined
% onto evenly spaced wavelengths with the span of cal.commonWls and
% the same number of samples as cal.commonWls.
%
% 4/14/2014 spitschan   Adapted.

%% Correct the wavelengths
omniCorrectedWls = omniRawWls-cal.wlShift;

%% Subtract the dark
darkSpectrum = IlluminationSpectraDataset_SpecCal_GetDarkSpectrum(cal, omniIntegrationTime/1000/1000, omniBoardTempInCelsius);
omniDarkSubtracted = (omniRawSpectrum-darkSpectrum) / omniIntegrationTime;

%% Interpolate to 1 nm spacing
omniCorrectedSpectrum = interp1(omniCorrectedWls, omniDarkSubtracted, cal.wls);

%% Filter the raw spectrum
omniFilteredSpectrum = medfilt1(omniCorrectedSpectrum, filterSpan);

% Get correction factor, end up in power per wavelength band.
omniAbsSpectrum = omniFilteredSpectrum.*cal.correctionFactor;

% Convert to radiance by multiplying with pi.
omniAbsSpectrum = pi*omniAbsSpectrum;

% Return the wls as well
wls = cal.wls;