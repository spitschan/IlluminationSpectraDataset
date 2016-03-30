function whichMeter = IlluminationSpectraDataset_SpecCal_GetWhichMeterFromWls(wls)
% whichMeter = IlluminationSpectraDataset_SpecCal_GetWhichMeterFromWls(wls)
%
% Returns which meter is used from the wavelength vector
%
% 3/26/2016 spitschan   Added to repository.

if any(wls > 1000)
    whichMeter = 2;
else
    whichMeter = 1;
end