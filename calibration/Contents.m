%% Contents
%
% IlluminationSpectraDataset_SpecCal_DarkCalibration -          Function for dark calibration of the spectrometers 
% IlluminationSpectraDataset_SpecCal_GetCalStruct -             Returns calibration struct
% IlluminationSpectraDataset_SpecCal_GetDarkSpectrum -          Gets the dark spectrum for a given spectrometer,
%                                   							temperature and integration time
% IlluminationSpectraDataset_SpecCal_GetWhichMeterFromWls -     Returns which spectrometer is being used
%                                                               given a wavelength vector
% IlluminationSpectraDataset_SpecCal_InitCalFile -              Performs the calibration
% IlluminationSpectraDataset_SpecCal_RawToAbsolute -            Converts raw uncalibrated spectra to
%                                   							absolute irradiance
% IlluminationSpectraDataset_SpecCal_SensitivityCalibration -   Runs the relative sensitivity
%                                   							calibration
% IlluminationSpectraDataset_SpecCal_WavelengthCalibration -    Runs the wavelength calibration