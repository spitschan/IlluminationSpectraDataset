%% Contents
%
% helpers/ - Various helper functions used in IlluminationSpectraDataset
%	bindataflex - function for data binning
%   corr_window - Calculates a correlation in a given window
%
% helpers/thirdparty - Third party tools
%   errorbarX - From BrainardLabToolbox [https://github.com/DavidBrainard/BrainardLabToolbox]
%	hex2rgb - From MATLAB file exchange [http://www.mathworks.com/matlabcentral/fileexchange/45727-hex2rgb]
%   GetInput - From BrainardLabToolbox [https://github.com/DavidBrainard/BrainardLabToolbox]
% 	mseb - From MATLAB file exchange [http://www.mathworks.com/matlabcentral/fileexchange/47950-mseb-x-y-errbar-lineprops-transparent-]
%   strsplit - From BrainardLabToolbox [https://github.com/DavidBrainard/BrainardLabToolbox]
%   subtightplot - From MATLAB file exchange [http://www.mathworks.com/matlabcentral/fileexchange/39664-subtightplot]
%
% measurements/ - Various programs used to control the Ocean Optics spectrometers
%   OO_ContinousDarkMeasurements - Takes continuous dark measurements
%   OO_Demo - Demo measurement program
%   OO_GetBoardTemperature - Get board temperature
%   OO_GetInfo - Assemble meta data
%   OO_GetSpectrum - Obtain the spectrum
%   OO_Init - Initialize the spectrometer
%   OO_IntervalMeasurements - Do timed measurements
%   OO_MeasureAndSave - Measure and save the spectrum
%   OO_SaveMeasAscii - Save to ASCII file
%   OO_SaveMeasMATLAB - Save to MATLAB file
%   OO_Scope - Online live demo