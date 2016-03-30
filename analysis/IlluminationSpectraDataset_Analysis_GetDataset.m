function [spd_dataset wls] = IlluminationSpectraDataset_Analysis_GetDataset(whichDataset, wls)
% [spd_dataset wls] = IlluminationSpectraDataset_Analysis_GetDataset(whichDataset, wls_in)
%
% Returns a dataset of daylight spectra.
%
% If no wavelength are passed into the function with the wls variable, it
% defaults to 150-1100 nm with 1 nm spacing.
%
% 12/10/2015  spitschan   Wrote it.

% Assume default if wls_in is not passed
if isempty(wls)
    wls = (150:1100)';
end

% Set the data path
dataPath = fullfile(fileparts(mfilename('fullpath')), 'data');

% Load and return the model
switch whichDataset
    case 'DiCarlo&Wandell'
        %% DiCarlo & Wandell (2000)
        tmp = load(fullfile(dataPath, 'spd_DiCarloWandell', 'DiCarloDaylightMeas.mat'));
        wls_DiCarlo = tmp.Wavelength;
        spd_DiCarlo = tmp.SunSpectra;
        
        % Spline it
        spd_dataset = SplineSpd(wls_DiCarlo, spd_DiCarlo, wls);
        
        % Find the limits of the dataset
        wl1 = min(wls_DiCarlo);
        wl2 = max(wls_DiCarlo);
        
        % Set everything outside of these limits to NaN
        spd_dataset(1:find(wls == wl1)-1, :) = NaN;
        spd_dataset(find(wls == wl2)+1:end, :) = NaN;
    case 'Granada'
        % Granada data set of 2600 spectra
        tmp = load(fullfile(dataPath, 'spd_Granada', 'Granada_daylight_2600_81.mat'));
        S_granada = [380 5 81]; % We know this from the website
        wls_granada = SToWls(S_granada);
        spd_granada = tmp.final;
        
        % Spline it;
        spd_dataset = SplineSpd(wls_granada, spd_granada, wls);
        
        % Find the limits of the Granada model
        wl1 = min(wls_granada);
        wl2 = max(wls_granada);
        
        % Set everything outside of these limits to NaN
        spd_dataset(1:find(wls == wl1)-1, :) = NaN;
        spd_dataset(find(wls == wl2)+1:end, :) = NaN;
end