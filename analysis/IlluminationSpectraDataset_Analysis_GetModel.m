function [B_model wls] = IlluminationSpectraDataset_Analysis_GetModel(whichModel, wls, returnNormalized)
% [B_model wls] = IlluminationSpectraDataset_Analysis_GetModel(whichModel, wls_in, returnNormalized)
%
% Returns the CIE model given in Wyszecki & Stiles (1982) or the Granada
% daylight model.
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
switch whichModel
    case 'CIE'
        tmp = load(fullfile(dataPath, 'B_CIE', 'B_CIE.csv'));
        wls_CIE = tmp(:, 1);
        B_CIE = tmp(:, 2:4);
        B_model = SplineSpd(wls_CIE, B_CIE, wls);
        
        % Find the limits of the CIE model
        wl1 = min(wls_CIE);
        wl2 = max(wls_CIE);
        
        % Set everything outside of these limits to NaN
        B_model(1:find(wls == wl1)-1, :) = NaN;
        B_model(find(wls == wl2)+1:end, :) = NaN;
    case 'Granada'
        B_granada = dlmread(fullfile(dataPath, 'B_Granada', 'daylight_380780_norm_eigenvectors.dat'));
        S_granada = [380 5 81];
        wls_granada = SToWls(S_granada);
        B_model = SplineSpd(wls_granada, B_granada, wls);
        
        % Find the limits of the Granada model
        wl1 = min(wls_granada);
        wl2 = max(wls_granada);
        
        % Set everything outside of these limits to NaN
        B_model(1:find(wls == wl1)-1, :) = NaN;
        B_model(find(wls == wl2)+1:end, :) = NaN;
end

% Normalize by the L2 norm if we want it normalized
if returnNormalized
    for i = 1:size(B_model, 2)
        B_model(:, i) = B_model(:, i)/norm(B_model(~isnan(B_model(:, i)), i));
    end
end