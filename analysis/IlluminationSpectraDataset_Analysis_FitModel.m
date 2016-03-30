function [spds_norm, spds_fit, w, gof, residuals, scalar] = IlluminationSpectraDataset_Analysis_FitModel(wls, spds, wls_model, B_model, wlRangeStart, wlRangeEnd)
% [spd_norm, spd_fit, w, gof, residuals] = IlluminationSpectraDataset_Analysis_FitModel(wls, spds, wls_model, B_model)
%
% Fits a model using linear regression
%
% 3/26/2016		spitschan	Added to repository

wlRangeStartIdxData = find(wls == wlRangeStart);
wlRangeEndIdxData = find(wls == wlRangeEnd);
validWlRngIdxData = wlRangeStartIdxData:wlRangeEndIdxData;

wlRangeStartIdxModel = find(wls_model == wlRangeStart);
wlRangeEndIdxModel = find(wls_model == wlRangeEnd);
validWlRngIdxModel = wlRangeStartIdxModel:wlRangeEndIdxModel;

% Loop over the spectra and fit them with the model in B
for s = 1:size(spds, 2)
    % Extract the wl range for which we can fit the model
    scalar(s) = norm(spds(validWlRngIdxData, s));
    spds_norm(:, s) = spds(:, s)/scalar(s);                                     % Extract the spd
    w(:, s) = B_model(validWlRngIdxModel, :) \ spds_norm(validWlRngIdxData, s); % Least-squares fit
    spds_fit(:, s) = B_model(validWlRngIdxModel, :)*w(:, s);                    % Reconstruct the spd
    gof(s) = corr(spds_norm(validWlRngIdxData, s), spds_fit(:, s))^2;           % Calculate R^2
    residuals(:, s) = spds_norm(validWlRngIdxData, s) - spds_fit(:, s);         % Calculate residuals
end