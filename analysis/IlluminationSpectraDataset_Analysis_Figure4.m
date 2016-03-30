%% IlluminationSpectraDataset_Analysis_Figure4
%
% Generates Figure 4.
%
% 12/10/2015  spitschan   Wrote it.

fprintf('\n>> Making Figure 4...');

%% sRGB rendering
% Load XYZ
load T_xyz1931
Snew = WlsToS((380:1:780)');
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,Snew);

for locIndex = [1 2]
    idx = data(locIndex).theIdx;
    spdsTmp = [data(locIndex).spd(:, idx)];
    solElsTmp = data(locIndex).solarElevationDeg(:, idx);
    
    % Bin
    for i = 1:nBands
        [~, ~, ~, ~, data(locIndex).spdsL2Binned(i, :)] = bindataflex(solElsTmp, spdsTmp(i, :), solarElevationLevels);
    end
    
    % Spline to be in 380:1:780 spacing
    data(locIndex).spdsL2Splined = NaN*ones(Snew(3), size(data(locIndex).spdsL2Binned, 2));
    for s = 1:size(data(locIndex).spdsL2Binned, 2)
        try
            data(locIndex).spdsL2Splined(:, s) = SplineSpd(S, data(locIndex).spdsL2Binned(:, s), Snew);
        end
    end
    
    for s = 1:size(data(locIndex).spdsL2Splined, 2)
        data(locIndex).spdsL2NormScalar(s) = norm(data(locIndex).spdsL2Splined(:, s));
        data(locIndex).spdsL2Norm(:, s) = data(locIndex).spdsL2Splined(:, s)/data(locIndex).spdsL2NormScalar(s);
    end
end

% Scale the spectra
scalars = [data(1).spdsL2NormScalar ; data(2).spdsL2NormScalar];
minScalar = min(log10(scalars(:)));
tmp = (-minScalar + log10(scalars));
maxScalar = max(tmp(:));
scaled = 0.3+0.7*(-minScalar + log10(scalars))/maxScalar;

for locIndex = [1 2]
    for s = 1:size(data(locIndex).spdsL2Norm, 2);
        data(locIndex).spdsL2NormRescaled(:, s) = scaled(locIndex, s)*data(locIndex).spdsL2Norm(:, s);
    end
end

theXYZ1 = T_xyz*data(1).spdsL2NormRescaled;
theXYZ2 = T_xyz*data(2).spdsL2NormRescaled;

% Convert to linear SRGB
SRGBPrimary1 = XYZToSRGBPrimary(theXYZ1);%
SRGBPrimary2 = XYZToSRGBPrimary(theXYZ2);%

% Scale to range 0-1.  This should be done
% with a single scale factor for everything
% that gets converted.
normalizingSRGB = [SRGBPrimary1 SRGBPrimary2];
SRGBPrimaryNorm1 = SRGBPrimary1/max(normalizingSRGB(:));
SRGBPrimaryNorm2 = SRGBPrimary2/max(normalizingSRGB(:));
SRGB1 = SRGBGammaCorrect(SRGBPrimaryNorm1,false)/255;
SRGB2 = SRGBGammaCorrect(SRGBPrimaryNorm2,false)/255;

% Make it into a square image patch and show it
clear theImage1 theImage2
for i = 1:length(solarElevationLevels)
    theImage1(1,i,:) = SRGB1(:, i);
    theImage2(1,i,:) = SRGB2(:, i);
end
%theImage1 = imresize(theImage1,255,'nearest');
imwrite(theImage1, fullfile(resultsPath, 'IlluminationSpectraDataset_Analysis_Figure4a.png'));

%theImage2 = imresize(theImage2,255,'nearest');
imwrite(theImage2, fullfile(resultsPath, 'IlluminationSpectraDataset_Analysis_Figure4b.png'));

fprintf('Done.');