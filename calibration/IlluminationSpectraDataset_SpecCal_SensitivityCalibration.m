function cal = IlluminationSpectraDataset_SpecCal_RelativeSensitivityCalibration(cal)
% cal = IlluminationSpectraDataset_SpecCal_RelativeSensitivityCalibration(cal)
%
% 3/31/15   spitschan   Commented.

switch cal.whichMeter
    case 'a'
        cal.whichMeterName = 'highSensitivity';
        yLimUpper = 0.035/2;
    case 'b'
        cal.whichMeterName = 'lowSensitivity';
        yLimUpper = 3.5/2;
end

% Set the median filter span
medianFilterSpan = 8;

switch cal.whichMeter
    case 'a'
        %%%% Slide projector %%%%
        fileToLoad = fullfile(cal.calDataPath, 'absoluteirradiance', '20141106T110342_sOO_USB2000+a_pLab.mat');
        load(fileToLoad, 'meas');
        meas_USB2000 = meas;
        meas_USB2000.meanSpd = mean(meas.spd, 2);
        wls_USB2000_slide = meas_USB2000.wls-cal.wlShift;
        
        % Dark subtract from the OO measurements
        dark = IlluminationSpectraDataset_SpecCal_GetDarkSpectrum(cal, meas_USB2000.integrationTime/1000/1000, mean(meas_USB2000.boardTempInCelsius));
        spd_USB2000_slide = (meas_USB2000.meanSpd-dark) / meas_USB2000.integrationTime;
        spd_USB2000_slide = medfilt1(interp1(wls_USB2000_slide, spd_USB2000_slide, cal.wls), medianFilterSpan);
    case 'b'
        %%%% Slide projector %%%%
        %% USB2000 measurements
        fileToLoad = fullfile(cal.calDataPath, 'absoluteirradiance', '20141106T110929_sOO_USB2000+b_pLab.mat');
        load(fileToLoad, 'meas');
        meas_USB2000 = meas;
        meas_USB2000.meanSpd = mean(meas.spd, 2);
        wls_USB2000_slide = meas_USB2000.wls-cal.wlShift;
        
        % Dark subtract from the OO measurements
        dark = IlluminationSpectraDataset_SpecCal_GetDarkSpectrum(cal, meas_USB2000.integrationTime/1000/1000, mean(meas_USB2000.boardTempInCelsius));
        spd_USB2000_slide = (meas_USB2000.meanSpd-dark) / meas_USB2000.integrationTime;
        spd_USB2000_slide = medfilt1(interp1(wls_USB2000_slide, spd_USB2000_slide, cal.wls), medianFilterSpan);
        
end

% Load in the data for the SL-x sources
%%%% SL-1 %%%%
whichSource = 'SL1_calRelSpec';
cd(fullfile(cal.calDataPath, 'stellarnet', cal.whichMeterName, whichSource));
% Load the files
theFiles = dir('*.mat');
spds = [];
for f = 1:length(theFiles)
    tmp = load(theFiles(f).name);
    tmpWls = tmp.meas.wls-cal.wlShift;
    spds = [spds mean(tmp.meas.spd, 2)];
    temperatureCelsius(f) = mean(tmp.meas.boardTempInCelsius);
    integrationTime(f) = tmp.meas.integrationTime;
end
% Remove saturated measurements
theIdxToRemove = max(spds, [], 1) == 65535;
spds(:, theIdxToRemove) = [];
temperatureCelsius(theIdxToRemove) = [];
integrationTime(theIdxToRemove) = [];
for s = size(spds, 2)
    darkSpd = IlluminationSpectraDataset_SpecCal_GetDarkSpectrum(cal, (integrationTime(s))/1000/1000, temperatureCelsius(s));
    spd_USB2000_sl1(:, s) = medfilt1(interp1(tmpWls, (spds(:, s)-darkSpd)/integrationTime(s), cal.wls), medianFilterSpan);
end
spd_USB2000_sl1 = mean(spd_USB2000_sl1, 2);

% The SL-1 measurements in all cases are at least 100x brighter than
% the slide projector measurements. We therefore scale the SL-1
% spectrum to match the slide projector meaurement where it is at 10%
% of the max.


%%%% SL-3 %%%%
whichSource = 'SL3_calRelSpec';
cd(fullfile(cal.calDataPath, 'stellarnet', cal.whichMeterName, whichSource));
% Load the files
theFiles = dir('*.mat');
spds = [];
for f = 1:length(theFiles)
    tmp = load(theFiles(f).name);
    tmpWls = tmp.meas.wls-cal.wlShift;
    spds = [spds mean(tmp.meas.spd, 2)];
    temperatureCelsius(f) = mean(tmp.meas.boardTempInCelsius);
    integrationTime(f) = tmp.meas.integrationTime;
end

% Remove saturated measurements
theIdxToRemove = max(spds, [], 1) == 65535;
spds(:, theIdxToRemove) = [];
temperatureCelsius(theIdxToRemove) = [];
integrationTime(theIdxToRemove) = [];
for s = size(spds, 2)
    darkSpd = IlluminationSpectraDataset_SpecCal_GetDarkSpectrum(cal, (integrationTime(s))/1000/1000, temperatureCelsius(s));
    spd_USB2000_sl3(:, s) = medfilt1(interp1(tmpWls, (spds(:, s)-darkSpd)/integrationTime(s), cal.wls), medianFilterSpan);
end
spd_USB2000_sl3 = mean(spd_USB2000_sl3, 2);

switch cal.whichMeter
    case 'a'
        [~, i] = min(abs(spd_USB2000_slide-max(spd_USB2000_slide)/10));
        spd_USB2000_sl1 = spd_USB2000_sl1*(spd_USB2000_slide(i)/spd_USB2000_sl1(i));
    case 'b'
        spd_USB2000_sl3 = spd_USB2000_sl3*(nansum(spd_USB2000_slide)/nansum(spd_USB2000_sl3));
        spd_USB2000_sl1 = spd_USB2000_sl1*(nansum(spd_USB2000_slide)/nansum(spd_USB2000_sl1));
end


%% Load in the calibration files
% SL-1
tmp = importdata(fullfile(cal.calDataPath, 'stellarnet', 'SL1CAL-15081402-ATPLANE.ICD'));
wls_cal_sl1 = tmp.data(:, 1);
spd_cal_sl1 = tmp.data(:, 2);

% SL-3
tmp = importdata(fullfile(cal.calDataPath, 'stellarnet', 'SL3CAL-15081437-ATPLANE.ICD'));
wls_cal_sl3 = tmp.data(:, 1);
spd_cal_sl3 = tmp.data(:, 2);

% Slide projector
fileToLoad = fullfile(cal.calDataPath, 'absoluteirradiance', '20141106T110508_1A_PR-670_Kodak_Carrousel_4400_Projector.mat');
load(fileToLoad, 'meas');
spd_cal_slide = meas.pr670.spd;
wls_cal_slide = SToWls(meas.pr670.S);

%% Plot the spectra
% Normalize at 555 nm
subplot(5, 1, 1)
RGB_sl1 = [255 160 0]/255;
RGB_sl3 = [147 112 219]/255;
RGB_slide = [255 0 0]/255;
plot(cal.wls, spd_USB2000_sl1/max(spd_USB2000_sl1), 'Color', RGB_sl1); hold on;
plot(cal.wls, spd_USB2000_sl3/max(spd_USB2000_sl3), 'Color', RGB_sl3);
plot(cal.wls, spd_USB2000_slide/max(spd_USB2000_slide), 'Color', RGB_slide);
xlabel('Wavelength [nm]');
ylabel('Normalized power');
pbaspect([1 1 1]); set(gca, 'TickDir', 'Out');
xlim([150 1100]);
ylim([0 2]);

allMeas = [spd_USB2000_sl1 spd_USB2000_sl3 spd_USB2000_slide];
% Find out which source has the most power where in the spectrum
theValidWls = ~isnan(spd_USB2000_sl1);
theInvalidWls = isnan(spd_USB2000_sl1);
theValidWlsIdx = find(theValidWls);
[~, theSpec] = max(allMeas, [], 2);

%% Find the transition points
transitionPoints(1) = find(cal.wls == 400);
transitionPoints(2) = find(cal.wls == 760);

weighting_sl1 = zeros(1, length(cal.wls));
weighting_sl3 = zeros(1, length(cal.wls));
weighting_slide = zeros(1, length(cal.wls));

weighting_sl3(1:transitionPoints(1)) = 1;
weighting_slide(transitionPoints(1):transitionPoints(2)) = 1;
weighting_sl1(transitionPoints(2):end) = 1;

weighting_sl1 = smooth(weighting_sl1, 40);
weighting_sl3 = smooth(weighting_sl3, 40);
weighting_slide = smooth(weighting_slide, 40);

wlTransition1 = find((weighting_sl3 == weighting_slide) & (weighting_sl3 > 0));
wlTransition2 = find((weighting_sl1 == weighting_slide) & (weighting_sl1 > 0));

weighting_sl1(theInvalidWls) = NaN;
weighting_sl3(theInvalidWls) = NaN;
weighting_slide(theInvalidWls) = NaN;

subplot(5, 1, 2);
plot(cal.wls, abs(corr_window(spd_USB2000_sl3, 20)), 'Color', RGB_sl3, 'LineWidth', 1); hold on
plot(cal.wls, abs(corr_window(spd_USB2000_slide, 20)), 'Color', RGB_slide, 'LineWidth', 1);
plot(cal.wls, abs(corr_window(spd_USB2000_sl1, 20)), 'Color', RGB_sl1, 'LineWidth', 1);

pbaspect([1 1 1]);
box off;
xlabel('Wavelength [nm]');
ylabel('Abs. corr. coeff.');
pbaspect([1 1 1]); set(gca, 'TickDir', 'Out');
xlim([150 1100]);
ylim([0 1.05]);
title('Correlation');

% Find out the correction factors
% First, we determine the absolute correction factor from the slide
% measurements
spd_cal_sl1 = medfilt1(SplineSpd(wls_cal_sl1, spd_cal_sl1, cal.wls), medianFilterSpan);
spd_cal_sl3 = medfilt1(SplineSpd(wls_cal_sl3, spd_cal_sl3, cal.wls), medianFilterSpan);
spd_cal_slide = medfilt1(SplineSpd(wls_cal_slide, spd_cal_slide, cal.wls), medianFilterSpan);

correct_slide = spd_cal_slide ./ spd_USB2000_slide;
correct_slide(correct_slide <= 0) = NaN;

correct_sl1 = spd_cal_sl1 ./ spd_USB2000_sl1;
correct_sl1(correct_sl1 <= 0) = NaN;

correct_sl3 = spd_cal_sl3 ./ spd_USB2000_sl3;
correct_sl3(correct_sl3 <= 0) = NaN;

subplot(5, 1, 3); % Plot the weighting functions
plot(cal.wls, weighting_sl1, 'Color', RGB_sl1); hold on;
plot(cal.wls, weighting_sl3, 'Color', RGB_sl3); hold on;
plot(cal.wls, weighting_slide, 'Color', RGB_slide); hold on;
ylim([-0.05 1.05]); xlim([150 1100]);
pbaspect([1 1 1]);


%% Find out where we are non-NaN
sl1_valid = (~isnan(correct_sl1) & ~isnan(correct_slide));
sl3_valid = (~isnan(correct_sl3) & ~isnan(correct_slide));

% Bring the correction factors into register by scaling at the transition
% points ±20 nm
subplot(5, 1, 4);
sl3_scalar = correct_sl3(transitionPoints(1)-20:transitionPoints(1)+20) \ correct_slide(transitionPoints(1)-20:transitionPoints(1)+20);
sl1_scalar = correct_sl1(transitionPoints(2)-20:transitionPoints(2)+20) \ correct_slide(transitionPoints(2)-20:transitionPoints(2)+20);
plot(cal.wls, sl3_scalar*correct_sl3, 'Color', RGB_sl3, 'LineWidth', 1.5); hold on;
plot(cal.wls, correct_slide, 'Color',RGB_slide, 'LineWidth', 1.5);
plot(cal.wls, sl1_scalar*correct_sl1, 'Color', RGB_sl1, 'LineWidth', 1.5);
ylim([0 yLimUpper]); xlim([150 1100]);
pbaspect([1 1 1]);
box off;
switch cal.whichMeter
    case 'a'
        set(gca, 'YTick', [0 0.01 0.02 0.03]);
        set(gca, 'YTickLabel', [0 1 2 3]);
    case 'b'
        set(gca, 'YTick', 100*[0 0.01 0.02 0.03]);
        set(gca, 'YTickLabel', [0 1 2 3]);
end
xlabel('Wavelength [nm]');
ylabel('Correction factor');
set(gca, 'TickDir', 'Out');
title('Correction factors');

weighting_sl1(isnan(weighting_sl1)) = 0; correct_sl1(isnan(correct_sl1)) = 0;
weighting_sl3(isnan(weighting_sl3)) = 0; correct_sl3(isnan(correct_sl3)) = 0;
weighting_slide(isnan(weighting_slide)) = 0; correct_slide(isnan(correct_slide))= 0;
% Replace NaN with 0 for a second

corrFac = sl1_scalar*weighting_sl1.*correct_sl1 + sl3_scalar*weighting_sl3.*correct_sl3 + weighting_slide.*correct_slide;
corrFac(corrFac == 0) = NaN;

switch cal.whichMeter
    case 'a'
        wlLow = 280;
        wlHigh = 840;
    case 'b'
        wlLow = 360;
        wlHigh = 840;
end

corrFacOrig = corrFac;
corrFac(1:(find(cal.wls == wlLow)-1)) = NaN;
corrFac((find(cal.wls == wlHigh)+1):end) = NaN;

subplot(5, 1, 5);
plot(cal.wls, corrFacOrig, '-r', 'LineWidth', 2); hold on;
plot(cal.wls, corrFac, '-k', 'LineWidth', 2);
ylim([0 yLimUpper]); xlim([150 1100]);
switch cal.whichMeter
    case 'a'
        set(gca, 'YTick', [0 0.01 0.02 0.03]);
        set(gca, 'YTickLabel', [0 1 2 3]);
    case 'b'
        set(gca, 'YTick', 100*[0 0.01 0.02 0.03]);
        set(gca, 'YTickLabel', [0 1 2 3]);
end

pbaspect([1 1 1]);
box off;
xlabel('Wavelength [nm]');
ylabel('Correction factor');
set(gca, 'TickDir', 'Out');

set(gcf, 'PaperPosition', [0 0 4 16])
set(gcf, 'PaperSize', [4 16]); %Set the paper to have width 5 and height 5.
saveas(gcf, fullfile(cal.plotPath, ['IlluminationSpectraDataset_Calibration_FigureS7_SensitivityCalibration_' cal.whichMeter '.pdf']), 'pdf');
close(gcf);

cal.effectiveWlRange = [wlLow wlHigh];
cal.correctionFactor = corrFac;