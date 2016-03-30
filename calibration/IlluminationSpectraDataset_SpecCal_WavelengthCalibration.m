function wavelengthCorrect = IlluminationSpectraDataset_SpecCal_WavelengthCalibration(cal)
% wavelengthCorrect = IlluminationSpectraDataset_SpecCal_WavelengthCalibration(cal)
%
% This script loads in wavelength measurements of Ar and Hg line sources
% and compares them with what they should. These were done for the two
% spectrometers.
%
% 11/4/2014     spitschan       Wrote it.
% 3/31/2015     spitschan       Updated paths

switch cal.whichMeter
    case 'a'
        %% Ar source measurements
        fileToLoad = fullfile(cal.calDataPath, 'linesources', '20140710_OO_USB2000+a_LineSourceAr.mat');
        load(fileToLoad, 'meas');
        meas_Ar = meas;
        meas_Ar.meanSpd = mean(meas_Ar.spd, 2);
        meas_Ar.meanSpd = (mean(meas_Ar.spd,2 )-mean(meas_Ar.spdDark, 2))/max(mean(meas_Ar.spd,2 )-mean(meas_Ar.spdDark, 2));
        
        %% Hg source measurements
        fileToLoad = fullfile(cal.calDataPath, 'linesources', '20140710_OO_USB2000+a_LineSourceHg.mat');
        load(fileToLoad, 'meas');
        meas_Hg = meas;
        meas_Hg.meanSpd = mean(meas_Hg.spd, 2);
        meas_Hg_spd_norm = (mean(meas_Hg.spd,2 )-mean(meas_Hg.spdDark, 2))/max(mean(meas_Hg.spd,2)-mean(meas_Hg.spdDark, 2));
        
    case 'b'
        %% Ar source measurements
        fileToLoad = fullfile(cal.calDataPath, 'linesources', '20140710_OO_USB2000+b_LineSourceAr.mat');
        load(fileToLoad, 'meas');
        meas_Ar = meas;
        meas_Ar.meanSpd = mean(meas_Ar.spd, 2);
        meas_Ar.meanSpd = (mean(meas_Ar.spd,2 )-mean(meas_Ar.spdDark, 2))/max(mean(meas_Ar.spd,2 )-mean(meas_Ar.spdDark, 2));
        
        %% Hg source measurements
        fileToLoad = fullfile(cal.calDataPath, 'linesources', '20140710_OO_USB2000+b_LineSourceHg.mat');
        load(fileToLoad, 'meas');
        meas_Hg = meas;
        meas_Hg.meanSpd = mean(meas_Hg.spd, 2);
        meas_Hg_spd_norm = (mean(meas_Hg.spd,2 )-mean(meas_Hg.spdDark, 2))/max(mean(meas_Hg.spd,2 )-mean(meas_Hg.spdDark, 2));
        
end

%% Ar source measurements

% Plot it
figure(1); clf; hold on
plot(meas_Ar.wls,meas_Ar.meanSpd,'b');
title('AR measurements');
drawnow;

% Now find the measured peak in each of the candidate
% peak locations
calibNominal = [];
calibMeasured = [];

arcandidatePeaks = [696.5 706.7 727.3 738.4 763.5];
wlSpread = 4;
for i = 1:length(arcandidatePeaks)
    theIndex = find(meas_Ar.wls > arcandidatePeaks(i)-wlSpread & ...
        meas_Ar.wls < arcandidatePeaks(i)+wlSpread);
    theWls = meas_Ar.wls(theIndex);
    theData = meas_Ar.meanSpd(theIndex);
    [~,peakIndex] = max(theData);
    peakWl = theWls(peakIndex);
    fprintf('Candidate wavelength %g, peak at %g\n',...
        arcandidatePeaks(i),peakWl);
    
    % Draw the peaks onto the spectrum
    figure(1); hold on
    plot([peakWl peakWl],[0 max(meas_Ar.meanSpd)],'g');
    plot([arcandidatePeaks(i) arcandidatePeaks(i)],[0 max(meas_Ar.meanSpd)],'r');
    
    hold off
    drawnow;
    
    % Make a blowup
    figure(i+1); clf;
    plot(theWls,theData,'b');
    hold on
    plot(theWls,theData,'b+');
    plot([peakWl peakWl],[0 max(theData)],'g');
    plot([arcandidatePeaks(i) arcandidatePeaks(i)],[0 max(theData)],'r');
    hold off
    xlabel('Wavelength (nm)');
    ylabel('Power');
    title(sprintf('AR nominal peak at %g nm',arcandidatePeaks(i)));
    drawnow;
    
    calibNominal = [calibNominal arcandidatePeaks(i)];
    calibMeasured = [calibMeasured peakWl];
end

%% Hg source measurements
% Plot it
figure(2+length(arcandidatePeaks)); clf; hold on
plot(meas_Hg.wls,meas_Hg.meanSpd,'b');
title('HG measurements');
drawnow;

% Now find the measured peak in each of the candidate
% peak locations
hgcandidatePeaks = [404.7 435.8 546.1 579];
wlSpread = 4;
for i = 1:length(hgcandidatePeaks)
    theIndex = find(meas_Hg.wls > hgcandidatePeaks(i)-wlSpread & ...
        meas_Hg.wls < hgcandidatePeaks(i)+wlSpread);
    theWls = meas_Hg.wls(theIndex);
    theData = meas_Hg.meanSpd(theIndex);
    [~,peakIndex] = max(theData);
    peakWl = theWls(peakIndex);
    fprintf('Candidate wavelength %g, peak at %g\n',...
        hgcandidatePeaks(i),peakWl);
    
    % Draw the found peaks onto the spectrum
    figure(2+length(arcandidatePeaks)); hold on
    plot([peakWl peakWl],[0 max(meas_Hg.meanSpd)],'g');
    plot([hgcandidatePeaks(i) hgcandidatePeaks(i)],[0 max(meas_Hg.meanSpd)],'r');
    hold off
    drawnow;
    
    % Make a blowup
    figure(i+ 2 + length(arcandidatePeaks)); clf;
    plot(theWls,theData,'b');
    hold on
    plot(theWls,theData,'b+');
    plot([peakWl peakWl],[0 max(theData)],'g');
    plot([hgcandidatePeaks(i) hgcandidatePeaks(i)],[0 max(theData)],'r');
    hold off
    xlabel('Wavelength (nm)');
    ylabel('Power');
    title(sprintf('HG nominal peak at %g nm',hgcandidatePeaks(i)));
    drawnow;
    
    calibNominal = [calibNominal hgcandidatePeaks(i)];
    calibMeasured = [calibMeasured peakWl];
end

figure; clf; hold on
wavelengthCorrect = mean(calibMeasured)-mean(calibNominal);
fprintf('Nominal wavelength correct (subtract from nominal wl to get actual wl) is %0.3g nm\n',wavelengthCorrect);
plot([150 1100],[150 1100],'-k','LineWidth',0.5);
plot([150 1100],[150 1100]+wavelengthCorrect,'--k','LineWidth',0.5, 'Color', 'r');
h1 = plot(calibNominal(1:5),calibMeasured(1:5),'sr','MarkerFaceColor','r');
h2 = plot(calibNominal(6:end),calibMeasured(6:end),'sb','MarkerFaceColor','b');
axis('square'); axis([150 1100 150 1100]);
set(gca, 'XTick', [150:100:1100]);
set(gca, 'YTick', [150:100:11001]);
set(gca, 'TickDir','Out')
xlabel('Nominal peak wavelength [nm]'); ylabel('Measured peak wavelength [nm]');
legend([h1 h2], 'Ar source', 'Hg source', 'Location', 'NorthWest'); legend boxoff;

set(gcf, 'PaperPosition', [0 0 4 4])
set(gcf, 'PaperSize', [4 4]); %Set the paper to have width 5 and height 5.
saveas(gcf, fullfile(cal.plotPath, ['IlluminationSpectraDataset_Calibration_FigureS6_WlCorrection_' cal.whichMeter '.pdf']), 'pdf');
close all;
%%
figure; clf; hold on
[~, orderWl] = sort(calibNominal);
calibNominalOrdered = calibNominal(orderWl);
calibMeasuredOrdered = calibMeasured(orderWl);
diffWls = calibMeasuredOrdered-calibNominalOrdered;
diffWlsCounter = 1:length(diffWls);
for i = 1:length([1:diffWlsCounter(end)+1])
    plot([-2 2], [i i], 'Color', [0.75 0.75 0.75], 'LineWidth', 0.2);
end
for i = 1:length(diffWls)
    if find(hgcandidatePeaks == calibNominalOrdered(i))
        plot(diffWls(i), diffWlsCounter(i), 'sb', 'MarkerFaceColor', 'b')
    elseif find(arcandidatePeaks == calibNominalOrdered(i))
        plot(diffWls(i), diffWlsCounter(i), 'sr', 'MarkerFaceColor', 'r')
    end
end
plot(mean(diffWls), length(diffWls)+1, 'sk', 'MarkerFaceColor', 'k');
errorbarX(mean(diffWls), length(diffWls)+1, std(diffWls), '-k');
set(gca, 'YTick', [diffWlsCounter diffWlsCounter(end)+1])
for i = 1:length(calibNominalOrdered)
    tickLabels{i} = [num2str(calibNominalOrdered(i)) ' nm'];
end
a = text(mean(diffWls), diffWlsCounter(end)+2, [num2str(mean(diffWls), '%.2f') '\pm' num2str(std(diffWls), '%.2f') ' nm']);
set(a, 'HorizontalAlignment', 'center');
plot([0 0], [0 diffWlsCounter(end)+2], '-k');
tickLabels{end+1} = 'Mean\pm1SD';
set(gca, 'YTickLabel', tickLabels);
set(gca, 'XTick', [-2 -1 0 1 2]);
xlim([-2 2]);
ylim([0 length(diffWls)+2]);
xlabel('Wavelength shift [nm]');
set(gca, 'TickDir','Out')
pbaspect([0.6 1 1]);
set(gcf, 'PaperPosition', [0 0 2.4 4])
set(gcf, 'PaperSize', [2.4 4]); %Set the paper to have width 5 and height 5.
saveas(gcf, fullfile(cal.plotPath, ['IlluminationSpectraDataset_Calibration_FigureS6_WlCorrectionInIlluminationSpectraDataset_' cal.whichMeter '.pdf']), 'pdf');
close all;